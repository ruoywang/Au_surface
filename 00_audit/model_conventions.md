# 00_audit/model_conventions.md

本文件回答 v2 文档第 1.1–1.3 节要求核对的问题。每条结论后标注来源文件:行号；标"[待定]"的是读源码仍不能唯一确定、需要小规模试算或进一步深挖才能回答的问题，不得当作已验证结论使用。

## 1.1 恒电势控制（LCEP）

**电子数由唯一一个控制器更新**：`src/main.F` 中挂在 `LCEP` 分支下的电子数调节逻辑（3196–3298行）。VASPsol++ 自带的 `update_NELECT`/`check_EFERMI`（`solvation.F:3623` 起）是否在别处也被调用、会不会和 `LCEP` 双重触发——搜索确认这两个函数在 `src/*.F` 中**未被任何地方调用**（`PUBLIC` 导出但无调用点），本地实现里唯一生效的电子数控制器就是 `main.F` 的 LCEP 块。[已确认，非待定]

**核心发现（解决了此前提出的关键疑点）**：电子数收敛到 `TARGETMU` 的迭代，对 `NESCHEME∈{1,2,3}` 是一个**内层 `DO...EXIT` 循环**（`main.F:3211-3243`），在同一个离子步（`NSTEP`）内反复：更新 `NELECT` → 重跑一次完整 SCF（`CALL ELECTRONIC_OPTIMIZATION`）→ 检查 `ABS(TARGETMU-EFERMI) < FERMICONVERGE`，不满足就继续。这个内层循环**与 `DYN%IBRION`/`NSW` 无关**，对 `IBRION=-1, NSW=0` 的静态单点计算同样会跑满多轮 SCF 直到收敛。因此：

- v2 模板用 `NESCHEME=3` 做静态标签是源码确认可行的，符合"逐构型收敛到目标电子化学势"的要求。
- 旧算例用的 `NESCHEME=5`（Nose-Hoover 链，`main.F:3260-3279`）和 `NESCHEME=4`（velocity-Verlet 型，`main.F:3247-3259`）只在每个**外层离子步**更新一次 `NELECT`，用 `DYN%POTIM*NEADJUST` 做时间积分——这是允许费米能围绕目标值涨落的 MD 方案，不能用于静态单点标签，**不要把旧例子的 `NESCHEME=5` 带进新工作流**，这一点 v2 文档已经写对了，本次核实只是确认它有源码依据。

外层 "ion" 循环的退出条件（`main.F:3121-3127`）：只有当 `DYN%IBRION∈{1,2,3}` 时才会附加"最后一次 SCF 的费米能已收敛到 `TARGETMU`"这个门槛；`IBRION=-1` 时纯粹按 `INFO%LSTOP` 退出。这对纯静态单点没有影响（内层循环已经保证收敛）。但**如果后续对代表性缺陷做正负偏压弛豫（v2 §4.1，`IBRION=2`）**，就会落入这个附加门槛：每一步弛豫都必须先让电子数收敛到 `FERMICONVERGE` 才能继续判断力是否收敛，这会显著增加每个弛豫步的 SCF 次数，做预算估计时要计入。

其它确认：
- `INIT_ECHANGE`、`ECHANGE_MAX` 也是通过 `PROCESS_INCAR` 从 INCAR 读取的独立标签（`src/reader.F`，与 `TARGETMU` 等一起声明），v2 文档 §4.1"保存全部解析后的默认参数"清单里目前没提到，需要补上。
- 每次更新电子数时，电荷密度用 `CHTOT = CHTOT * (NELECT_new/NELECT_old)` 做整体等比缩放作为下一次 SCF 的初猜（`main.F:3212-3213`），不是重新分布到具体能带——这是一个简单缩放而非变分最优的重启方式，记录下来供后续解释 SCF 收敛速度用。
- Grand potential 诊断量 `GCE = (TOTEN - EENTROPY/(2+NORDER)) + (CNNE-old_nele)*TARGETMU`（`main.F:3787/3954/4170`，写入 OSZICAR 风格的 17 号单元）。这与 v2 §9.3 的 `Omega_N`/`Omega_q` 定义是否严格对应、`CNNE` 与 `old_nele` 的具体差异含义，**[待定]**——留给 §5/§9.3 的能量自洽核对阶段用小体系有限差分直接验证，这里不下结论。

## 1.2 电子—空腔—偶极—溶剂闭环

在 `src/pot.F`（VASP 总势组装例程）中确认了完整链条：

1. `CVHAR := POTHAR(CHTOT)`（纯电子 Hartree 势，`pot.F:363-367`）
2. `CVHAR := CVHAR + POTION(...)`（叠加离子局域赝势，`pot.F:376, 405-408`）→ **此时的 `CVHAR` 就是 v2 文档里的 `phi_explicit = POTHAR+POTION`，源码确认，未混入 XC 或 DENCOR。**
3. 若 `DIP%LCOR_DIP`（即 `LDIPOL`）为真且 `ICALL>1`（第一次调用时还没有可用的电荷/偶极矩，跳过），直接调用官方 `CDIPOL`（`dipol.F`）**就地修改 `CVTOT`**（`pot.F:335-336`）。`CDIPOL` 内部读取 `solvation_moments` 模块的 `Qsol_cache`/`Dsol_cache`（由 `CALC_SOLVENT_MOMENTS(n_b+n_ion, ...)` 算出，`solvation.F:1244`），把溶剂的净电荷/偶极矩并入官方多极矩计算，再喂给 `EWALD_DIPOL` 生成锯齿势场（`dipol.F` ~671-712）。
4. 紧接着把这次 `CDIPOL` 造成的 `CVTOT` 增量单独抽出来存成 `CVDIP := CVTOT_after − CVTOT_before`（`pot.F:332-338`），**目的仅仅是把这个已经生效的场同步告诉溶剂/离子求解器**（传入 `SOL_Vcorrection`，`pot.F:419`），并不会再把 `CVDIP` 加回 `CVTOT` 第二次。

**结论：本地实现里只有一条偶极修正代码路径——官方 `LDIPOL/IDIPOL/DIPOL` 机制本身已经被扩展为"溶剂电荷感知"的版本。不存在需要叠加/去重的第二套修正；`CVDIP` 只是一份记账副本，不是独立的物理修正项。** v2 §1.2/§2.3 提醒的"不要另加第二个 slab dipole correction"，落到这份源码上就是：**只要打开 `LDIPOL=.TRUE.`，闭环就已经包含溶剂偶极反馈，不需要、也不存在另一个开关去做"真正的" CEP-DIP 偶极修正。**

5. `SOL_Vcorrection`（`solvation.F:3435`）按 `ISOL` 分派到 `LPCM_solver`（`ISOL=1`，线性，`solvation.F:1174`）或 `NLPCM_solver`（`ISOL=2`，非线性，对应 [R2] 文献模型，`solvation.F:1733`）。两个分支都求解得到：
   - `phi`（= `phi_sol + phi_solv`，即 `phi_explicit + phi_solv`，也就是 v2 文档的 `phi_used_by_ions`）
   - `V_corr`：真正加回 `CVTOT` 的修正项（`pot.F:424-426: CVTOT = CVTOT + CVHAR + V_corr`）
   - `n_b`（RHOB，介电束缚电荷）、`n_ion`（RHOION，移动离子净电荷）

**重要的、v2 文档没写到的细节**：`ISOL` 决定了输出文件 `PHI` 的实际含义——
- `ISOL=2`（`NLPCM_solver`，本项目和旧算例都用这个）：`PHI` 文件写的是 `phi`（**总电势**，`phi_explicit+phi_solv`，`solvation.F:1846`）。
- `ISOL=1`（`LPCM_solver`）：`PHI` 文件写的却是 `phi_solv`（`solvation.F:1251`），跟 `PHISOLV` 文件内容完全相同。

即**同一个文件名 `PHI` 在两种 `ISOL` 下含义不同**，这正是 v2 文档反复强调的"数组名不能望文生义，需核对源码"的一个具体例子。**本项目必须显式写 `ISOL=2`**（`reader.F:3407` 显示 VASP 内部默认 `ISOL=1`，不写就会静默用错误分支）。

`phi_explicit` 本身（纯 `POTHAR+POTION+官方偶极场`，溶剂反馈之前）**在当前版本里没有单独的输出文件**——磁盘上只有 `PHI`(=phi_explicit+phi_solv，仅 ISOL=2 时)、`PHISOLV`(=phi_solv)、`VSOLV`(=V_corr)。如果后续分析（v2 §6/§7）确实需要单独的 `phi_explicit`，优先在后处理里用 `PHI − PHISOLV` 重建（零风险），不必先改源码新增输出。

`ELOC`、`P` 两个输出文件分别是局域电场、极化强度**仅 z 分量**（`E_loc(:,3)`、`P(:,3)`，`solvation.F:1849-1850`），不是三维矢量场。v2 §7.3 如果需要横向（x/y）场分量，目前的只读输出补丁还没有覆盖，需要作为新增只读输出另外加（按 v2 §0 的规则：新增只读输出走独立分支+回归测试，不改变解）。

空腔/掩码输出核对：`SVDW`、`SION`、`SSOLV`、`SCAV`、`SDIEL`（`solvation.F:2048/2067/2085/2104/2124`，均为 `ISOL=2` 分支）。v2 §6 只列了 `SVDW, SDIEL, SION` 三个，实际还有 `SSOLV`、`SCAV` 两个额外掩码，补充进字段清单。

## 1.3 单侧标签映射

- `LVAC`（逻辑）、`SOL_Z0`、`SOL_Z1`、`SOL_SIGMA`（`solvation.F:1160-1161`，INCAR 读取见 1503-1541行）：erfc 型平滑窗口，作用在 `S_diel`（介电掩码，通过 `M_SOL` 相乘，`solvation.F:2123`）。`SOL_SIGMA` 默认 0.8 Å。
- **离子窗口是单独派生的、更窄的一个窗口**，不是与介电窗口共用：`ION_Z0 = SOL_Z0 + D_STERN`，`ION_Z1 = SOL_Z1 - D_STERN`，`ION_SIGMA = SOL_SIGMA`（`solvation.F:2008-2010`），`D_STERN` 默认 2.0 Å（`solvation.F:1722`）。`M_ION` 相乘作用在 `S_ion`(SION)、`S_cav`、`SWORK`(SVDW 相关) 上（`solvation.F:2066/2084/2103`），`M_SOL` 只作用在 `S_diel`(SDIEL) 上（`solvation.F:2123`）。**这是 v2 文档目前没写清楚的一点：单侧窗口对介电和对离子实际是两个不同宽度的窗口，相差两倍 `D_STERN`（Stern 层厚度），不能当成同一个窗口处理。**
- 位置单位：`SOL_Z0`/`SOL_Z1` 用 `RDATAB(...,'F',...)` 按普通 `REAL` 读入，源码注释明确写"Angstrom, Cartesian along slab normal"（`solvation.F:1514`）——**确认是笛卡尔 Å，不是分数坐标**。
- `DIPOL`（即 `POSCEN`）**确认是分数/直接坐标**：读入后与 `T_INFO%POSION`（VASP 内部分数坐标）用 `MOD(POSION-POSCEN+10.5,1.0)-0.5` 做周期wrap（`dipol.F:118-120,152,360-362`）。**`DIPOL` 只是计算多极矩的参考中心，不是修正势跳跃面的位置**——跳跃面的位置由官方 `CDIPOL`/`EWALD_DIPOL` 的锯齿场构造方式决定，固定在沿 `IDIPOL` 方向的周期晶胞边界上，与 `DIPOL` 取值无关（`dipol.F` ~930-1047 的锯齿场构造未见依赖 `POSCEN` 之外的位置输入）。**这直接回答 v2 §2.3 的核心问题：移动 `DIPOL` 不会把跳跃面移进/移出真空，只有改变真空缓冲厚度或晶胞大小才能做到。** 单侧液层很长时的风险仍然存在——不是因为 `DIPOL` 选错了会把跳跃面带进液区，而是因为如果真空缓冲本身不够厚（晶胞沿该方向的周期边界离活跃液区太近），跳跃面（在周期边界）自己就会落进液区。
- v2 文档提到的 `VAC_Z0/VAC_Z1` **在源码中不存在**，确认是占位猜测，不是真实标签；单侧窗口只有 `SOL_Z0/SOL_Z1`（及派生的 `ION_Z0/ION_Z1`）。已在 v2 文档更新清单里去掉这个假设。

## 待后续核实的开放项（不是已验证结论）

1. `GCE`/`TOTEN`/`CNNE` 与 v2 §9.3 `Omega_N`/`Omega_q` 定义的严格对应关系——留给能量-力一致性有限差分阶段核实。
2. 编译时实际使用的 `arch/makefile.include.*` 具体是哪一份（现有多个 intel 变体），不影响功能正确性，只影响可复现性记录，非阻塞项。
3. `phi_explicit` 用 `PHI-PHISOLV` 重建是否在数值上精确等于源码内部真正的 `phi_sol`（理论上应该精确相等，因为都是同一浮点数组的线性运算，但尚未跑一次算例做数值验证）。
4. 正负偏压弛豫（`IBRION=2`）下 CP 内层循环与几何优化外层循环交替收敛的实际迭代次数/成本，需要 pilot 小算例实测，不能只靠读代码估计。
