# 00_audit/parameter_map.md

按补充说明第3节要求的映射表整理；详细依据见 `model_conventions.md`，本文件只做速查。全部来自 `/anvil/projects/x-che190065/rywang/CEP-HALF` 源码核对，标注文件:行号。

## A. 核对内容 → 当前机器结论

| 要核对的内容 | 当前机器结论 |
|---|---|
| 哪个模块更新电子数 | 唯一控制器：`main.F` 的 `LCEP` 分支（3196-3298行）。VASPsol++ 自带的 `update_NELECT`/`check_EFERMI`（`solvation.F:3623`起）未被调用，不存在双重控制器冲突。 |
| 目标化学势如何定义/收敛 | `NESCHEME∈{1,2,3}`：内层 `DO...EXIT` 循环，同一离子步内反复调 SCF 直到 `\|TARGETMU-EFERMI\|<FERMICONVERGE`（`main.F:3211-3243`），与 `IBRION/NSW` 无关，静态单点可用。`NESCHEME=4/5`：每个外层离子步只更新一次，按 MD 时间积分，费米能允许涨落，只适合 MD（`main.F:3247-3279`）。 |
| 单侧窗口标签、单位、位置、作用对象 | `LVAC/SOL_Z0/SOL_Z1/SOL_SIGMA`（Å，笛卡尔，沿 slab 法向；作用于 `S_diel`）+ 派生 `ION_Z0=SOL_Z0+D_STERN, ION_Z1=SOL_Z1-D_STERN`（作用于 `S_ion/S_cav/SWORK`），`D_STERN` 默认 2.0 Å。两个窗口不同宽，相差 2×D_STERN。（`solvation.F:1160-1161,1503-1541,1717-1722,2008-2124`） |
| 偶极修正：参考中心、修正位置、溶剂是否参与反馈 | `DIPOL`=分数坐标，仅为多极矩参考中心，不是跳跃面位置；跳跃面固定在周期晶胞边界（与 `DIPOL` 无关）。溶剂电荷/偶极矩（`Qsol_cache/Dsol_cache`，来自 RHOB+RHOION）已并入官方 `CDIPOL`（`dipol.F`），只有一条修正路径，`CVDIP` 只是记账副本不重复叠加。（`pot.F:323-338,415-426`；`dipol.F` 多处） |
| 输出场：哪个数组真正驱动离子响应的电势 | `phi`(=phi_explicit+phi_solv) 是驱动 PB/离子求解的势，`ISOL=2` 时写入文件 `PHI`；`ISOL=1` 时文件 `PHI` 实际写的是 `phi_solv`（与 `PHISOLV` 相同）——**必须显式设 `ISOL=2`**，默认值是 1。真正加回 KS 哈密顿量的是 `V_corr`（文件 `VSOLV`），不是 `PHI`。（`solvation.F:1174-1877,3435-3578`；`reader.F:3407`） |
| 输出能量：是否已含电子储库/偶极修正 | `TOTEN` 本身未见含 `(N-N0)*mu` 项；`GCE=(TOTEN-EENTROPY/(2+NORDER))+(CNNE-old_nele)*TARGETMU` 是另外单独算出、写到 17 号单元（OSZICAR 风格）的诊断量。与 v2 §9.3 `Omega_N/Omega_q` 的精确对应关系**[待定]**，需小体系数值核对。（`main.F:3786-3790,3953-3965,4169-4176`） |
| 执行文件与源码对应关系 | `bin/vasp_std`（2025-10-21 构建，`sha256=5e7152b7...a702c41`）对应本次读取的 `src/*.F`（`main.F/reader.F/solvation.F/dipol.F/pot.F`）；`makefile.include` 具体版本未逐一比对，非阻塞项。 |

## B. INCAR 标签速查（本地源码确认存在，含默认值/单位/文件位置）

| 标签 | 默认值 | 单位/类型 | 作用 | 源码位置 |
|---|---|---|---|---|
| `LCEP` | `.FALSE.` | 逻辑 | 开启恒电势 | reader.F:862 |
| `TARGETMU` | `-4.6` | eV | 目标电子化学势 | reader.F:847 |
| `FERMICONVERGE` | `0.05` | eV | 费米能收敛阈值 | reader.F:859 |
| `CAP_MAX` | `2.0` | 电子/eV | 电子数-电势响应"电容"上限 | reader.F:856 |
| `NESCHEME` | `1` | 整数(1-5) | 电子数调节方案，见上表 | reader.F:868 |
| `NEADJUST` | `1` | 整数(离子步) | 每几步调整一次电子数 | reader.F:865 |
| `INIT_ECHANGE` | 需查INCAR初值 | 电子 | 第一步电子数变化步长 | reader.F(与TARGETMU同段声明) |
| `ECHANGE_MAX` | 需查INCAR初值 | 电子 | 单次电子数变化上限 | reader.F(同上) |
| `LSOLOFF` | `.FALSE.`(IBRION=0时自动True) | 逻辑 | 关闭溶剂响应 | reader.F:871-873 |
| `T_eta`/`eta_length`/`M_eta` | `TEBEG`/1/... | K/整数 | NESCHEME=5的Nose-Hoover链参数，静态标签不用 | reader.F:881-888 |
| `LSOL` | — | 逻辑 | 开启隐式溶剂 | solvation.F |
| `ISOL` | `1` | 整数(1/2) | 1=线性PCM，2=非线性（本项目必须显式设为2） | solvation.F:3407 |
| `C_MOLAR` | — | mol/L | 体相盐浓度 | solvation.F |
| `R_ION` | — | Å | 离子半径 | solvation.F:2957 |
| `R_SOLV`/`R_DIEL` | — | Å | 溶剂/介电有效半径 | solvation.F:2924,2946 |
| `NC_K`/`SIGMA_K`/`A_K` | — | 1/Å³ / Å / Å | 空腔截断密度/宽度/平滑长度 | solvation.F:1451-1468,2891-2919 |
| `LVAC` | `.FALSE.` | 逻辑 | 开启单侧窗口 | solvation.F:1503-1511 |
| `SOL_Z0`/`SOL_Z1` | `0/0` | Å(笛卡尔) | 介电/整体窗口边界 | solvation.F:1514-1530 |
| `SOL_SIGMA` | `0.8` | Å | 窗口平滑宽度 | solvation.F:1535 |
| `D_STERN` | `2.0` | Å | 离子窗口相对介电窗口的内缩量 | solvation.F:1722 |
| `LDIPOL` | — | 逻辑 | 开启偶极修正（含溶剂反馈，唯一路径） | dipol.F:124 |
| `IDIPOL` | — | 整数(0-4) | 修正方向 | dipol.F:141 |
| `DIPOL` | `(-100,-100,-100)` | 分数坐标 | 多极矩参考中心（非跳跃面位置） | dipol.F:118-120,152 |

## C. 输出场文件速查（`ISOL=2` / `NLPCM_solver` 分支，即本项目实际使用的分支）

| 文件名 | 实际内容 | 源码位置 |
|---|---|---|
| `PHI` | `phi = phi_explicit + phi_solv`（总电势，驱动离子响应） | solvation.F:1846 |
| `PHISOLV` | `phi_solv`（溶剂反馈电势，单独项） | solvation.F:1848 |
| `VSOLV` | `V_corr`（真正加回 KS 哈密顿量/CVTOT 的修正势） | solvation.F:1832 |
| `RHOB` | `n_b`，介电束缚电荷 | solvation.F:1851 |
| `RHOION` | `n_ion`，移动离子净电荷 | solvation.F:1852 |
| `ELOC` | 局域电场 **仅 z 分量** | solvation.F:1849 |
| `P` | 极化强度 **仅 z 分量** | solvation.F:1850 |
| `SVDW`/`SION`/`SSOLV`/`SCAV`/`SDIEL` | 五个空腔/可达性掩码（v2 §6 只列了前三个，实际还有 SSOLV/SCAV） | solvation.F:2048/2067/2085/2104/2124 |
| `phi_explicit`（单独） | **当前版本无此文件**，需要时用 `PHI-PHISOLV` 后处理重建 | — |

## D. 需要更新回 v2 文档正文的地方（清单）

1. §1.1：补充"NESCHEME=3 的收敛是内层 SCF 循环，与 IBRION 无关"这一确认结论；补充 `INIT_ECHANGE`/`ECHANGE_MAX` 到默认参数清单。
2. §1.2：明确"只有一条偶极修正路径，`CVDIP` 是记账副本非独立修正"，去掉"不要叠加第二个修正"里隐含的"可能存在两条路径"的不确定语气。
3. §1.3/§2.3：改正 `VAC_Z0/VAC_Z1` 不存在；补充离子窗口相对介电窗口内缩 `D_STERN` 的机制；明确 `DIPOL` 不是跳跃面位置、跳跃面固定在周期边界。
4. §4.1 模板：必须显式加 `ISOL=2`（否则默认值1会让 `PHI` 文件语义完全不同）；`NESCHEME=3` 已确认可用，但要注明若后续做恒电势弛豫（IBRION=2）需要把内层SCF收敛计入预算。
5. §6 输出字段清单：补充 `SSOLV`、`SCAV` 两个掩码；注明 `phi_explicit` 需后处理重建；`ELOC`/`P` 目前只有 z 分量，如需横向场需新增只读输出。
6. §2.2 液区/真空厚度："40Å/20Å"起点无任何推导依据，不是显式水类比也不是德拜长度推算，应改为按1M德拜长度（λ_D≈3.0Å）重新给出液区20Å/真空10Å量级的起点，并保留必须靠§5平台判据实测确认的措辞。
7. §5 验收判据：补充"选平台窗口时要分别对每个场（SDIEL/SION/PHI）核实,不能用同一个z范围套所有场"——SION 的物理空腔开放位置系统性地比 SDIEL 靠外（约一个 R_ION+cavity smoothing 的距离），用同一个窗口去算统计量会把 SION 的近表面过渡尾算进"平台"，得到虚高的标准差,不是真实的不收敛。

## E. 2026-09-22 最小验证/首个真实PZC计算的操作性发现（并行、内存、MPI）

这些是跑通第一个真正收敛的中性T参考点过程中踩出来的坑，跟审计本身的源码结论不同，是**操作层面**的经验，production阶段会反复用到，一并记在这里：

1. **`KPAR>1` 在这套代码里不划算**：`KPAR` 会让每个k点组各自复制一份完整的三维网格（包括溶剂/空腔求解用的额外大数组），实测 `KPAR=4`（32核）连续两次卡在~59-62GB附近被OOM杀掉；改 `KPAR=2` 内存问题缓解但每步耗时没有变化（跟 `KPAR=1` 基本一样，8-9分钟/步），怀疑是`NLPCM`非线性PB求解器没有对KPAR做感知，各k点组各自重复求解同一个（与k点无关的）电荷/空腔问题，抵消了k点并行本该带来的收益。**结论：这套代码目前不要用 `KPAR>1`，把核数都留给单个k点组内的band/FFT并行（`NCORE`）。**
2. **真正决定速度的是总核数，不是 KPAR/NCORE 怎么分配**：32核时无论 `KPAR/NCORE` 怎么调都是~8分钟/DAV步；换成128核（`NCORE=8, NPAR=16, KPAR=1`）稳定在~2分钟/步，跟核数提升到4倍基本对应。用户已有生产MD算例（`surface_charge/6-Au/1-32_water_GCE`）用的是5节点×128核=640核（`NPAR=16/NCORE=8/KPAR=5`），比我们这次128核又快5倍量级，符合线性预期。**做真实产出计算时核数不能按"排队快"随便砍，要按这个体系规模留够总核数。**
3. **Intel MPI 2019.5 在单节点高核数（≥48核）时有已知 SHM collective bug**：`Assertion failed in ch4_shm_coll.c ... node_info->numa_num <= (...)`，在AMD多NUMA域节点上、单节点核数较高时（我们在48核和128核都触发过，16-32核没触发）必现崩溃，VASP还没打印任何东西就在 `MPI_Barrier` 层面段错误。**规避方法：`export I_MPI_COLL_INTRANODE=pt2pt`（强制节点内collective走点对点而不是共享内存），加了之后128核/单节点可以正常跑。这个环境变量应该固化进所有生产 `job-run` 脚本，不止是这次的小测试。**
4. ~~**`ALGO=Fast` 在这套代码里比 `ALGO=Normal` 明显更慢**（实测慢约6倍）~~ **【2026-09-26 撤回，见K节】**：当年的测试同时改了ALGO和KPAR、且只跑了4步Davidson（未进入RMM-DIIS阶段）；重测同体系每步Fast 53~55 s vs Normal 60~67 s，Fast反而快~10%。"溶剂PB求解被触发更多"的解释也错——POTLOK只占每步6%。第1条KPAR结论方向正确但原因需修正（溶剂网格工作按KPAR组重复计算，非OOM本身），第2条"总核数决定速度"仍成立。
5. **Anvil 分区实测差异**（都是同款128核/257GB节点，区别只是调度策略）：
   - `wholenode`/`standard`：要求整节点空闲，繁忙时段排队可能是十几小时（QOS `part-standard` 允许多节点，node上限16）。
   - `shared`：允许部分节点分配，但内存**严格按核数配额**（`MaxMemPerCPU=1896MB`），QOS `part-shared` 单作业上限 `node=1`；集群整体排队深度会让"部分空闲"优势打折扣。
   - `debug`：节点池很小（17节点），QOS `part-debug` 单作业上限 `node=2`、墙钟硬顶2小时；空的时候排队几乎瞬间，忙的时候不比其他分区快，是个易变的"赌一把"选项。
   - `highmem`：**内存按核数配额是 `7.95GB/核`**（比 shared/debug 的1.9GB/核高4倍），QOS `part-highmem` 单作业上限 `node=1`，**计费权重 `CPU=4.0`**（4倍SU消耗）。这次能在128核跑起来主要靠它内存给得多，不需要跟KPAR内存复制的问题死磕。**代价是贵，别当成默认分区，只在需要大内存单节点时用。**
   - 所有分区都有**每用户同时排队上限**（`MaxSubmitPU`，`part-highmem`=4），批量提交时要检查排队数量，不能假设可以无限堆。
6. **第一个真正收敛的中性T（q_e=0）算例**（job 20854473，128核highmem+`I_MPI_COLL_INTRANODE=pt2pt`+`ALGO=Normal`，2:35:24完成，70个DAV步，`dE`最终 6.7E-8<EDIFF）：`E-fermi=-4.9411 eV`；`PHI`场的体相液相平台（z≈23-33Å）极平（mean≈0, std=3.2e-5 eV，远好于<1mV判据）；`SDIEL`/`SION` 平台（各自正确窗口内）都精确等于1.000000；`RHOB+RHOION`积分闭合到-3.3e-6 e。**`mu0 = EFERMI - phi_bulk ≈ -4.941 eV`**，是本项目第一个通过全部§5验收判据的真实数值，不是占位值。真空区存在约1e-4 V/Å量级的微小残余斜率（卡在v2 §7.3给的检查尺度边缘），记录但未判定为拦截性问题。

## F. 补充：highmem 并发上限

`part-highmem` QOS 除了 `MaxSubmitPU=4`（排队总数上限）还有 **`MaxJobsPU=2`**（同时真正在跑的作业数上限）——两者独立，即使节点大量空闲（实测9个idle），本账号在 `highmem` 上也只能2个作业同时RUNNING，其余排队等位。批量提交时按此吞吐量估算总耗时（9点约需 ceil(9/2)×单点耗时）。

## G. 2026-09-22 项目决定：展宽方案改为 ISMEAR=1/SIGMA=0.2

用户决定：以后统一用 `ISMEAR=1, SIGMA=0.2`（Methfessel-Paxton，跟旧MD算例一致），不再用 v2 文档模板默认的 `ISMEAR=-1, SIGMA=0.05`（Fermi-Dirac）。

背景：排查速度差距（我们T结构~130秒/DAV步 vs 旧MD算例4-8秒/离子步）时发现两者展宽方案不同；MP展宽在费米面附近能带分辨要求更低，Davidson/RMM-DIIS每次迭代通常更便宜。这是用户明确的取舍决定，接受的代价：v2文档§4.1原本要求Fermi-Dirac是为了熵项/自由能定义有明确物理温度含义，改用MP展宽后，`TOTEN`/`EENTROPY`等的物理解释和后续§9.3能量-力一致性检查要按新展宽方案重新核对，不能直接套用文档里假设Fermi-Dirac的表述。

## H. 2026-09-22 项目决定：标准配置改为 4层 + 3×3×1 k点

排查速度问题的过程（KPAR/ALGO/ISMEAR都测过，均非主因；6层96原子T结构在128核highmem上稳定~2-2.3分钟/DAV步，跟旧MD算例对比有约13-14倍标准缩放律解释不了的差距，原因未完全查清）中，最终发现**层数从6降到4 + k点从4×4×1降到3×3×1**给出了干净、可用缩放律解释的~4-5倍提速（4L-3×3×1完整收敛只需30分23秒，vs 6L-4×4×1约2.5小时）。

用户决定：**以后统一用4层Au(111) + 3×3×1 k点作为标准生产配置**，不再用v2文档§3.2默认的6层起点。

**明确记录的风险**：v2文档§3.2原话——"四层仅可作低成本预试，不直接当正式厚度"——是对金属表面电子结构/功函数收敛性的担忧，4层是否物理上足够尚未验证。这是用户为了计算速度做出的、明确覆盖文档警告的决定，不代表4层已经过收敛检验。**待办**：有空应做4/6/8层收敛对比（v2 §3.3/§8本来就要求的敏感性测试），确认4层对目标量（EFERMI、富集积分等）是否足够，不能假设已经验证。

新标准配置汇总：4层Au(111)，3×3×1 Γ-centered k点，`ISMEAR=1/SIGMA=0.2`（Methfessel-Paxton，见G节），`ALGO=Normal`，`KPAR`不设（=1），`NCORE=8`，128核highmem分区，`I_MPI_COLL_INTRANODE=pt2pt`。**【2026-09-26 补充，见K节】这里的`ALGO=Normal`/`PREC=Accurate`/128纯MPI不再是速度上的推荐：实测PREC=Normal + ALGO=Fast + 16进程×8线程 + NPAR=16每步快4.5倍。改生产配置前需先做一次Normal vs Accurate的精度对照（K节第8条）。**

## I. 2026-09-23 9点pilot系统性QC发现

全部9点（T/V1/A1_fcc × ΔU=-0.2/0/+0.2）跑完后做的系统性检查（不只是spot-check CP loop自己报的mu_e）：

1. **CP loop内部收敛**：全部9点`|mu_e - TARGETMU| < 0.003 eV`，在`FERMICONVERGE=0.01`容差内，一致性好。
2. **SDIEL/SION平台**：全部9点在体相窗口内都精确等于1.000000，无异常。
3. **体相PHI平台平坦度随|ΔU|变差**：ΔU=0的点（T_dUp00/V1_dUp00）平台std在1e-4~1e-6eV量级（很好）；ΔU=±0.2的点平台std跳到1e-3~5e-3eV量级（差了1-3个量级）。怀疑跟带电态下双电层的实际屏蔽长度有关（20Å液区窗口对中性态够用，对~2-3µC/cm²的带电态可能不够），也可能跟下面第4点是同一个根源，**未确认**，先记录不下结论。
4. **`RHOB+RHOION`场直接积分做电荷闭合的独立校验，在带电点（ΔU≠0）上给出离谱结果**：积分出几百到上千个电子的"净电荷"，跟CP loop自己算出的~0.15-0.24电子完全不匹配；但ΔU=0的点（T_dUp00残差-0.025e，跟中性参考的-7.8e-7e同量级）完全正常。查了`RHOION`的z剖面（T_dUm02）：金属区和真空区精确为零，液区内呈现物理上合理的"表面附近大、指数衰减进入体相"双电层形状（不是文件损坏/重复写入），但绝对幅度（表面附近~2 electron/格点单位）比按1M体相密度换算的量级大了约3个数量级。~~**开放问题，未解决**~~

   **2026-09-25 已解决，根因找到并确认，不是猜测。** 直接查了本地`solvation.F`源码（1375/1384/2339/2396行等）：`n_b`/`n_ion`在写出前都乘了一次`LATT_CUR%OMEGA`（晶胞体积），紧接着直接`WRITE_TO_FILE_RC('RHOB',n_b)`/`WRITE_TO_FILE_RC('RHOION',n_ion)`写出，没有再除回去——RHOB/RHOION文件里的原始格点值本来就是`V_cell·ρ(r)`,不是`ρ(r)`本身。而`pilot_qc.py`原来的`closure = (rhob.sum()+rhoion.sum())*dV`（`dV=V/ngrid`）又多乘了一次体积,总共错了`V_cell`这个因子（对5342.6 Å³的胞,几百到上千电子的离谱数字正好对应上~0.1-0.3电子的真实闭合乘以体积后的量级）。

   修正为`closure = (rhob.sum()+rhoion.sum())/ngrid`（不再多乘`V`）后,重新对全部9点验证:残差从原来的几百电子降到**1e-6~1e-5 e量级**,跟CP loop自己的电子数变化(`dn_cp`)在符号上是`closure ≈ -dn_cp`(原脚本这个符号约定本身没错),幅度匹配到小数点后5-6位。第3条提到的"体相PHI平台随|ΔU|变差"是否同源尚未专门复查,但这条本身已经**解决**：`RHOB`/`RHOION`场现在可以正确使用（除以`V_cell`即可），PHI/SDIEL/SION不受影响（它们的写出路径没有这个体积因子，不能套用同样的除法）。修复已落到`scripts/pilot_qc.py`；v2 §10.2的K_D/N_excess计算现在可以用修正后的场，但仍需先完成第3条和几何弛豫等其他前置项。

## J. 9点pilot的电容/PZC分析（用CP loop自己报的电子数,不依赖上面第I节存疑的场积分）

用σ=(N_neutral-N_e)·e/A_cell对每个结构的3个ΔU点做线性拟合：

| 结构 | dσ/dΔU (µC/cm²/V) | 截距 σ(ΔU=0) (µC/cm²) |
|---|---|---|
| T | -11.08 | -0.01 |
| V1 | -11.28 | -0.15 |
| A1_fcc | -11.89 | +0.82 |

三者斜率彼此相差<7%——说明这套连续介质模型里"双电层电容"主要由电解液本身（1M、R_ION=4Å）决定,缺陷种类对它影响很小。真正体现缺陷差异的是**截距**：换算成PZC偏移(=-截距/斜率)，A1-fcc相对T偏移约+0.069V，V1只偏移约0.014V——即"加原子"这类凸起缺陷对局域PZC/功函数的扰动明显比"空位"缺陷大。这是本项目第一个跨结构、有统计支撑（3点线性拟合而非单点对比）的物理信号，但仍建立在**未弛豫的理想几何**上（见下条）,不是最终结论。

## K. 2026-09-26 速度问题重查:E节第4条(ALGO)、H节(KPAR)结论撤回并修正

起因:Step-16x1_dUp02(72 Au,ALGO=Normal)6小时墙钟超时未收敛;用户质疑"Fast慢6倍"与其200原子MD经验矛盾。逐项对比OUTCAR计时后,**E节第4条"ALGO=Fast比Normal慢6倍"的结论是错的,撤回**:

1. **当年的测试是坏的**:`04_speedtest/T_ismear1_algofast_kpar4`把ALGO=Fast和KPAR=4两个变量同时改,总共只跑了4步SCF、每步142 s(同体系Normal为134 s),而ALGO=Fast前5步本来就是Davidson——**那4步根本没进入RMM-DIIS阶段**,却被用来下"Fast慢"的结论。同一节还写"Fast触发更多次昂贵的PB求解",也是错的:溶剂求解(POTLOK)只占每步时间的6%。
2. **每步时间的真实去向**(Step-16x1_muref第30步,128纯MPI):LOOP 60.1 s = EDDAV 51.8(86%)+ POTLOK 3.9 + CHARGE 3.5 + FINALIZE 3.7。瓶颈是本征求解器,不是溶剂。
3. **单变量计时测试**(同一Step-16x1输入,NELM=30,拿到12~14步即scancel;每步SCF墙钟,单节点128核):

| 配置 | 每步 | 说明 |
|---|---|---|
| 基线:PREC=Accurate, ALGO=Normal, 128纯MPI, NCORE=8, pt2pt | 52~75 s(均值~64) | 生产配置 |
| ALGO=Fast(RMM阶段) | 53~55 s | 快~10%,不是慢6倍 |
| KPAR=7(112进程,每k点16进程) | 93~103 s | **更慢**:EDDAV 64 s(每k点16进程≈8×7.4 s,128进程强扩展本已近理想),POTLOK 3.9→36 s、FINALIZE 3.7→34 s——溶剂网格工作不按k点并行,每个KPAR组各算一遍。用户MD的KPAR=5之所以有效,是因为配了5个节点、一个k点独占一节点 |
| 混合并行16进程×8线程(NCORE=8→NPAR=2,无pt2pt) | 37~61 s(均值~48) | 快~20%;EDDAV仍51 s。pt2pt/进程布局不是主因。注意需`OMP_STACKSIZE=512m`,否则初始化段错误 |
| **PREC=Normal**(其余同基线) | **21~32 s(均值~26)** | **2.5倍**。ROPT -2.5e-4→-5e-4,实空间投影算符格点1052→334,细网格22.6M→8.6M;EDDAV 52→30 s,POTLOK/CHARGE/FINALIZE各降3倍。ncg不变,即每次H|ψ>作用便宜了 |
| **组合:PREC=Normal + ALGO=Fast + 16×8线程 + NPAR=16** | **Davidson 18~21 s,RMM-DIIS 14~15 s** | **4.3~4.5倍**;与用户MD的速度设置一致(`distr: one band on NCORE=1 cores, 16 groups`)。第9步分解:LOOP 15.0 = RMM-DIIS 8.8 + EDDIAG 2.8 + POTLOK 1.4 + FINALIZE 1.3 + CHARGE 0.3。ncg与基线相同(7000~7800)。核数归一1900核·秒/步,已低于MD的4100核·秒/步 |

4. **PREC=Accurate的来源**是v2文档§4.1模板(第161行),该模板自注"待核对";用户MD用PREC=normal(VASP默认)。**PREC是每步时间的主因**,不是溶剂、不是ALGO、不是KPAR、不是MPI集合通信。
5. **与用户MD的"步数"差距是另一回事**:静态单点冷启动、EDIFF=1e-7、3轮CP需200~320步SCF;MD热启动、EDIFF=1e-5每离子步~16步。muref逐步计时显示每步时间∝ncg,每轮CP更新电子数后ncg重回7000~10000再磨~100步。
6. **dUp02(Normal)超时的真实原因**:第一轮CP从中性电子数起步,与muref完全相同的问题,muref 46步收敦到-220.451 eV,而dUp02在第50步停在-218.65 eV(高1.8 eV)——SCF从随机初始波函数走上了错误轨迹(47Å×45Å长胞的电荷晃动),不是"慢",也不是ALGO问题;Fast重试(job 20915640)第一轮正常收敛到相同能量。**长胞需要slab型混合参数(AMIX/BMIX/AMIN)来固化,不能靠运气。**
7. **与MD INCAR的其余差异(影响物理不影响速度,用户2026-09-26决定:不需与MD一致,保留PBE、无IVDW、默认TAU)**:MD用GGA=RP、IVDW=12+LVDW_EWALD、TAU=9e-3(本地默认8.79e-4)、LORBIT=11。D3不进KS哈密顿量,对固定几何的密度/PHI/RHOION/电子数零影响;TAU通过空腔项进电势,有微小电子效应。
8. **待办**:采用PREC=Normal前做一次精度对照(同一Step-8x1_muref点,Normal vs Accurate的N_e、TOTEN、K_D),确认物理量不变后再改生产配置;混合并行需在job-run加`ulimit -s unlimited; export OMP_STACKSIZE=512m`。
