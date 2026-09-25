# 给 Claude 的完整执行指令 v2
## 单侧 CEP-DIP：Au(111) 缺陷表面的三维电解液响应、非特异性阴离子富集与 ML 数据集

本版完整替代 v1。默认采用已有 CEP-DIP 的“单侧隐式溶剂＋包含溶剂电荷/偶极的自洽偶极修正”，不再默认构造双面对称 slab。

这是一份实施方案，不是计算结果。文中的尺寸、参数扫描、预算和阈值都是第一版建议起点，必须经目标体系验证。当前已知用户的 CEP-DIP 支持单侧溶剂和 dipole correction；本文没有访问服务器当前二进制和全部源码，因此具体接口由你在本地核对后填写，不能猜测。

## 0. 任务与执行边界

请在当前授权环境中完成：环境审计 → 初始结构 → 小规模 CP-DFT 验证 → 三维场与物种密度导出 → 富集分析 → 第一批可用于训练的标签。

研究问题：

> 在相同电极电势、温度和体相盐浓度下，Au(111) 表面缺陷如何改变金属电子响应、三维电势/电场、溶剂极化，以及离子可达区域中的平均阴离子富集？

POSCAR 只包含 Au；不添加显式水、Cl、K、H、OH 或其他吸附物。使用连续介电水和对称 1:1 电解质。阴离子可以解释为指定电荷及有效尺寸下 Cl− 的非特异性近似，不能声称模型已包含 Cl 的特异性成键、脱水或电荷转移。

本阶段不计算 Au–Cl 化学吸附能、fcc/top 化学吸附占据、扩散系数、吸附速率、捕获轨迹或腐蚀。固定表面计算描述给定构型的连续介质平衡响应，不自动包含 Au 构型的有限温度统计。

优先沿用现有已验证代码和解析器，不重写成熟的 CEP-DIP 主链。建立独立工作目录；不覆盖旧结果、不修改生产二进制、不公开授权 VASP 源码/POTCAR。新增只读输出在独立分支实现并做回归测试。先 dry-run 和小规模试算，再按预算逐阶段提交；未获授权的昂贵队列任务只生成输入。

必须接受零结果：局部电场增强可能只存在于离子无法进入的区域，或在离子可达位置已低于可分辨水平。不要为了产生热点而调整参数。

## 1. 本地接口审计与已有路径复用

交付 `00_audit/environment.md`、`model_conventions.md`、`parameter_map.md`。

记录 CEP-DIP/CP-VASP/VASPsol++/VASP 版本、commit、局部补丁、二进制路径和 SHA256，编译/MPI 配置、队列预算、Au PAW 的泛函/ZVAL/ENMAX/版本和校验和。

旧项目记录曾使用 `/work/08384/tg876840/stampede3/CEP-DIP`。只把它作为待确认线索，不假设路径仍有效。优先复用当前成功算例和用户已有 `slide_workflow_cal1.md` 等审计记录。

### 1.1 恒电势控制

公开 CP-VASP v2 固定结构示例使用 `LCEP`、`NESCHEME=1`、`TARGETMU`、`FERMICONVERGE`、`CAP_MAX` [R1]。核对本地实际含义、单位、目标电势零点和更新循环。

由一套控制器更新电子数。不要同时触发 CP-VASP 和 VASPsol++ 的 `EFERMI_ref` 控制器。静态参考标签要求逐构型收敛到目标电子化学势，不使用允许费米能涨落的 MD 模式，也不使用冻结/延迟溶剂响应代替自洽标签。

### 1.2 保护完整的电子—空腔—偶极—溶剂链

核实并保留：

```text
CHGCAR / CHTOT                -> n_val
DENCOR                       -> partial-core density
n_cavity_input                = n_val + DENCOR
q_e                           = NELECT - sum(ZVAL)
phi_explicit                  <- POTHAR + POTION + CVDIP
phi_explicit + cavity + q_e    -> converged solvent response
RHOB + RHOION                 -> Qsol, Dsol -> updated CVDIP
```

这里的数组名、卷积、总势组装和符号以本地源码为准。不能用 CHGCAR 代替全部空腔输入，不能遗漏 POTION/DENCOR/CVDIP，不能用任意高斯核代替已验证的赝势路径。[R8]

本版要求继续使用包含溶剂电荷及偶极矩的闭环；不要另外叠加第二个普通 slab dipole correction，造成重复修正。输出补丁前后，用同一算例比较能量、力、电势、RHOB/RHOION；新增输出不得改变解。

### 1.3 单侧标签映射

逐个核对当前源码是否使用 `LVAC`、`SOL_Z0/SOL_Z1`、`VAC_Z0/VAC_Z1`、`SOL_SIGMA/ION_SIGMA`、`D_STERN`、`IDIPOL/LDIPOL/DIPOL` 或其他名称。记录：

- 标签是否存在、类型和默认值；
- 位置单位是 Å、分数坐标还是其他内部定义；
- 窗口作用于介电、离子还是两者；
- 如何与密度空腔相乘/组合，过渡区宽度如何定义；
- 偶极中心、修正势不连续面、Qsol/Dsol 的实际计算与更新位置；
- 能量和力是否包含相应修正及必要响应项。

不要沿用旧体系的窗口数值；也不要未经核对把公开原版 VASP 的带电修正限制机械套到已修改的、总电荷由电解液补偿的 CEP-DIP 上。[R3–R5]

## 2. 唯一默认边界：上侧电解液、下侧真空

### 2.1 结构与区域

取表面法向为 +z，第三晶格矢量沿法向，面内保留周期性。默认结构为：

```text
周期连接的真空缓冲（修正势跳跃面如存在，应位于安全真空区域）
远端溶剂—真空平滑窗口
充分恢复体相的连续电解液
电双层与非局域离子/介电空腔
带缺陷的 Au(111) 上表面
Au slab（先六层）
平整、无吸附物的 Au 背面
周期连接的真空缓冲
```

缺陷只放在上表面，背面不复制缺陷。先固定底部两层，其余 Au 可弛豫。所有力标签保留原始力及约束掩码，不把固定原子的力人为改成零。

使用真实三维密度构建的空腔描述目标表面粗糙形状；单侧窗口只负责屏蔽背面和限定远端储库。允许时，将近端开关放在金属内部，使它在全部目标界面附近已经为 1；禁止用一条水平截断面把凹坑/台阶上方的离子响应削平。若本地窗口定义不允许这样做，展示其实际影响后再确定位置。

### 2.2 真空和液层分别定义

配置中分别保存：slab 原子厚度、从最高 Au 原子到远端液相平台/窗口的距离、窗口过渡宽度、全周期中连通真空区总厚度。不要把真空和电解液长度相加后只称“vacuum”。

1 M 起点：目标侧最高 Au 以上约 40 Å 的液区尺度，连通真空总厚度约 20 Å；过渡宽度首先沿用已有稳定设置。分别检查：液区 30/40/60 Å，真空 15/20/30 Å，slab 六/八层。四层仅可作低成本预试，不直接当正式厚度。

这些长度是预算与收敛起点，不是已验证参数。添加凸起后按实际最高原子重新计算有效距离。比较同一批缺陷时采用共同 z 布置和足够覆盖最大凸起的窗口，不给每种缺陷随意选择不同近端截断。

较低浓度重新验收液层：0.1 M 可先尝试 60/80/100 Å，最终由平台与目标量收敛决定。第一版不默认 0.01 M。

### 2.3 DIPOL 不等于“跳跃面位置”

不能机械设置 `DIPOL=0.5 0.5 0.5`，也不能默认原子质心必然给出正确位置。公开 VASP 中 DIPOL 是矩计算参考中心，采用分数坐标；修正势跳跃面与它的关系需要核对具体实现。[R5]

单侧液层很长时，照搬原子质心可能让修正势跳跃面进入液区。让边界生成器根据本地实际公式定位跳跃面，并在密度/电势剖面上确认它位于无显式电荷、无活跃溶剂响应的真空区域。若采用没有显式跳跃面的其他实现，记录对应边界条件，不人为添加一个。

### 2.4 单侧专用验收

无额外施加远场时，真空远端应无残余伪斜率，液相内部应恢复体相。真空平台和溶液平台可以不同，不得强行拉齐；电双层的真实电势降必须保留。[R4]

分别延长真空、移动远端窗口、加厚液层和 slab，确认目标界面 ψ、Q、n−、富集积分稳定。窗口平移测试不能把目标界面的实际物理溶剂层删去。

再做整体平移回归：同步移动所有原子、窗口、参考中心和相关几何掩码，并保持周期关系，检查物理结果不变。不能只移动原子而把外部窗口留在原位后要求不变。

双面对称体系仅是可选交叉验证，非默认流程、非前置门槛；如已有充分回归证据，不要求再构造。不得将两种边界的数据无标签混入训练。

## 3. 初始结构生成

用可复现脚本生成 Au fcc slab，保存 parent_id、生成规则、种子、原子映射、层号、缺陷中心、约束、周期镜像距离和构型图；不手工修改坐标而不留记录。

### 3.1 体相与母结构

优先延续目标模型既定泛函/PAW。若无既定 Au 参考，先用 PBE＋标量相对论 Au PAW 求体相平衡晶格常数 a0（5–7 个体积点或可靠优化）。记录 d_nn=a0/√2、d_111=a0/√3，不直接锁定实验晶格常数后将残余应变解释为缺陷效应。

起点为 4×4 Au(111) 表面原胞、六层，共 96 Au；6×6、六层为 216 Au。这里指非正交原始表面网格；若构建库使用矩形胞，重新核对原子数和实际面积。保持正确 ABC 堆垛，不人为镜像半个晶体。

### 3.2 八类第一批结构（全部只修改上表面）

| ID | 定义 | 未弛豫原子数/建议起点 |
|---|---|---|
| T | 未重构平整台面 | 4×4×6，96 |
| V1 | 顶层删除一个 Au | 95 |
| V2 | 顶层删除两个最近邻 Au | 94 |
| V3 | 顶层删除三个构成紧凑三角形的 Au | 93 |
| V7 | 顶层删除中心及六个最近邻；不删除第二层 | 6×6×6，209 |
| A1-fcc | 在 fcc 堆垛延续位置加一个 Au | 97 |
| A1-hcp | 在 hcp hollow 上加一个 Au | 97 |
| A3 | 在三个相邻 fcc 延续位置加紧凑三原子岛 | 99 |

新加 Au 初始高度用下一 (111) 层的几何高度，再弛豫。A1-fcc/hcp 若收敛为同一结构，去重并说明实际 parent 数减少，不人为固定以保持预定数量。

V7 使用 6×6 时，必须同时生成相同 6×6 平板参考，不能用 4×4 平板做逐网格相减。任何新尺寸、层数、应变或边界设置都配相应 T 参考。

第二阶段用约 8×4 slab 顶面半个区域增加完整的一层 Au，形成条带和两个台阶。记录两条边的方向、配位和宽度；增加台阶间距须真实改变台面宽度，不是简单复制同一高指数晶胞。

本版不构造 Au(111) herringbone 重构：这是受控未重构形貌模型，不是对真实平衡表面的普遍断言。[R9]

### 3.3 几何检查与横向收敛

检查原子数、ABC 堆垛、连通性、PBC 最短距离和背面平整性。可用 d_min>0.75d_nn 作初步碰撞筛查，不代替实际结构判断。

4×4 缺陷首先代表周期阵列；用同一缺陷的 4×4→6×6，必要时更大胞检验孤立极限。不能在没有远端台面的晶胞中选一块区域称为“未受缺陷影响的台面”。

## 4. DFT 模板、溶剂参数和公共电势标尺

### 4.1 起点模板

以下是待核对模板，不是可不加检查直接提交的 INCAR。生成脚本必须拒绝剩余占位符、未知标签和重复控制器。

```text
SYSTEM = Au_EDL_v2_single_sided
PREC = Accurate
ENCUT = <max(500 eV, 1.3*ENMAX), then convergence-tested>
ISPIN = 1
ISYM = 0
LREAL = Auto
ALGO = Normal
EDIFF = 1E-7
NELM = 200
ISMEAR = -1
SIGMA = 0.05
IBRION = -1
NSW = 0
LCHARG = .TRUE.
LWAVE = .FALSE.
LVHAR = .TRUE.

LSOL = .TRUE.
ISOL = 2
C_MOLAR = 1.0
R_ION = 4.0

LCEP = .TRUE.
NESCHEME = 3
TARGETMU = <validated_target_mu_eV>
FERMICONVERGE = 0.01
CAP_MAX = 2.0

# 下列常见偶极标签仅在本地 CEP-DIP/CEP-HALF 按此接口工作时使用
IDIPOL = 3
LDIPOL = .TRUE.
DIPOL = 0.5 0.5 0.5
LVAC      = .TRUE.
SOL_Z0    = 8.0  #放在金属层内
SOL_Z1    = 35.0

```

`prepare_runs.py` 将电子、溶剂、CP、单侧窗口配置合并为一个实际 INCAR；不要假定 VASP 支持 YAML 或自定义 include。若本地 CEP-DIP 的触发方式与上面的 LDIPOL/IDIPOL 不同，按源码生成并记录差异，不叠加重复修正。

CAP_MAX 可先审阅官方小体系示例中的 2.0，但应核对单位、面积依赖和已有收敛经验，不让其失控更新电子数。[R1]

使用完整三维非线性介电/离子与非局域空腔。[R2,R3] R_ION=4 Å 是敏感性起点，不是已验证 Cl 水合半径。保存全部解析后的默认参数，包括 D_ION、R_SOLV、R_DIEL、NC_K、A_K、窗口平滑尺度、介电参数等。空腔平滑参数的真实标签必须核对，不得与电子 SIGMA 混淆。

4×4 slab 的 k 网格先用 Γ-centered 4×4×1，与 5×5×1 对照，必要时加密；大胞按倒空间密度缩放并检查，不能无验证改 Γ-only。

电子展宽与 SOLTEMP 分开。代表点检查 SIGMA=0.0257/0.05/0.10 eV，并记录所用电子自由能与熵项。能量标签应与输出力对应，不默认将零展宽外推能量和有限展宽力配对。[R6]

弛豫固定晶胞、底部两层，初始可用 IBRION=2、ISIF=2、EDIFFG=-0.02 eV/Å；实际需确认 CP 模式中优化器使用正确势。所有最终标签另做 NSW=0 自洽静态计算。必要时对小体系核查 SOC/自旋敏感性，不混用不同设置的标签。

### 4.2 建立统一电势参考

先用平整 T、1 M、q_e=0、相同单侧边界求电子化学势，按液相体相对齐得到 mu0。此时关闭电子数调节但保持隐式溶剂自洽。将其作为公共的“中性 slab 参考”，不要未经证明等同于溶液侧单一界面的 PZC：单侧不对称 slab 的总净电荷为零不必然意味着每个表面分区的电荷为零。

定义：

```text
mu_target[eV] = mu0[eV] - ΔU[V]
```

前提是已核对 TARGETMU 的内部零点和能量单位。pilot 用 ΔU=−0.2,0,+0.2 V；主扫描用 −0.4,−0.2,0,+0.2,+0.4 V。范围过强导致模型拥挤/异常时先缩小，不强行跑完。

所有缺陷和尺寸采用同一组物理电子化学势；不要对每个缺陷以自身中性点重新设零。更换盐浓度后保留同一标尺，记录可能的模型参考差异。不同体系相同 TARGETMU 字符串不保证相同实际标尺，必须检查体相对齐。

分析以液相体相为参照，背面真空平台只用于边界诊断。未完成模型特定参照标定，不直接标精确 V vs SHE/RHE。[R3]

建议参考几何流程：中性 T 弛豫并定义 mu0 → 在公共 ΔU=0 下弛豫各缺陷 → 固定该几何扫描电势 → 对代表性缺陷在正负偏压再次弛豫，区分纯充电响应与结构响应。

## 5. 单侧数值验收与 9 点 pilot

先做平板的截断、k 点、网格、slab/液区/真空厚度和窗口检查，再运行 T、V1、A1-fcc × 三个电势，共 9 个核心点。边界/收敛附加点另列预算，不把“9 点”当作总成本。

液相平台必须在窗口内部，远离 Au 界面和远端溶剂—真空过渡。建议有至少约 5 Å 厚的区域满足：

- SION、SDIEL>0.999，n+/n_bulk 和 n−/n_bulk 均在 1% 内；
- 平面平均电势变化<1 mV，残余场可先以 1e−4 V/Å 为检查尺度；
- 该区域不与修正势跳跃面、截断过渡区重叠。

体相验证同时查看横向起伏，不只看平面平均。真空平台在独立掩码内检查；真空和液相平台不要求相等。禁止扣除拟合斜率来“修复”失败平台。

每帧必须检查：

1. 电子 SCF、PB、自洽偶极和 CP 目标分别收敛，最终所有场来自同一收敛态。
2. 密度积分对应 NELECT；物种重建对应 RHOION；包含必要溶剂项的全胞总物理电荷闭合。
3. 建议总电荷残差<1e−3 e；用完整定义检查，既不漏掉 RHOB 的边界贡献，也不重复计数。
4. 浓度非负、体相占据合法、没有金属内部/背面真空的异常溶剂泄漏。
5. 目标界面与窗口、真空、slab 厚度分离；改变这些数值设置不改变结论。
6. 同一几何从不同初始电子数重启可复现；失败或亚稳解需保留诊断。

建议主要指标收敛：K_D 变化<3%，目标界面电势差变化约<1–2 mV，非零电荷变化<2%；Q 接近零时可先用 |ΔQ|<0.005 e/胞再按精度需求收紧。实际信号小于这些阈值时继续收紧或报告无法分辨，不能宣称微弱排序可靠。

正式标签优先将 FERMICONVERGE 收紧到 0.001 eV；最终目标根据误差敏感性确定。并发默认不超过 2 个作业，每项自动重试不超过 2 次；更大并发须在预算配置中明确批准。不要在登录节点直接跑重计算。

交付 pilot_summary：实际算例、设置、误差、耗时、内存、磁盘和失败原因。pilot 未通过不得批量生成生产标签。

## 6. 必须导出的参考量

保存结构/标量于 manifest、extxyz 或等效格式；网格存 HDF5/NPZ，记录格矢、原点、维度、存储顺序、单位和积分权重。先估算每帧与整个数据集的磁盘成本。保留真值网格；压缩或下采样必须额外验证，不能静默替换参考标签。

```text
geometry:
  cell, positions, species, layer_ids, fixed_mask, parent_id, seed
  defect_centers/types, boundary_mode="single_sided_CEP_DIP"
state:
  NELECT, NELECT_neutral, q_e, Q_phys
  mu_target, raw_EF, aligned_mu, bulk_reference
energy:
  all_raw_terms, electronic_entropy, F_effective, Omega_and_definition
  solvent_terms, dipole_terms, raw_forces_all_atoms
fields:
  n_val, DENCOR_or_exact_reconstruction_data, n_cavity_input
  phi_explicit, phi_solv, phi_used_by_ions
  CVDIP_or_exact_reconstruction_data, Qsol, Dsol
  SVDW, SDIEL, SION, single_side_windows
  RHOB_raw, RHOION_raw, physical_charge_conversions
  n_plus, n_minus
metadata:
  XC/PAW/code hashes, solvent parameters, temperature/concentration
  dipole reference/jump geometry, vacuum/liquid/bulk masks
  SCF/PB/CP/dipole residuals, input hashes, status and reasons
```

LVHAR 官方含义是 ionic＋Hartree 势，不等于可不审计地把任意 LOCPOT 当成物种响应的驱动势。代码可能另有 PHI/PHI_SOLV、卷积、常数位移或空腔导数势；明确每个文件对应哪个内部量，离子不能使用含 XC 的总 KS 局域势替代静电势。[R3,R7]

CHGCAR 的归一化和数组顺序用电子数积分校验；需要 PAW 相关字段时保留合法重建信息。后续 MD 抽帧一律静态重算，不能把外推密度当自洽标签。[R7]

## 7. 统一电荷、电势、电场和浓度

### 7.1 符号与参考

```text
q_e = NELECT - sum(ZVAL)             # excess electrons
Q_phys/e = -q_e                     # explicit slab net physical charge
psi = phi_physical - phi_bulk       # V
E = -grad(psi)                      # V/Å
rho_ion_phys/e = n_plus - n_minus    # Å^-3
```

保留 raw 和 converted。VASPsol++ 使用 electron-positive 约定；若 raw potential 是电子势能，须按真实定义变号/换单位。[R2,R3] 先验证正向增加 ΔU 对电子数的方向，再解释局部分布，不预设任意缺陷一定富集阴离子。

### 7.2 物种浓度

优先直接导出内部正负离子占据。若没有输出，在独立只读补丁添加，或用完全相同的本构关系重算，必须重建出原来的 RHOION。

对本地实现确实对应的等尺寸、对称 1:1 格点气体，可写：

```text
u = e*psi/(kB*T)
D = 1 + (2*n_bulk/n_max)*(cosh(u)-1)
n_plus  = SION*n_bulk*exp(-u)/D
n_minus = SION*n_bulk*exp(+u)/D
```

若有额外平滑或电化学势位移，采用内部实际驱动变量，而不是直接套未平滑 PHI。n_max、D_ION 与 R_ION 的关系从代码获得，不自行假定。[R2,R3]

数值上使用稳定表达式防溢出，不随意截断浓度。换算：

```text
n[Å^-3] = 6.02214076e-4 * c[mol/L]
```

不能令 n_minus=−RHOION/e，也不能把 RHOB+RHOION 的负区域全部解释成 Cl−。

### 7.3 电场

按完整晶格矩阵/倒格矢求 E=−∇psi，非正交面内坐标不可当作笛卡尔坐标。

特别注意带跳跃的 dipole potential：不要对包含不连续面的完整数组直接做三维周期 FFT 导数，再把振铃当热点。可在光滑子区域差分；或按实际解析形式分离修正项，对光滑部分求导，再加回修正势在目标区域的解析导数。用两种方法交叉检查目标区域。

输出横向/法向场、方向、电势差和可达区场强百分位数。不以原子核附近 E_max 排缺陷；修正面与远端窗口附近的数值场不进入目标界面统计。

## 8. 主参考矩阵及参数敏感性

pilot 后运行至多八类唯一 parent × 五电势，约 40 个主点，外加必需的匹配尺寸平板参考。先固定各 parent 的公共参考几何，再对 2–3 类代表缺陷做正负偏压弛豫。高成本 V7、台阶和大胞后置；实际结构数和重复去重后重新记账。

敏感性先选 T、V1、A1/台阶三类，不做所有参数全排列：

- 盐浓度 1.0/0.3/0.1 M；每种浓度重新验收液层和电势参考。
- R_ION=3/4/5 Å，作为模型有效接近距离的不确定性测试，不作为实测 Cl 半径。
- 能独立控制时，固定 D_ION 改 R_ION 检查可达性；再独立改拥挤尺寸。若参数联动，说明同时改变了什么。
- 必须满足体相 2*n_bulk<n_max；否则参数组合无效，不能靠 SCF 强行计算。
- 选少量空腔平滑/网格测试，分清数值误差与物理模型敏感性；不为热点调参。

## 9. ML 数据集：第一版约 300–600 个合格标签

### 9.1 环境条件

基础训练集固定 1 M、298 K、一套溶剂/边界模型。盐浓度、离子尺寸等变化独立归档；只有模型通过输入或显式 PB 模块接收这些条件时，才合并训练。

模型对应函数应明确依赖 R、Z、q、cell 和环境边界 B。不能让同一未包含环境条件的输入对应互相冲突的标签。窗口/参考中心/外部法向的处理属于模型定义；坐标增强时必须同步处理，不能随意旋转原子但保持 z 定向边界不变。

### 9.2 构型采样

每个唯一 parent：五个电势参考点；可动 Au 的高斯随机位移，每分量标准差 0.02/0.05/0.10 Å，每档 3 个种子；每个扰动几何做 −0.2/0/+0.2 V 静态标签。

底部固定层不扰动；不要求上下镜像；窗口和 slab 锚定方式保持一致。筛掉原子碰撞/非预期断裂。未弛豫扰动构型不应先完全最小化后才打标签，否则会丢失力训练需要的非平衡覆盖。

每 parent 约 5+9×3=32 点，8 个 parent 约 256 点。加入匹配平板、少量 ±1%/±2% 面内应变、不同超胞、双缺陷和主动学习难例后，控制在 300–600 个合格标签。不同应变需记录固定 q 还是固定面电荷密度，不能混淆。

优先“相同几何、多种 q/电势”和“相似局部结构、不同远端缺陷排列”。必要时在相同 R 上增加固定电子数邻点以检查充电曲率；重新求电子/溶剂自洽。金属-only MD 仅用于几何候选采样，抽帧后重新静态打标签，不报告真实水界面动力学。

### 9.3 能量—力一致性和 Legendre 定义

从源码确认当前 TOTEN 是否已经包含电子储库修正。VASPsol++ 自带 CP 与 CEP-DIP/CP-VASP 的报告规则未必相同，不凭名字判断。[R3]

同时保留有效固定电子数自由能 F、mu_e、NELECT、q_e、所有原始项与明确的 Omega：

```text
Omega_N = F - mu_e*NELECT
Omega_q = F - mu_e*q_e
```

二者对固定组成在相同 mu 下只差常数，跨 Au 原子数或跨 mu 比较时不能混同。以 (R,q) 为输入的能量模型通常训练对应 F，不把已做 Legendre 变换的 Omega 误当 F，也不重复减 mu*q。不额外再加一次训练目标已包含的溶剂能。

选择 T/V1/A1，在代表性 Au 的 x/z 方向做 ±0.01、±0.02 Å 有限差分：固定 q 比 F 的导数，固定 mu 比一致定义 Omega 的导数。边界条件保持物理一致，不能每次位移任意重定窗口/参考面而不计响应。

建议力差约<0.01 eV/Å，依据实际精度需求收紧。未经能量—力校验的数据只能做场分析，不能自动进入保守 MLFF 训练。固定层仍保留原始力；禁止靠全局净力强制归零、漂移消除或恒温器掩盖修正项错误。未验证溶剂应力，不训练 stress。

### 9.4 划分和评价

同一 parent 的电势/扰动/相邻帧按组划分，禁止近重复几何泄漏。八个 parent 不足以声称广泛泛化，可采用留一缺陷族验证加独立大胞/双缺陷/台阶挑战集。记录实际覆盖范围。

除 E/F/q/mu 误差外，还报告离子可达区的势误差、RHOB/RHOION 的积分与空间误差、K_D 和作用范围误差。全空间 RMSE 可能被巨大真空体积稀释，分别给目标界面、体相、真空和不可达区指标；不能只给一个全网格误差。

## 10. 最终分析指标与图组

### 10.1 单侧 QC 图

平面平均 n_val、RHOB、RHOION、psi、n+、n−；叠加 slab、SION/SDIEL、液相平台、窗口过渡、真空和修正面标记。另给 q/mu/PB/偶极闭环收敛历史与 Q–ΔU 曲线。

图中零点与坐标统一；不为了让曲线对齐而独立改变物理参考。相同 z 截面与相对离子可达界面相同距离的截面都给出，分辨几何高度效应。凹坑/台阶的等值面多值时，用三维最近距离/符号距离，不强行定义单值 z(x,y)。

### 10.2 富集及离子过量

在看到结果之前预定义缺陷邻域 Omega_D：

```text
V_access,D = integral_D SION dV
N_minus,D  = integral_D n_minus dV
K_D = N_minus,D/(n_bulk*V_access,D)
N_excess,D = integral_D (n_minus - n_bulk*SION) dV
```

这是本项目的可达体积归一化统计定义，不自动等于某个唯一实验 Gibbs 表面过量。报告 K_D、K_D/K_T、V_access、N_minus 和 N_excess。V_access 太小时标不可判定；不能用近零分母放大热点。仅对目标金属/电解液界面积分，不包含背面真空或远端人工液—真空窗口。

初始可按横向半径 0–3/3–6/6–9 Å 和距离子可达界面 0–5/5–10 Å 分区；必须小于周期无歧义尺度。在小胞中缩小区域或扩大胞，不能把跨周期距离当独立空间。台阶按距边缘的条带统计并按长度归一化。

### 10.3 全局充电与局部效应

相同 mu 下不同缺陷可以有不同总 Q。相对体相的富集与同一表面内部“缺陷/远端台面”富集分开报告，缺少真正远端时用匹配尺寸平板及周期阵列定义。

电子差分优先比较同几何不同电势，不能将赝势核心峰当表面额外电荷。不同原子数的图需清楚说明参考和对齐。必要时少量匹配总电荷对照，但单独标记 fixed_Q，不与 fixed_mu 的排序混用。

### 10.4 单侧面积与电容

A_proj=|a1×a2|。总电极充电响应可记录 sigma_cell=Q_phys/A_proj，不除以 2。

但 sigma_cell 不自动等于溶液侧局域电荷。通过同几何充电差分、累积电荷/高斯面或收敛的前后分区检查背面贡献与薄 slab 耦合。绝对 Bader 分区不是唯一表面电荷定义，说明方法。不能简单令金属正面电荷等于负的近界面 RHOION 积分而忽略其他贡献。

若背面增量响应可忽略且厚度收敛，再将 C=d(Q/A_proj)/dU 解释为目标界面的充电电容；否则明确称整个单侧计算模型的响应。给导数窗口和拟合误差；不把静态连续模型结果称为包含 Cl 化学吸附赝电容的实验预测。

### 10.5 作用范围和双缺陷

扩大横向胞测 psi/n− 起伏随高度和距离的衰减，给出在数值噪声/参数不确定性之上可分辨的范围，不叫动力学捕获半径。

用相同 cell、mu、边界、网格、原子参考坐标比较 AB、A、B、0：

```text
delta_n   = n_AB   - n_A   - n_B   + n_0
delta_psi = psi_AB - psi_A - psi_B + psi_0
```

先做冻结几何匹配四态，再讨论弛豫影响。对不同 SION 同时报告共同可达区和各自实际积分。非零 delta_n 可能来自 PB 非线性、空腔重叠或总充电差，不是模型未包含的离子关联或化学协同。

### 10.6 简单模型基线

粗糙度影响电双层已有理论基础 [R10]。本项目要量化原子级电子响应与三维空间分辨率的额外贡献。

先做横向平均诊断，明确它不是自洽 1D 求解。资源允许时做等势导体＋同一电解质模型的几何基线：只在平板上标定一次参考，不逐缺陷拟合；匹配边界、浓度、可达性和拥挤定义，列清实际改变了什么。不存在可用求解器时如实列为后续任务，不把真空电势后处理冒充这个基线。

## 11. 预期结果、结论范围与近似

数值上必须看到：体相恢复、正确符号、物种电荷重建、总电荷闭合、单侧无伪远场、重复收敛和尺寸稳定。上下界面不再要求相同，因为本模型有意不对称。

科学上允许：单原子缺陷几乎不影响可达区；更宽形貌才有差异；凸起/凹坑的排序取决于电势和有效尺寸；缺陷影响近似独立或明显非加和；简单几何已足够或原子级电子响应仍必需。不能提前指定结论。

最终报告明确三个层次：

1. 数值事实：代码/网格/边界是否收敛，输出是否自洽。
2. 模型内规律：指定连续介质、电势和表面构型下的富集与响应。
3. 真实体系外推：分子水、离子特异性、脱水、真实表面重构等未包含，需要独立验证。

禁止下列推断：K_D>1 就形成 Au–Cl 键；E_max 最大所以 c 最大；静态富集范围就是捕获半径；同隐式 DFT 一致就已验证真实 Cl/水；不同 Au 原子数总能量更低就证明缺陷更稳定。若研究缺陷形成热力学，另定义 Au 储库和一致的恒电势势能，不能偷用本阶段原始总能量排序。

## 12. 工程交付与停止条件

实现：

```text
00_audit/{environment,model_conventions,parameter_map,output_patch_regression}.md
configs/v2_single_sided.yaml
scripts/build_structures.py
scripts/build_boundary.py
scripts/prepare_runs.py
scripts/parse_fields.py
scripts/analyze_enrichment.py
scripts/build_dataset.py
tests/test_units_and_signs.py
tests/test_charge_and_species.py
tests/test_single_sided_boundary.py
tests/test_field_derivatives.py
tests/test_force_energy.py
manifests/{structures,calculations}.jsonl
reports/{pilot_summary,convergence,reference_results,dataset_card}.md
RUN_STATUS.json
```

最小配置应包含：

```yaml
boundary_mode: single_sided_CEP_DIP
explicit_elements: [Au]
slab_layers: 6
fixed_bottom_layers: 2
base_supercell: [4, 4]
solvent_temperature_K: 298
bulk_salt_molarity: 1.0
R_ION_A: 4.0
delta_U_pilot_V: [-0.2, 0.0, 0.2]
delta_U_main_V: [-0.4, -0.2, 0.0, 0.2, 0.4]
liquid_extent_initial_A: 40.0
connected_vacuum_initial_A: 20.0
pilot_geometries: [T, V1, A1_fcc]
seed: 20260921
max_concurrent_jobs: 2
max_retries_per_case: 2
production_label_cap: 600
submit_jobs: false  # 仅在已有授权和预算确认后开启
```

上面是研究配置，不是 INCAR；所有实际单侧标签、DIPOL、窗口宽度、目标 mu、运行路径/队列由已验证生成器补齐。输入占位符/未知标签/来源不明默认值一律拒绝提交。

程序支持 dry-run、幂等性、固定种子、合法续算和失败隔离。几何/cell/XC/PAW/溶剂或边界条件不同的重启文件须验证兼容，不无差别复制。保留失败算例与变更记录；禁止无限重试。

停止生产的情形：无液相平台；存在修正面穿过活跃液区；符号/物种重建错误；总电荷闭合失败；窗口/厚度误差大于研究信号；能量标签定义不清或力验证失败。场标签合格但力失败时可独立存为场分析数据，不混入能量保守训练集。

第一阶段的合格交付是：单侧 T/V1/A1 的 9 个核心点及必要附加收敛测试可复现，从电子数、CVDIP、三维溶剂场到物种浓度、富集指标和能量/力标签都可追溯。不要以目录数量或正常退出作完成标准。

## 参考资料与使用范围

[R1] CP-VASP 官方 README 和 Version 2 / Example Inputs / Fixed Structure / INCAR。仅用于公开接口对照，本地 CEP-DIP 修改版仍以实际源码/日志为准。

[R2] S. M. R. Islam, F. Khezeli, S. Ringe, C. Plaisance. An implicit electrolyte model for plane wave density functional theory exhibiting nonlinear response and a nonlocal cavity definition. J. Chem. Phys. 159, 234117 (2023), DOI: 10.1063/5.0176308. 公开预印本 arXiv:2307.04551v5。用于非局域空腔、非线性电解质和场/浓度定义。

[R3] VASPsol/VASPsol 官方 README：溶剂参数、R_ION/D_ION、输出能量和场文件说明。不要将其输出规则未经验证套用到 CEP-DIP。

[R4] L. Bengtsson. Dipole correction for surface supercell calculations. Phys. Rev. B 59, 12301 (1999), DOI: 10.1103/PhysRevB.59.12301。用于不对称 slab 的周期伪场处理原则，不是用户改版所有实现细节的证明。

[R5] VASP Wiki: LDIPOL, IDIPOL, DIPOL；GPAW dipole-layer correction tutorial。用于一般偶极接口和边界诊断；本地带电补偿改版的适用范围另行审计。

[R6] VASP Wiki: Smearing technique。用于能量—力与电子展宽一致性。

[R7] VASP Wiki: LVHAR, CHGCAR。用于势文件语义、密度归一化和 MD 外推密度警告。

[R8] 用户项目记录 slide_workflow_cal1.md。记录了 CHGCAR/POTION/DENCOR/CVDIP 与 Qsol/Dsol 自洽闭环的单构型复现；不等于当前所有 Au 体系已经验证。

[R9] F. Hanke, J. Björk. Structure and local reactivity of the Au(111) surface reconstruction. Phys. Rev. B 87, 235422 (2013), DOI: 10.1103/PhysRevB.87.235422。用于未重构模型的适用范围。

[R10] L. I. Daikhin, A. A. Kornyshev, M. Urbakh. Double-layer capacitance on a rough metal surface. Phys. Rev. E 53, 6192 (1996), DOI: 10.1103/PhysRevE.53.6192。用于粗糙电双层的几何理论基线。

公开地址：

```text
https://github.com/yuanyue-liu-group/CP-VASP
https://github.com/yuanyue-liu-group/CP-VASP/blob/main/Version%202/Example%20Inputs/Fixed%20Structure/INCAR
https://arxiv.org/html/2307.04551v5
https://github.com/VASPsol/VASPsol
https://doi.org/10.1103/PhysRevB.59.12301
https://vasp.at/wiki/index.php/LDIPOL
https://vasp.at/wiki/index.php/IDIPOL
https://vasp.at/wiki/index.php/DIPOL
https://gpaw.readthedocs.io/tutorialsexercises/electrostatics/dipole_correction/dipole.html
https://vasp.at/wiki/index.php/Smearing_technique
https://vasp.at/wiki/index.php/LVHAR
https://vasp.at/wiki/index.php/CHGCAR
https://doi.org/10.1103/PhysRevB.87.235422
https://doi.org/10.1103/PhysRevE.53.6192
```
