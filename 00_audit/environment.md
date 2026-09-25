# 00_audit/environment.md

核对时间：2026-09-21。核对方式：直接读取 `/anvil/projects/x-che190065/rywang/CEP-HALF` 源码/二进制、`/anvil/projects/x-che190065/rywang/AutoVASP2.0/potential/Au` POTCAR、旧算例 `/anvil/projects/x-che190065/rywang/surface_charge/6-Au/1-32_water_GCE`，均为只读操作，未修改这些目录下任何文件。

## 1. 代码与二进制

- 位置：`/anvil/projects/x-che190065/rywang/CEP-HALF`（本项目语境下这就是"CEP-DIP"）。
- 基础代码：VASP 6.3 develop 分支，`src/version.F` 内嵌日期字符串 `27Jun22`。
- 本地补丁：
  - `vaspsol++-vasp_6.3.2.patch`（VASPsol++ 官方补丁）
  - `src/cp-vaspsol++.patch`（CP-VASP 恒电势补丁）
  - 用户自加的单侧窗口扩展：`solvation.F` 中多处标注 `! LVAC Ruoyu`，确认 `LVAC/SOL_Z0/SOL_Z1/SOL_SIGMA` 及派生的 `ION_Z0/ION_Z1/D_STERN` 是本地在 VASPsol++ 基础上追加的，非官方/非 CP-VASP 上游代码。
  - `dipol.F` 中官方 `CDIPOL` 例程被修改为读取 `solvation_moments` 模块的 `Qsol_cache/Dsol_cache`，把溶剂电荷/偶极矩并入官方偶极修正（细节见 `model_conventions.md` 1.2 节）。
  - `difflog` 是相对上游 `release.6.3` 的通用代码差异记录（elphon/fileio/fock/greens_real_space/hamil_lr/nonlr/rhfatm/version 等文件的小改动），与溶剂化/CP 功能无关，未深入分析。
- 二进制：`bin/vasp_std`，构建时间 2025-10-21，
  `sha256 = 5e7152b78fcb27332156d7dacdb8495e38dd6736499f47460b2ec5032a702c41`。
- 编译配置：`arch/` 下有 intel/gnu/aocc/nvhpc 等多套 `makefile.include`；当前 `job-run` 实际使用 `module load intel/19.0.5.281 impi/2019.5.281 intel-mkl hdf5 python`，与 Anvil 上现有 `bin/vasp_std` 的构建环境一致（对应 `arch/makefile.include.intel*` 系列，未逐字节核对 build 时使用的具体是哪一份，属于低风险待办）。
- `build/std/` 下有一份完整构建目录（.o/.mod 等），未展开核对。

## 2. Au PAW

- 路径：`/anvil/projects/x-che190065/rywang/AutoVASP2.0/potential/Au/{POTCAR,PSCTR}`。
- `TITEL = PAW_PBE Au 04Oct2007`，`LEXCH = PE`（PBE），`ZVAL = 11.000`，`ENMAX = 229.943 eV`，`ENMIN = 172.457 eV`，`RPACOR = 2.33 Å`（partial core，对应文档中 `DENCOR` 的来源半径），`RCORE = 2.500`，`LPAW = T`。
- POTCAR 文件头自带官方校验码：`SHA256 = d0044ae04e2bdce24051b198fc5c053d722a5bc6fe3c3b100514a13fc5d2db88`（VASP 官方发布值，确认这是未经修改的标准 PAW_PBE Au 势）。
- 本地文件整体 `sha256sum` = `ec1ee5aca8475a9214c670a1c2127468662b62a6adf2550f15324cbe139af665`（含文件头注释，故与上面官方内嵌值不同，属正常现象，两者都已记录备查）。
- ENCUT 起点按 v2 文档 `max(500, 1.3*ENMAX)` = `max(500, 298.9) = 500 eV`，与模板一致。

## 3. 作业脚本与队列

- 当前目录 `job-run`：`--account=CHE190065`，`--partition=wholenode`，5 节点×16 任务×8 线程，`time=96:00:00`，直接 `mpirun ... $PW/vasp_std`，`PW` 指向上面确认的 `CEP-HALF/bin`。可作为新任务脚本模板，但 **96 小时 wholenode×5 节点是相当大的资源申请，pilot 阶段不应照搬，需要单独定义小规模测试用的资源请求**。
- `sacctmgr` 确认账号 `CHE190065` 下有 `cpu`、`gpu`、`ai` 三个 QOS；未见任何 SU/预算上限的机器可读记录，**具体已批准的计算预算仍需你确认**，不是代码审计能回答的问题。

## 4. 旧参考算例

- `/anvil/projects/x-che190065/rywang/surface_charge/6-Au/1-32_water_GCE`：真实生产计算，输出完整（`PHI/PHISOLV/RHOB/RHOION/SION/SDIEL/SSOLV/SCAV/SVDW/VSOLV/ELOC/P` 等场文件齐全，另有 `CHGCAR/LOCPOT/OUTCAR/vasprun.xml/vaspout.h5` 等标准 VASP 输出）。
- 其 `INCAR` 用的是 `ISOL=2`、`LSOL=.TRUE.`、`LCEP=.TRUE.`、`NESCHEME=5`、`LVAC=.TRUE., SOL_Z0=8.0, SOL_Z1=45.0`、`LDIPOL=.TRUE., IDIPOL=3, DIPOL=0.5 0.5 0.5`、`IBRION=0`（MD）——**这是一次 MD 恒电势计算，不是静态单点**，`NESCHEME=5`（Nose-Hoover 链）是它专属于 MD 的合理选择，不能被 v2 的静态工作流照搬（见 `model_conventions.md` 1.1 节）。
- 按你的补充说明，这个算例只验证了"这一个构型"的电荷/空腔/偶极闭环，不代表当前 Au 缺陷体系、不代表能量-力一致性已验证——本审计仅把它当作标签含义和输出格式的核对样本，不作为当前项目的正式回归基准。

## 5. 计算环境（Python 等）

- `module load intel/19.0.5.281 impi/2019.5.281 intel-mkl hdf5 python` 后得到 `python/3.9.5`（Anvil spack 安装），自带 `numpy 1.19.5`、`scipy 1.6.0`、`matplotlib 3.9.3`、`ase 3.23.1b1`、`h5py 3.7.0`。
- 没有 `pymatgen`；ASE 已够用于本项目结构生成/IO，如后续确实需要 pymatgen 会在当前目录本地装虚拟环境，不改动全局模块。

## 6. 存储与配额

- `/anvil/scratch/x-rywang`：308.9 GB / 100 TB，空间充裕。
- `/anvil/projects/x-che190065`：4.1 TB / 5.0 TB（82.5%），**文件数配额已 100%（1.0M/1.0M inode 已用满）**——这意味着该 projects 目录下无法再新建任何文件；由于我们只在 `/anvil/scratch/x-rywang/Au_Cl` 下工作、只读访问 projects 下的代码/POTCAR/旧算例，不受影响，但如果之后要把大规模数据集写回 projects 目录需要另外处理配额问题。
