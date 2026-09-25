#!/usr/bin/env python3
"""Write the 9-point pilot INCARs (T/V1/A1_fcc x DeltaU=-0.2/0/+0.2), using
mu0 from the corrected T_neutral_v3 run (4-layer, 3x3x1, ISMEAR=1/SIGMA=0.2 --
project standard config decided 2026-09-22, see 00_audit/parameter_map.md H).
Call with --mu0 <value> once known."""
import argparse
import json

TEMPLATE = """# Pilot point: {tag}, DeltaU={du:+.1f} V, TARGETMU = mu0({mu0}) + DeltaU = {targetmu:.4f} eV
# Structure: {defect_type}, a0=4.158 A (locally fitted, 02_bulk_eos), shared
# single-sided window across this T/V1/A1_fcc batch (v2 Sec.2.2).
# Proven-fast config: 128 cores (highmem), NCORE=8/NPAR=16/KPAR=1,
# I_MPI_COLL_INTRANODE=pt2pt (job-run), ALGO=Normal.
SYSTEM = Au_{tag}_dU{du_tag}

NCORE = 8

PREC   = Accurate
ENCUT  = 500
ISPIN  = 1
ISYM   = 0
LREAL  = Auto
ALGO   = Normal
EDIFF  = 1E-7
NELM   = 200
ISMEAR = 1
SIGMA  = 0.2
IBRION = -1
NSW    = 0
LCHARG = .TRUE.
LWAVE  = .FALSE.
LVHAR  = .TRUE.

LSOL    = .TRUE.
ISOL    = 2
C_MOLAR = 1.0
R_ION   = 4.0

LCEP          = .TRUE.
NESCHEME      = 3
TARGETMU      = {targetmu:.4f}
FERMICONVERGE = 0.01
CAP_MAX       = 2.0
NEADJUST      = 1

IDIPOL = 3
LDIPOL = .TRUE.
DIPOL  = 0.5 0.5 0.5

LVAC      = .TRUE.
SOL_Z0    = {sol_z0}
SOL_Z1    = {sol_z1}
SOL_SIGMA = 0.8
"""

JOBRUN = """#!/bin/bash
#SBATCH --job-name=Au_{tag}_dU{du_tag}
#SBATCH -o myjob.o%j
#SBATCH -e myjob.e%j
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --partition=highmem
#SBATCH --account=CHE190065

module purge
module load intel/19.0.5.281 impi/2019.5.281
module load intel-mkl hdf5 python
module list

export I_MPI_COLL_INTRANODE=pt2pt
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
PW=/anvil/projects/x-che190065/rywang/CEP-HALF/bin
runvasp="mpirun -np $SLURM_NTASKS $PW/vasp_std  > log.out"

eval $runvasp
"""

KPOINTS = """Automatic mesh
0
Gamma
3 3 1
0 0 0
"""

DU_MAP = {"m02": -0.2, "p00": 0.0, "p02": 0.2}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mu0", type=float, required=True)
    args = ap.parse_args()

    import os
    for tag in ["T", "V1", "A1_fcc"]:
        manifest = json.load(open(f"03_pilot/structures/{tag}/manifest.json"))
        for du_tag, du in DU_MAP.items():
            d = f"03_pilot/{tag}_dU{du_tag}"
            os.makedirs(d, exist_ok=True)
            if not os.path.exists(f"{d}/POSCAR"):
                import shutil
                shutil.copy(f"03_pilot/structures/{tag}/POSCAR", f"{d}/POSCAR")
            if not os.path.exists(f"{d}/POTCAR"):
                import shutil
                shutil.copy("/anvil/projects/x-che190065/rywang/AutoVASP2.0/potential/Au/POTCAR", f"{d}/POTCAR")
            targetmu = args.mu0 + du
            incar = TEMPLATE.format(tag=tag, du=du, du_tag=du_tag, mu0=args.mu0,
                                     targetmu=targetmu, defect_type=manifest["defect_type"],
                                     sol_z0=manifest["SOL_Z0_A"], sol_z1=manifest["SOL_Z1_A"])
            with open(f"{d}/INCAR", "w") as f:
                f.write(incar)
            with open(f"{d}/job-run", "w") as f:
                f.write(JOBRUN.format(tag=tag, du_tag=du_tag))
            with open(f"{d}/KPOINTS", "w") as f:
                f.write(KPOINTS)
            print(d, "TARGETMU=", targetmu)


if __name__ == "__main__":
    main()
