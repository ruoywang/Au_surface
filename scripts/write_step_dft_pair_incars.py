#!/usr/bin/env python3
"""INCAR/KPOINTS/job-run for the Step-8x1 vs Step-8x2 numerical-consistency
check. Reuses the exact template/config from write_pilot_incars.py (the
already-validated T/V1/A1-fcc batch) with three deliberate changes:
  - FERMICONVERGE 0.01 -> 0.001 (tighter target for this specific comparison)
  - TARGETMU reuses mu0=-4.9071 as-is -- an internal CP-DFT reference value
    already used throughout the pilot, NOT re-labeled as a calibrated
    1.0 V vs RHE value (SHE/RHE calibration is a separate, not-yet-done task)
  - KPOINTS paired by consistent reciprocal-space density (not copied
    unchanged): baseline is the pilot's own 3x3x1 for an ~11.76 A cell
    (N*|a|=const=~35.3). a1=23.52 A (both cells, unchanged) -> N1=2.
    a2: Step-8x1=2.94 A -> N2=12; Step-8x2=5.88 A -> N2=6.
"""
import json
import os
import shutil

TEMPLATE = """# Numerical consistency check: {tag}, same structure as its x{mult_other} pair,
# repeated x{mult} along the step direction. TARGETMU reuses the pilot's internal
# mu0 reference (-4.9071 eV) -- NOT a calibrated V vs RHE value (see defect_plan.md).
# Static, no ionic relaxation, no potential sweep, no perturbation -- single question:
# does the physical result change when the same periodic structure is repeated x{mult}?
SYSTEM = {tag}_muref

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
TARGETMU      = -4.9071
FERMICONVERGE = 0.001
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
#SBATCH --job-name={tag}_muref
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

KPOINTS_TEMPLATE = """Automatic mesh
0
Gamma
{n1} {n2} 1
0 0 0
"""

PAIRS = [
    ("Step-8x1", 1, 2, 2, 12),
    ("Step-8x2", 2, 1, 2, 6),
]

for tag, mult, mult_other, n1, n2 in PAIRS:
    manifest = json.load(open(f"03_pilot/structures/{tag}/manifest.json"))
    d = f"03_pilot/{tag}_muref"
    os.makedirs(d, exist_ok=True)
    shutil.copy(f"03_pilot/structures/{tag}/POSCAR", f"{d}/POSCAR")
    shutil.copy("/anvil/projects/x-che190065/rywang/AutoVASP2.0/potential/Au/POTCAR", f"{d}/POTCAR")
    incar = TEMPLATE.format(tag=tag, mult=mult, mult_other=mult_other,
                             sol_z0=manifest["SOL_Z0_A"], sol_z1=manifest["SOL_Z1_A"])
    with open(f"{d}/INCAR", "w") as f:
        f.write(incar)
    with open(f"{d}/job-run", "w") as f:
        f.write(JOBRUN.format(tag=tag))
    with open(f"{d}/KPOINTS", "w") as f:
        f.write(KPOINTS_TEMPLATE.format(n1=n1, n2=n2))
    print(d, f"KPOINTS={n1} {n2} 1", "TARGETMU=-4.9071 FERMICONVERGE=0.001")
