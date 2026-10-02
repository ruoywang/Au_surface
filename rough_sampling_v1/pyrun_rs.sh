#!/bin/bash
# Launcher for rough_sampling_v1 scripts (2026-10-02). Same base interpreter and the same user-site exclusion as
# scripts/pyrun.sh; the only difference is that rough_sampling_v1/env/pylib (dscribe 2.1.2 and the numpy 2.0.2 /
# ase 3.26.0 it was installed against, pip --target, project-local) is put first on PYTHONPATH. The cluster's
# inherited PYTHONPATH is dropped for the reason recorded in scripts/pyrun.sh.
PY=/anvil/projects/x-che190065/rywang/anaconda3/bin/python3
export PYTHONNOUSERSITE=1
export PYTHONPATH=/anvil/scratch/x-rywang/Au_Cl/rough_sampling_v1/env/pylib:/anvil/scratch/x-rywang/Au_Cl/.pyshim
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
exec "$PY" "$@"
