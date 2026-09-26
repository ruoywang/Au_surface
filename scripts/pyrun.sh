#!/bin/bash
# Python launcher for the analysis scripts (2026-09-26). The user-site
# ~/.local/lib/python3.9/site-packages currently has a broken numpy (dist-info
# for 1.22.4 and 1.26.0, no package) which makes its matplotlib fail against the
# conda base numpy 1.21.5. Nothing outside Au_Cl is modified: user-site is
# disabled and only ASE (pure python, needed and absent from conda base) is
# reached through the .pyshim symlink inside this project directory.
PY=/anvil/projects/x-che190065/rywang/anaconda3/bin/python3
export PYTHONNOUSERSITE=1
export PYTHONPATH=/anvil/scratch/x-rywang/Au_Cl/.pyshim${PYTHONPATH:+:$PYTHONPATH}
exec "$PY" "$@"
