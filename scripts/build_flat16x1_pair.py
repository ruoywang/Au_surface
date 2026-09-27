#!/usr/bin/env python3
"""Flat 64-atom reference slab in the Step-16x1 cell, at the two Step-16x1 work points.

Approved 2026-09-27 (user, item (a)): remove the 8 added strip Au atoms from the
Step-16x1 production input, keep the 64 base atoms at their ORIGINAL positions, keep the
cell, the single-sided solvent window (SOL_Z0/SOL_Z1), DIPOL, k-points (1x12x1), PREC=Normal,
ENCUT, smearing and the run configuration unchanged. No relaxation. The neutral electron
count follows from the 64-atom POSCAR (NELECT is not set; CP loop starts from neutral and
sets the final charge). Two static points: TARGETMU = -4.9071 and -5.1071. Nothing else.

Purpose: a same-numerics baseline for  d(sigma)_Step16 / d(sigma)_flat16  and
d(Gamma_-)(s) / d(Gamma_-)_flat16  -- the current flat reference (T, 4x4 cell, 3x3x1,
PREC=Accurate) differs in cell, k-mesh and PREC and is not fit for few-percent comparisons.

Text-based POSCAR filtering (keeps coordinate strings and selective-dynamics flags
byte-identical); ASE is used only for validation.

Usage (from Au_Cl/):  scripts/pyrun.sh scripts/build_flat16x1_pair.py
"""
import json
import os
import shutil

import numpy as np
from ase.io import read
from scipy.spatial import cKDTree

SRC = "03_pilot/Step-16x1_muref_fastcfg"
Z_CUT = 13.5            # everything above this is the added strip (top layer at 14.603)
POINTS = [("Flat16x1_muref", "-4.9071", "Step-16x1_muref_fastcfg"),
          ("Flat16x1_dUp02", "-5.1071", "Step-16x1_dUp02_fastcfg")]

# ---- POSCAR: drop the strip atoms, keep everything else verbatim ----
lines = open(f"{SRC}/POSCAR").read().splitlines()
assert lines[7].strip().lower().startswith("selective") and lines[8].strip().lower().startswith("cartesian")
header, coords = lines[:9], lines[9:]
coords = [l for l in coords if l.split()]
assert len(coords) == int(lines[6].split()[0]) == 72
keep = [l for l in coords if float(l.split()[2]) < Z_CUT]
drop = [l for l in coords if float(l.split()[2]) >= Z_CUT]
assert len(keep) == 64 and len(drop) == 8, (len(keep), len(drop))
new_header = header[:]
new_header[6] = "  64"
new_poscar = "\n".join(new_header + keep) + "\n"

# ---- validation with ASE ----
os.makedirs("03_pilot/_tmp_flat_check", exist_ok=True)
open("03_pilot/_tmp_flat_check/POSCAR", "w").write(new_poscar)
flat = read("03_pilot/_tmp_flat_check/POSCAR")
step = read(f"{SRC}/POSCAR")
shutil.rmtree("03_pilot/_tmp_flat_check")
assert len(flat) == 64
assert np.allclose(flat.get_cell(), step.get_cell(), atol=1e-10), "cell changed"
base = step[step.get_positions()[:, 2] < Z_CUT]
assert len(base) == 64
# bijective match flat <-> step base atoms, unrounded coordinates
d, idx = cKDTree(base.get_positions()).query(flat.get_positions())
assert len(set(idx)) == 64 and d.max() < 1e-8, f"positions changed: max {d.max():.2e}"
# flags preserved?
flags_src = {l.split()[0] + l.split()[1] + l.split()[2]: l.split()[3:6] for l in coords}
for l in keep:
    p = l.split()
    assert flags_src[p[0] + p[1] + p[2]] == p[3:6]
z = np.round(flat.get_positions()[:, 2], 3)
layers = sorted(set(z))
# minimum distance incl. periodic images (cell is 2.94 A along a2, so images matter)
from ase.neighborlist import neighbor_list
i, j, dd = neighbor_list("ijd", flat, 3.5)
dmin = float(dd.min())
val = dict(n_atoms=64, layers_z=layers, n_per_layer=[int((z == l).sum()) for l in layers],
           cell=flat.get_cell().tolist(), metal_top_z=float(z.max()),
           min_Au_Au_distance_periodic=dmin, removed_atoms_z=sorted({round(float(l.split()[2]), 3) for l in drop}),
           base_atoms_match_step16x1_max_displacement=float(d.max()),
           selective_dynamics_flags_preserved=True, source=SRC)
print(json.dumps(val, indent=1))
assert layers == [5.0, 7.401, 9.802, 12.202] or len(layers) == 4, layers
assert all(n == 16 for n in val["n_per_layer"])
assert abs(dmin - 2.940) < 0.01

# ---- write run dirs ----
for tag, mu, src_run in POINTS:
    d = f"03_pilot/{tag}"
    os.makedirs(d, exist_ok=True)
    open(f"{d}/POSCAR", "w").write(new_poscar)
    shutil.copy(f"{src_run}/KPOINTS".replace(src_run, f"03_pilot/{src_run}"), f"{d}/KPOINTS")
    shutil.copy(f"03_pilot/{src_run}/POTCAR", f"{d}/POTCAR")
    inc = open(f"03_pilot/{src_run}/INCAR").read().splitlines()
    # replace the header comment and SYSTEM; every numerical tag stays byte-identical
    body = [l for l in inc if not l.startswith("#")]
    body = [f"SYSTEM = {tag}" if l.startswith("SYSTEM") else l for l in body]
    assert any(l.replace(" ", "").startswith(f"TARGETMU={mu}") for l in body), "TARGETMU mismatch with source run"
    head = [f"# Flat 64-atom reference slab in the Step-16x1 cell: the 8 strip atoms of {src_run}",
            f"# removed, base atoms and cell untouched. TARGETMU={mu}. All numerical tags, k-points,",
            "# window (SOL_Z0/SOL_Z1), DIPOL and run configuration identical to the Step-16x1",
            "# production runs. Static, no relaxation. Baseline for d(sigma) and d(Gamma_-) comparisons."]
    open(f"{d}/INCAR", "w").write("\n".join(head + body) + "\n")
    jr = open(f"03_pilot/{src_run}/job-run").read().replace(f"--job-name={src_run}", f"--job-name={tag}")
    assert f"--job-name={tag}" in jr
    open(f"{d}/job-run", "w").write(jr)
    print(f"wrote {d}: POSCAR(64) INCAR(TARGETMU={mu}) KPOINTS POTCAR job-run")

# geometry library copy + validation record
open("03_pilot/all_defect_structures/Flat-16x1.poscar", "w").write(new_poscar)
json.dump(val, open("03_pilot/report_assets/batch1/flat16x1_build_validation.json", "w"), indent=1)
print("library copy: 03_pilot/all_defect_structures/Flat-16x1.poscar")
