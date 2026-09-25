#!/usr/bin/env python3
"""Reduced-cell Step series: same construction as build_step_structure.py
(registry-safe, single-build fcc-continuation site), but with ny=1 instead
of ny=4 along the step direction. A straight, unperturbed step has no
physical periodicity requirement beyond the primitive translation along
its own length (kinks/perturbations would need a longer period, but this
is the idealized-straight-step model) -- ny=4 was carrying x4 redundant
along-step repeats for zero physical reason.

Terrace width (nx, across-step) is UNCHANGED -- only the along-step
redundancy is cut. Validates by tiling the reduced cell back up x4 along
a2 and checking it reproduces the ny=4 build atom-for-atom.
"""
import os
import numpy as np
from ase.build import fcc111
from ase.io import write, read

A0 = 4.158
NLAYERS = 4
OUT = "03_pilot/all_defect_structures"
os.makedirs(OUT, exist_ok=True)

SIZES = {"Step-8x1": (8, 1), "Step-16x1": (16, 1), "Step-24x1": (24, 1),
          "Step-8x2": (8, 2), "Step-16x2": (16, 2), "Step-24x2": (24, 2)}


def base_slab(nx, ny, nlayers=NLAYERS):
    slab = fcc111('Au', size=(nx, ny, nlayers), a=A0, vacuum=None, orthogonal=False, periodic=True)
    pos = slab.get_positions()
    pos[:, 2] += (5.0 - pos[:, 2].min())
    slab.set_positions(pos)
    return slab


def build(nx, ny):
    base = base_slab(nx, ny)
    posb = base.get_positions()
    frac = base.get_scaled_positions()
    zb = np.sort(np.unique(np.round(posb[:, 2], 2)))
    d111 = zb[-1] - zb[-2]
    registry_z = zb[NLAYERS - 3]
    registry_idx = np.where(np.isclose(posb[:, 2], registry_z, atol=0.05))[0]
    keep = [i for i in registry_idx if frac[i, 0] < 0.5]
    from ase import Atom
    for i in keep:
        x, y, _ = posb[i]
        base.append(Atom('Au', position=[x, y, zb[-1] + d111]))
    return base, len(keep)


for tag, (nx, ny) in SIZES.items():
    s, n_added = build(nx, ny)
    write(f"{OUT}/{tag}.poscar", s, format="vasp", direct=False, sort=True)
    expected_total = int(nx * ny * (NLAYERS + 0.5))
    print(tag, len(s), "atoms (expected", expected_total, ") —", n_added, "added")

print()
print("=== validation: tile ny=1 cell x4 along a2, compare to the existing ny=4 build ===")
for nx in [8, 16, 24]:
    small, _ = build(nx, 1)
    tiled = small * (1, 4, 1)
    big = read(f"{OUT}/Step-{nx}x4.poscar")

    pos_t = np.round(tiled.get_positions(), 3)
    pos_b = np.round(big.get_positions(), 3)
    # compare as sets of (x,y,z) since atom ordering/index may differ after write(sort=True)
    set_t = set(map(tuple, pos_t))
    set_b = set(map(tuple, pos_b))
    match = (len(tiled) == len(big)) and (set_t == set_b)
    print(f"Step-{nx}x1 tiled x4  vs  Step-{nx}x4:  N_tiled={len(tiled)} N_big={len(big)}  "
          f"coordinate sets match: {match}")
    if not match:
        print("  MISMATCH -- extra in tiled:", len(set_t - set_b), " extra in big:", len(set_b - set_t))
