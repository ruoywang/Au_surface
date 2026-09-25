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

def strict_validate(tiled, big, tag, tol=1e-6):
    """Strict one-to-one atom correspondence, not a rounded-coordinate set
    comparison (which loses duplicate counts and can't distinguish 'every
    atom has a match' from 'every atom has SOME match, possibly reused').
    Uses unrounded coordinates; fails loudly (raises) on any mismatch."""
    from scipy.spatial import cKDTree

    if len(tiled) != len(big):
        raise AssertionError(f"{tag}: atom count mismatch, tiled={len(tiled)} big={len(big)}")
    if set(tiled.get_chemical_symbols()) != {"Au"} or set(big.get_chemical_symbols()) != {"Au"}:
        raise AssertionError(f"{tag}: unexpected species present")

    cell_t = tiled.get_cell()[:2, :2]
    cell_b = big.get_cell()[:2, :2]
    if not np.allclose(cell_t, cell_b, atol=1e-6):
        raise AssertionError(f"{tag}: in-plane cell vectors differ: tiled={cell_t.tolist()} big={cell_b.tolist()}")

    pt = tiled.get_positions()
    pb = big.get_positions()
    tree = cKDTree(pb)
    dist, idx = tree.query(pt, k=1)

    # bijection check: every big-atom index used at most once (no duplicate
    # coordinates silently collapsed onto the same target)
    if len(set(idx.tolist())) != len(big):
        raise AssertionError(f"{tag}: nearest-neighbor mapping is not a bijection "
                              f"({len(set(idx.tolist()))} distinct targets for {len(big)} atoms) "
                              "-- duplicate/collapsed coordinates likely")
    max_err = float(dist.max())
    if max_err > tol:
        raise AssertionError(f"{tag}: max atom displacement {max_err:.3e} A exceeds tolerance {tol:.0e} A")
    return max_err


print()
print("=== strict validation: tile ny=1 cell x4 along a2, compare to the existing ny=4 build ===")
validation_results = {}
for nx in [8, 16, 24]:
    small, _ = build(nx, 1)
    tiled = small * (1, 4, 1)
    big = read(f"{OUT}/Step-{nx}x4.poscar")
    tag = f"Step-{nx}x1 tiled x4 vs Step-{nx}x4"
    max_err = strict_validate(tiled, big, tag)
    validation_results[f"Step-{nx}x1"] = max_err
    print(f"{tag}: PASSED, max atom displacement = {max_err:.3e} A (bijective match, unrounded coordinates)")

import json
with open(f"{OUT}/../report_assets/batch1/step_reduction_validation.json", "w") as f:
    json.dump(validation_results, f, indent=2)
print("\nwrote max-displacement results to step_reduction_validation.json")
