#!/usr/bin/env python3
"""Batch 2: Island-7/Island-19 and Pit-7/Pit-19, all in a COMMON 8x8 base
cell, so the 7- vs 19-atom size comparison isn't confounded by also
changing the cell (defect_plan.md Correction 3 / the follow-up review's
explicit requirement: put both sizes in the same, big-enough cell).
Compact shapes only this round; elongated deferred to a later pass."""
import os
import numpy as np
from ase.build import fcc111
from ase.io import write
from ase import Atom

A0 = 4.158
NLAYERS = 4
NX, NY = 8, 8
OUT = "03_pilot/all_defect_structures"
os.makedirs(OUT, exist_ok=True)


def base_slab(nx, ny, nlayers=NLAYERS):
    slab = fcc111('Au', size=(nx, ny, nlayers), a=A0, vacuum=None, orthogonal=False, periodic=True)
    pos = slab.get_positions()
    pos[:, 2] += (5.0 - pos[:, 2].min())
    slab.set_positions(pos)
    return slab


def nearest_n_sites(base, idx_pool, center_idx, n):
    d = base.get_distances(center_idx, list(idx_pool), mic=True)
    order = np.argsort(d)
    pool = np.array(list(idx_pool))
    return pool[order[:n]]


# ---- shared registry/candidate setup ----
base = base_slab(NX, NY)
posb = base.get_positions()
zb = np.sort(np.unique(np.round(posb[:, 2], 2)))
d111 = zb[-1] - zb[-2]
top_z = zb[-1]

registry_z = zb[NLAYERS - 3]  # fcc-continuation registry, same convention as elsewhere
registry_idx = np.where(np.isclose(posb[:, 2], registry_z, atol=0.05))[0]
top_idx = np.where(np.isclose(posb[:, 2], top_z, atol=0.05))[0]

# pick a center near the cell centroid, in BOTH the registry set and the top set
cand_xy = posb[registry_idx][:, :2]
centroid = cand_xy.mean(axis=0)
d2 = ((cand_xy - centroid) ** 2).sum(axis=1)
center_registry_idx = registry_idx[np.argsort(d2)[0]]

top_xy = posb[top_idx][:, :2]
d2t = ((top_xy - centroid) ** 2).sum(axis=1)
center_top_idx = top_idx[np.argsort(d2t)[0]]

for n_atoms, tag in [(7, "Island-7-8x8"), (19, "Island-19-8x8")]:
    s = base_slab(NX, NY)
    chosen = nearest_n_sites(base, registry_idx, center_registry_idx, n_atoms)
    for i in chosen:
        x, y, _ = posb[i]
        s.append(Atom('Au', position=[x, y, top_z + d111]))
    write(f"{OUT}/{tag}.poscar", s, format="vasp", direct=False, sort=True)
    print(tag, len(s), "atoms —", len(chosen), "island atoms")

for n_atoms, tag in [(7, "Pit-7-8x8"), (19, "Pit-19-8x8")]:
    s = base_slab(NX, NY)
    chosen = nearest_n_sites(base, top_idx, center_top_idx, n_atoms)
    for i in sorted(chosen, reverse=True):
        del s[i]
    write(f"{OUT}/{tag}.poscar", s, format="vasp", direct=False, sort=True)
    print(tag, len(s), "atoms —", len(chosen), "removed")
