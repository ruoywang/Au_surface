#!/usr/bin/env python3
"""Batch 1 new geometry: Au(211), Au(221) true crystallographic vicinal
slabs, and a properly-sized Island-7 on a 6x6 base (per defect_plan.md
Correction 3). All geometry-only, no DFT input prep (single-sided vacuum
window handled separately if/when these become DFT candidates)."""
import os
import numpy as np
from ase.build import bulk, surface, fcc211, fcc111
from ase.io import write
from ase import Atom

A0 = 4.158
NLAYERS = 4
OUT = "03_pilot/all_defect_structures"
os.makedirs(OUT, exist_ok=True)


def pad_z(atoms, vacuum=15.0):
    """ase.build.surface()/fcc211(vacuum=None) leave the cell's c-vector
    degenerate (zero length) -- fine for these builders' own bookkeeping,
    but unwritable and unusable for any real distance/PBC check. Pad it."""
    pos = atoms.get_positions()
    cell = atoms.get_cell()
    zmin, zmax = pos[:, 2].min(), pos[:, 2].max()
    cell[2] = [0, 0, (zmax - zmin) + 2 * vacuum]
    atoms.set_cell(cell)
    pos[:, 2] += vacuum - zmin
    atoms.set_positions(pos)
    return atoms

# --- Au(211): A-type step ({100} microfacet), 3(111)x(100) unit ---
au211 = fcc211('Au', size=(3, 4, 4), a=A0, vacuum=None, orthogonal=True)
au211 = pad_z(au211)
write(f"{OUT}/Au211.poscar", au211, format="vasp", direct=False, sort=True)
print("Au211", len(au211), "atoms, thickness",
      au211.get_positions()[:, 2].max() - au211.get_positions()[:, 2].min(), "A")

# --- Au(221): B-type step ({111} microfacet) ---
au_bulk = bulk('Au', 'fcc', a=A0)
au221 = surface(au_bulk, (2, 2, 1), layers=8, vacuum=None, periodic=True)
au221 = au221 * (3, 3, 1)
au221 = pad_z(au221)
write(f"{OUT}/Au221.poscar", au221, format="vasp", direct=False, sort=True)
print("Au221", len(au221), "atoms, thickness",
      au221.get_positions()[:, 2].max() - au221.get_positions()[:, 2].min(), "A")


# --- Island-7 on 6x6 (corrected base cell per Correction 3) ---
def base_slab(nx, ny, nlayers=NLAYERS):
    slab = fcc111('Au', size=(nx, ny, nlayers), a=A0, vacuum=None, orthogonal=False, periodic=True)
    pos = slab.get_positions()
    pos[:, 2] += (5.0 - pos[:, 2].min())
    slab.set_positions(pos)
    return slab


base = base_slab(6, 6)
posb = base.get_positions()
zb = np.sort(np.unique(np.round(posb[:, 2], 2)))
d111 = zb[-1] - zb[-2]

# fcc-continuation registry site set, read from THIS build (same fix as Step script)
registry_z = zb[NLAYERS - 3]
cand_idx = np.where(np.isclose(posb[:, 2], registry_z, atol=0.05))[0]
cand_xy = posb[cand_idx][:, :2]

centroid = cand_xy.mean(axis=0)
d2 = ((cand_xy - centroid) ** 2).sum(axis=1)
c_local = np.argsort(d2)[0]  # index into cand_idx/cand_xy, nearest the cell centroid
center_atom_idx = cand_idx[c_local]

d = base.get_distances(center_atom_idx, list(cand_idx), mic=True)
order = np.argsort(d)
chosen = cand_idx[order[:7]]  # center + 6 nearest = compact 7-site island

for i in chosen:
    x, y, _ = posb[i]
    base.append(Atom('Au', position=[x, y, zb[-1] + d111]))

write(f"{OUT}/Island-7-6x6.poscar", base, format="vasp", direct=False, sort=True)
print("Island-7-6x6", len(base), "atoms —", len(chosen), "island atoms added")
