#!/usr/bin/env python3
"""Build all 8 parent defect types from v2 doc Sec.3.2 (T/V1/V2/V3/V7/A1-fcc/A1-hcp/A3)
for illustration in the report. T/V1/A1-fcc reuse the already-computed geometry;
V2/V3/V7/A1-hcp/A3 are geometry-only (not yet run through DFT), 4 layers, same a0=4.158 A,
except V7 which needs the 6x6 supercell per the doc to avoid self-interaction."""
import os
import numpy as np
from ase.build import fcc111
from ase.constraints import FixAtoms
from ase.io import write
from ase import Atom

A0 = 4.158
NLAYERS = 4
BOTTOM_FIXED = 2
OUT = "03_pilot/all_defect_structures"
os.makedirs(OUT, exist_ok=True)


def base_slab(nx, ny):
    slab = fcc111('Au', size=(nx, ny, NLAYERS), a=A0, vacuum=None, orthogonal=False, periodic=True)
    pos = slab.get_positions()
    pos[:, 2] += (5.0 - pos[:, 2].min())
    slab.set_positions(pos)
    return slab


def top_layer_idx(slab, layer_z):
    pos = slab.get_positions()
    return np.where(np.isclose(pos[:, 2], layer_z[-1], atol=0.05))[0]


def nearest_center_idx(slab, idxs):
    pos = slab.get_positions()[idxs][:, :2]
    centroid = pos.mean(axis=0)
    d2 = ((pos - centroid) ** 2).sum(axis=1)
    return idxs[np.argmin(d2)]


def neighbor_order(slab, center_i, idxs):
    """Return idxs (excluding center) sorted by in-plane periodic distance to center_i."""
    d = slab.get_distances(center_i, [j for j in idxs if j != center_i], mic=True)
    order = np.argsort(d)
    cand = np.array([j for j in idxs if j != center_i])
    return cand[order], d[order]


def finalize(tag, slab, note):
    write(f"{OUT}/{tag}.poscar", slab, format="vasp", direct=False, sort=True)
    print(tag, len(slab), "atoms —", note)


# ---- T, V1, A1-fcc: reuse already-computed 4x4 geometry ----
slab_T = base_slab(4, 4)
pos = slab_T.get_positions()
layer_z = np.sort(np.unique(np.round(pos[:, 2], 2)))
fixed_mask = np.array([np.any(np.isclose(z, layer_z[:BOTTOM_FIXED], atol=0.05)) for z in pos[:, 2]])
slab_T.set_constraint(FixAtoms(mask=fixed_mask.tolist()))
finalize("T", slab_T, "flat terrace, no defect")

top_idx = top_layer_idx(slab_T, layer_z)
center = nearest_center_idx(slab_T, top_idx)
nbrs, dists = neighbor_order(slab_T, center, top_idx)

# V1: remove center
s = base_slab(4, 4)
del s[center]
finalize("V1", s, "1 vacancy")

# V2: remove center + its single nearest neighbor
s = base_slab(4, 4)
del s[sorted([center, nbrs[0]], reverse=True)[0]]
del s[sorted([center, nbrs[0]], reverse=True)[1]]
finalize("V2", s, "2 nearest-neighbor vacancies")

# V3: remove center + two neighbors that are mutually nearest-neighbors too (compact triangle)
# among the first-shell neighbors (6, all ~d_nn), pick two adjacent ones (60 deg apart)
first_shell = nbrs[:6]
pos_c = slab_T.get_positions()[center][:2]
angles = []
for j in first_shell:
    v = slab_T.get_positions()[j][:2] - pos_c
    angles.append(np.arctan2(v[1], v[0]))
order_ang = first_shell[np.argsort(angles)]
tri = [center, order_ang[0], order_ang[1]]  # adjacent pair -> compact triangle with center
s = base_slab(4, 4)
for i in sorted(tri, reverse=True):
    del s[i]
finalize("V3", s, "3 vacancies, compact triangle")

# V7: center + all 6 first-shell neighbors, on 6x6 (per doc, avoid self-interaction)
slab_T66 = base_slab(6, 6)
pos66 = slab_T66.get_positions()
layer_z66 = np.sort(np.unique(np.round(pos66[:, 2], 2)))
top_idx66 = top_layer_idx(slab_T66, layer_z66)
center66 = nearest_center_idx(slab_T66, top_idx66)
nbrs66, _ = neighbor_order(slab_T66, center66, top_idx66)
remove7 = [center66] + list(nbrs66[:6])
s = base_slab(6, 6)
for i in sorted(remove7, reverse=True):
    del s[i]
finalize("V7", s, "7 vacancies (center + 6 neighbors), 6x6 cell")

# A1-fcc / A1-hcp: adatom at fcc-continuation vs hcp-hollow registry
d111 = layer_z[-1] - layer_z[-2]
top_z = layer_z[-1]


def add_adatom(nx, ny, registry_layer_offset):
    s = base_slab(nx, ny)
    p = s.get_positions()
    frac = s.get_scaled_positions()
    lz = np.sort(np.unique(np.round(p[:, 2], 2)))
    ref_idx = np.where(np.isclose(p[:, 2], lz[-1 - registry_layer_offset], atol=0.05))[0][0]
    ref_xy = frac[ref_idx, :2]
    new_cart = s.get_cell().T @ [ref_xy[0], ref_xy[1], 0.0]
    s.append(Atom('Au', position=[new_cart[0], new_cart[1], lz[-1] + d111]))
    return s, new_cart[:2]


s, xy_fcc = add_adatom(4, 4, 2)  # registry of layer TWO below top -> continues ABC -> fcc site
finalize("A1_fcc", s, "adatom, fcc-continuation hollow")

s, xy_hcp = add_adatom(4, 4, 1)  # registry of layer ONE below top -> breaks ABC -> hcp site
finalize("A1_hcp", s, "adatom, hcp hollow")

# A3: three adjacent adatoms at fcc-continuation sites forming a compact triangular island
s = base_slab(4, 4)
p = s.get_positions()
frac = s.get_scaled_positions()
lz = np.sort(np.unique(np.round(p[:, 2], 2)))
ref_idx = np.where(np.isclose(p[:, 2], lz[-3], atol=0.05))[0]
ref_xy_all = frac[ref_idx, :2]
centroid = ref_xy_all.mean(axis=0)
d2 = ((ref_xy_all - centroid) ** 2).sum(axis=1)
c_i = np.argsort(d2)[0]
cart_all = np.array([s.get_cell().T @ [x, y, 0.0] for x, y in ref_xy_all])
dists_2d = np.linalg.norm(cart_all - cart_all[c_i], axis=1)
order = np.argsort(dists_2d)
tri_xy = cart_all[order[:3]]  # center + 2 nearest fcc-registry sites -> compact trimer
for xy in tri_xy:
    s.append(Atom('Au', position=[xy[0], xy[1], lz[-1] + d111]))
finalize("A3", s, "3-atom island, fcc-continuation sites")
