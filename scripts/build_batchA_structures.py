#!/usr/bin/env python3
"""Batch-A structures of dataset_plan_v1 rev 2 that were not yet in the library (2026-09-27):

  Flat-8x2                              : Step-8x2 with its 8 strip atoms removed (64)      -- same recipe as Flat-16x1
  Step-8x2_edge-vacancy_plus_foot-adatom: Step-8x2 with one edge2 strip atom moved to the fcc hollow of the
                                          lower terrace adjacent to edge2 (72; atom number and ORDER preserved,
                                          so it is the end point of the detachment path with identity mapping)
  R1-hcp-terminated                     : T-4x4 with the whole top layer translated by the fcc->hcp registry
                                          vector taken from the slab itself (|d| = a/sqrt6 = 1.70 A)

Every build is validated (atom count, min distance incl. periodic images, registry conditions) and the
record written to 03_pilot/report_assets/batch1/batchA_build_validation.json.  R2 is built separately.
Usage (from Au_Cl/):  scripts/pyrun.sh scripts/build_batchA_structures.py
"""
import json

import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.neighborlist import neighbor_list

A0 = 4.158
LIB = "03_pilot/all_defect_structures"
rec = {}


def mic_xy(d, cell):
    """minimum-image in-plane vector for a Cartesian difference d (2D cell rows)."""
    C = cell[:2, :2]
    f = np.linalg.solve(C.T, d[:2])
    f -= np.round(f)
    return f @ C


def nl(atoms, cutoff):
    """neighbour list with z-periodicity OFF: the library templates have a short c vector (vacuum=None),
    so z-images would produce fake 2.4 A contacts between the bottom and top layers."""
    a = atoms.copy(); a.set_pbc((True, True, False))
    return neighbor_list("ijd", a, cutoff)


def min_dist(atoms):
    i, j, d = nl(atoms, 3.5)
    return float(d.min()) if len(d) else float("inf")


def layers(atoms, tol=0.3):
    z = atoms.get_positions()[:, 2]
    lv = []
    for v in np.sort(z):
        if not lv or abs(v - lv[-1]) > tol: lv.append(v)
    return lv


# ---------- Flat-8x2 ----------
s82 = read(f"{LIB}/Step-8x2.poscar")
z = s82.get_positions()[:, 2]
flat = s82[z < 13.5]
assert len(flat) == 64 and len(s82) == 72
write(f"{LIB}/Flat-8x2.poscar", flat, format="vasp", direct=False, sort=True)
rec["Flat-8x2"] = dict(n_atoms=len(flat), min_AuAu_A=min_dist(flat), layers_z=[round(v, 3) for v in layers(flat)], source="Step-8x2.poscar minus strip")
print("Flat-8x2:", rec["Flat-8x2"])

# ---------- Step-8x2_edge-vacancy_plus_foot-adatom ----------
at = s82.copy()
pos = at.get_positions(); cell = at.get_cell().array
strip = np.where(pos[:, 2] > 13.5)[0]
assert len(strip) == 8
# fractional coordinate along a1 (across-step); edge2 = strip atoms with the largest u (the u ~ 0.42-0.46 side)
frac = at.get_scaled_positions(wrap=True)
u = frac[strip, 0]
rows_u = np.sort(np.unique(np.round(u, 4)))
assert len(rows_u) == 4, rows_u
edge2 = strip[np.isclose(np.round(u, 4), rows_u[-1])]           # 2 atoms (ny = 2)
edge1 = strip[np.isclose(np.round(u, 4), rows_u[0])]
assert len(edge2) == 2 and len(edge1) == 2
# row translation vector: from the row rows_u[-2] to rows_u[-1] within the same ny index (nearest pair)
row_prev = strip[np.isclose(np.round(u, 4), rows_u[-2])]
mover = edge2[0]
d_candidates = [pos[mover] - pos[k] for k in row_prev]
d_candidates = [np.array([*mic_xy(d, cell), 0.0]) for d in d_candidates]
row_vec = min(d_candidates, key=np.linalg.norm)
assert abs(np.linalg.norm(row_vec) - A0 / 2) < 0.05 or abs(np.linalg.norm(row_vec) - A0 / np.sqrt(2)) < 0.05, np.linalg.norm(row_vec)
# the foot site: continue the strip lattice one row beyond edge2 (same B-registry lattice => fcc hollow of the lower terrace)
new_xy = pos[mover] + row_vec
new_xy[2] = pos[mover, 2]
pos2 = pos.copy(); pos2[mover] = new_xy
at.set_positions(pos2); at.wrap(pbc=(True, True, False), eps=1e-8)
# validation: same atom count/order, min distance, the moved atom sits in a hollow of layer 4 (3 neighbours at 2.94 in layer 4)
i, j, d = nl(at, 3.1)
nn_mover = [(jj, dd) for ii, jj, dd in zip(i, j, d) if ii == mover]
below = [dd for jj, dd in nn_mover if abs(at.get_positions()[jj, 2] - 12.2) < 0.3]
assert len(below) == 3 and all(abs(x - A0 / np.sqrt(2)) < 0.02 for x in below), (len(below), below)
assert len(at) == 72 and min_dist(at) > 2.9
write(f"{LIB}/Step-8x2_edge-vacancy_plus_foot-adatom.poscar", at, format="vasp", direct=False, sort=False)
rec["Step-8x2_edge-vacancy_plus_foot-adatom"] = dict(n_atoms=72, moved_atom_index=int(mover), from_row_u=float(rows_u[-1]),
        displacement_A=[float(x) for x in row_vec], min_AuAu_A=min_dist(at), foot_site_neighbours_in_layer4=len(below),
        note="atom order identical to Step-8x2.poscar (sort=False) so the 3 detachment images interpolate with identity mapping; "
             "Step-8x2.poscar itself is written sorted -- the path generator must use this file's order for both ends")
print("edge-vacancy+foot-adatom:", rec["Step-8x2_edge-vacancy_plus_foot-adatom"])
# also store an order-matched copy of the ideal Step-8x2 as the path start (same atom order as the end point)
write(f"{LIB}/Step-8x2_pathstart_ordered.poscar", s82, format="vasp", direct=False, sort=False)

# ---------- R1-hcp-terminated ----------
T = read(f"{LIB}/T.poscar")
pos = T.get_positions(); cell = T.get_cell().array
lz = layers(T); assert len(lz) == 4
top = np.where(np.abs(pos[:, 2] - lz[3]) < 0.3)[0]
l2 = np.where(np.abs(pos[:, 2] - lz[1]) < 0.3)[0]
l3 = np.where(np.abs(pos[:, 2] - lz[2]) < 0.3)[0]
# fcc->hcp registry vector: from the first top atom to the nearest layer-2 xy (minimum image); must be a/sqrt6
cands = [np.array([*mic_xy(pos[k] - pos[top[0]], cell), 0.0]) for k in l2]
shift = min(cands, key=np.linalg.norm)
assert abs(np.linalg.norm(shift) - A0 / np.sqrt(6)) < 0.02, np.linalg.norm(shift)
posR = pos.copy(); posR[top] += shift
R1 = Atoms("Au64", positions=posR, cell=cell, pbc=True); R1.wrap(pbc=(True, True, False), eps=1e-8)
# acceptance: every top atom above a layer-2 atom (xy within 0.02 A, MIC) AND in a threefold hollow of layer 3 (3 neighbours at a/sqrt2)
pR = R1.get_positions()
for t in top:
    dxy = min(np.linalg.norm(mic_xy(pR[t] - pR[k], cell)) for k in l2)
    assert dxy < 0.02, dxy
i, j, d = nl(R1, 3.1)
for t in top:
    nb = [dd for ii, jj, dd in zip(i, j, d) if ii == t and abs(pR[jj, 2] - lz[2]) < 0.3]
    assert len(nb) == 3 and all(abs(x - A0 / np.sqrt(2)) < 0.02 for x in nb), (t, nb)
    atop = [dd for ii, jj, dd in zip(i, j, d) if ii == t and abs(pR[jj, 2] - lz[2]) < 0.3 and dd < 2.5]
    assert not atop
write(f"{LIB}/R1-hcp-terminated.poscar", R1, format="vasp", direct=False, sort=True)
rec["R1-hcp-terminated"] = dict(n_atoms=64, shift_vector_A=[float(x) for x in shift], shift_length_A=float(np.linalg.norm(shift)),
        expected_a_over_sqrt6=A0 / np.sqrt(6), min_AuAu_A=min_dist(R1),
        acceptance="every top atom within 0.02 A (xy, MIC) of a layer-2 atom AND 3 layer-3 neighbours at a/sqrt2; no atop position")
print("R1-hcp-terminated:", rec["R1-hcp-terminated"])

json.dump(rec, open("03_pilot/report_assets/batch1/batchA_build_validation.json", "w"), indent=1)
print("wrote 03_pilot/report_assets/batch1/batchA_build_validation.json")
