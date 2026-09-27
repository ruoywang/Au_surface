#!/usr/bin/env python3
"""R2-stripe-wall (dataset_plan_v1 rev 2): constrained straight-domain-wall / compressed-stacking approximation
of the Au(111) stripe reconstruction, 2026-09-27.

Cell: the Flat-16x1 cell (a1 = 47.04 A along the [1-10] nearest-neighbour direction, a2 = 2.94 A). Layers:
3 base layers (48 atoms, z = 5.0/7.4/9.8 A, fcc) + a top layer of 17 atoms on the 16-site row (6 % uniaxial
compression along a1) = 65 atoms, 4 layers in total (the top layer REPLACES the 4th fcc layer).

Construction (1-D soliton / Frenkel-Kontorova picture): along the row the cumulative displacement of atom i is
    U(i) = h_x [ S(i; c1, w) + S(i; c2, w) ],   V(i) = h_y [ S(i; c1, w) - S(i; c2, w) ],   S = (1 + tanh((i - c)/w)) / 2
where h = (h_x, h_y) is the fcc->hcp registry vector with |h| = a/sqrt6 = 1.70 A chosen with h_x < 0 (its component
along the row is -a/(2 sqrt2) = -1.47 A, perpendicular +-0.85 A). Wall 1 (centre c1) takes the layer from fcc to
hcp registry, wall 2 (centre c2) back to fcc; the two walls together absorb exactly one site (2 x 1.47 = 2.94 A),
which is what makes 17 atoms close periodically on 16 sites. w = 2.3 atoms keeps every neighbour distance >= 2.6 A.

Acceptance (all checked, recorded in report_assets/batch1/R2_build_validation.json):
  * 65 atoms; min Au-Au distance (in-plane periodic, z open) >= 2.6 A;
  * fcc-domain atoms within 0.3 A (xy, minimum image) of an fcc site (= layer-1 xy); hcp-domain atoms within 0.3 A
    of an hcp site (= layer-2 xy); wall atoms are the remaining ones and are listed;
  * periodic closure: the spacing between the last top atom and the image of the first equals the mean spacing.
This is NOT a model of the 22 x sqrt3 herringbone; it is the straight-wall, single-period approximation the plan
asks for, to be relaxed with the bottom two layers fixed. ny = 1 keeps the walls straight.
Usage (from Au_Cl/):  scripts/pyrun.sh scripts/build_R2_stripe_wall.py
"""
import json

import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.neighborlist import neighbor_list

A0 = 4.158
NN = A0 / np.sqrt(2)
LIB = "03_pilot/all_defect_structures"
N_TOP, N_SITES = 17, 16
W, C1, C2 = 2.3, 4.25, 12.75


def mic_xy(d, cell):
    C = cell[:2, :2]; f = np.linalg.solve(C.T, d[:2]); f -= np.round(f)
    best = None
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            v = (f + np.array([i, j])) @ C
            if best is None or np.linalg.norm(v) < np.linalg.norm(best): best = v
    return best


def S(i, c, w):
    return 0.5 * (1 + np.tanh((i - c) / w))


flat = read(f"{LIB}/Flat-16x1.poscar")
pos = flat.get_positions(); cell = flat.get_cell().array
lz = sorted(set(np.round(pos[:, 2], 3)))
assert len(lz) == 4 and abs(cell[0][1]) < 1e-6, (lz, cell[0])
l1 = pos[np.abs(pos[:, 2] - lz[0]) < 0.3]; l2 = pos[np.abs(pos[:, 2] - lz[1]) < 0.3]; l4 = pos[np.abs(pos[:, 2] - lz[3]) < 0.3]
base = flat[pos[:, 2] < lz[3] - 0.3]
assert len(base) == 48 and len(l4) == 16
# fcc sites of the top layer = the removed 4th-layer positions, ordered along a1 (x); hcp registry vector from layer 2
fcc = l4[np.argsort(l4[:, 0])]
z_top = float(lz[3])
cands = [mic_xy(q - fcc[0], cell) for q in l2]
cands = [v for v in cands if abs(np.linalg.norm(v) - A0 / np.sqrt(6)) < 0.02]
h = min([v for v in cands if v[0] < -1.0], key=lambda v: v[0])        # the hcp vector pointing backwards along the row
assert abs(h[0] + NN / 2) < 0.02 and abs(abs(h[1]) - A0 / (2 * np.sqrt(6))) < 0.02, h
# 17 atoms: atom i starts from fcc site i (site 16 = image of site 0) and is displaced by (U_i, V_i)
top = []
for i in range(N_TOP):
    site = fcc[i % N_SITES, :2] + (i // N_SITES) * cell[0][:2]
    U = h[0] * (S(i, C1, W) + S(i, C2, W)); V = h[1] * (S(i, C1, W) - S(i, C2, W))
    top.append([site[0] + U, site[1] + V, z_top])
top = np.array(top)
R2 = Atoms("Au%d" % (48 + N_TOP), positions=np.vstack([base.get_positions(), top]), cell=cell, pbc=True)
R2.wrap(pbc=(True, True, False), eps=1e-8)

# ---- validation ----
a = R2.copy(); a.set_pbc((True, True, False))
i_, j_, d_ = neighbor_list("ijd", a, 3.5)
dmin = float(d_.min())
p = R2.get_positions(); topidx = np.where(np.abs(p[:, 2] - z_top) < 0.3)[0]
def reg(q, ref):
    return min(np.linalg.norm(mic_xy(q - r, cell)) for r in ref)
labels = []
for k in topidx:
    dfcc, dhcp = reg(p[k], l1), reg(p[k], l2)
    labels.append("fcc" if dfcc < 0.3 else "hcp" if dhcp < 0.3 else "wall")
# closure: spacing between consecutive top atoms along x incl. the wrap
xs = np.sort(top[:, 0]); gaps = np.diff(np.append(xs, xs[0] + cell[0][0]))
val = dict(n_atoms=len(R2), layers_z=[float(v) for v in lz[:3]] + [z_top], top_atoms=N_TOP, sites=N_SITES, compression=1 - N_SITES / N_TOP,
           hcp_vector_A=[float(x) for x in h], wall_width_atoms=W, wall_centres=[C1, C2], min_AuAu_A=dmin,
           registry_labels_along_row=labels, n_fcc=labels.count("fcc"), n_hcp=labels.count("hcp"), n_wall=labels.count("wall"),
           row_gaps_A=[round(float(g), 3) for g in gaps], gap_min_A=float(gaps.min()), gap_max_A=float(gaps.max()),
           closure_total_length_A=float(gaps.sum()), cell_a1_A=float(cell[0][0]),
           validation="PASS" if (len(R2) == 65 and dmin >= 2.6 and labels.count("fcc") >= 3 and labels.count("hcp") >= 3
                                 and abs(gaps.sum() - cell[0][0]) < 1e-6) else "FAIL")
print(json.dumps({k: v for k, v in val.items() if k != "row_gaps_A"}, indent=1))
print("row gaps:", val["row_gaps_A"])
assert val["validation"] == "PASS", "R2 failed acceptance"
write(f"{LIB}/R2-stripe-wall.poscar", R2, format="vasp", direct=False, sort=True)
json.dump(val, open("03_pilot/report_assets/batch1/R2_build_validation.json", "w"), indent=1)
print(f"wrote {LIB}/R2-stripe-wall.poscar and 03_pilot/report_assets/batch1/R2_build_validation.json")
