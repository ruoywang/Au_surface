#!/usr/bin/env python3
"""Batch-B structures of dataset_plan_v1 rev 2 not yet in the library (2026-09-27):

  Kink-edge1 / Kink-edge2 : Step-8x3 strip (fcc111 8x3x4 = 96 + 12 strip atoms) + one extra strip atom outward of edge1 / edge2
                            in ny-row 0 -> 109 atoms. ny = 3 is the minimum period for one kink per period.
  Island-7-elongated      : 6x6 slab (144) + 7 atoms as a 3+4 two-row zigzag chain in fcc hollows (B sites) -> 151
  Pit-7-trench            : 6x6 slab with 7 top-layer atoms removed in the same 3+4 zigzag footprint -> 137
  C1-island-near-step     : library Step-8x4 (144) + compact 7-atom island on the lower terrace -> 151. A 10 A terrace has 4 rows
                            and a hexagon spans 3, so the island is necessarily attached to one step foot row; it is placed on
                            the edge2 side (rows foot+0..+2) and classified by actual connectivity.
  C2-island+pit           : the plan's 6x6 cell cannot hold two 3-row features with a 2-row rim gap (3+2+3 > 6 rows), so C2 is
                            built in the 8x8 cell (256 atoms net): island rows 0-2 and pit rows 4-6, i.e. ONE vacant row between
                            the rims on both sides (3+1+3+1 = 8) plus a 4-site a2 offset; rim-to-rim distance recorded (6.12 A).
                            Recorded as a plan adjustment (§6 rule a).

Site conventions (verified on the existing strips): fcc hollows of the top layer = xy of layer-2 atoms ("B sites");
top-layer atoms = xy of layer-1 atoms. All builds: a = 4.158 A, layers at 5.0/7.4/9.8/12.2 A, added atoms at 14.603 A.
Validation per structure (atom count, min Au-Au distance with in-plane periodicity, registry of added atoms = 3 neighbours
in the layer below at a/sqrt2, connectivity, periodic-image separation) -> report_assets/batch1/batchB_build_validation.json
Usage (from Au_Cl/):  scripts/pyrun.sh scripts/build_batchB_structures.py
"""
import json
import sys

import numpy as np
from ase import Atoms
from ase.build import fcc111
from ase.io import read, write
from ase.neighborlist import neighbor_list

sys.path.insert(0, "scripts")
from geom_check import periodic_image_isolation  # noqa: E402

A0 = 4.158
NN = A0 / np.sqrt(2)
D111 = A0 / np.sqrt(3)
LIB = "03_pilot/all_defect_structures"
rec = {}


def nl(atoms, cutoff):
    a = atoms.copy(); a.set_pbc((True, True, False)); return neighbor_list("ijd", a, cutoff)


def min_dist(atoms):
    i, j, d = nl(atoms, 3.5); return float(d.min())


def base_slab(nx, ny):
    s = fcc111("Au", size=(nx, ny, 4), a=A0, vacuum=None, orthogonal=False, periodic=True)
    p = s.get_positions(); p[:, 2] += 5.0 - p[:, 2].min(); s.set_positions(p)
    c = s.get_cell().array.copy(); c[2] = [0, 0, 44.603]; s.set_cell(c, scale_atoms=False)
    return s


def layer_idx(atoms, k):                 # k = 0 bottom ... 3 top
    z = atoms.get_positions()[:, 2]; lz = sorted(set(np.round(z, 2)))
    return np.where(np.abs(z - lz[k]) < 0.3)[0]


def frac_rows(atoms, idx):
    """row index along a1 (0..nx-1) and site index along a2 (0..ny-1) for atoms idx, from fractional coordinates."""
    f = atoms.get_scaled_positions(wrap=True)[idx]
    nx = len(np.unique(np.round(f[:, 0], 4))); ny = len(np.unique(np.round(f[:, 1], 4)))
    u = np.round(f[:, 0] * nx - f[:, 0].min() * nx).astype(int) % nx
    v = np.round(f[:, 1] * ny - f[:, 1].min() * ny).astype(int) % ny
    return u, v, nx, ny


def add_atoms(slab, xy_list, z=14.603):
    pos = np.vstack([slab.get_positions(), [[x, y, z] for x, y in xy_list]])
    a = Atoms("Au%d" % len(pos), positions=pos, cell=slab.get_cell(), pbc=True); a.wrap(pbc=(True, True, False), eps=1e-8)
    return a


def registry_ok(atoms, idx_added, z_below):
    """each added atom has exactly 3 neighbours at NN distance in the layer below (threefold hollow)."""
    i, j, d = nl(atoms, 3.1); p = atoms.get_positions(); out = {}
    for k in idx_added:
        below = [dd for ii, jj, dd in zip(i, j, d) if ii == k and abs(p[jj, 2] - z_below) < 0.3]
        inplane = [dd for ii, jj, dd in zip(i, j, d) if ii == k and abs(p[jj, 2] - p[k, 2]) < 0.3]
        out[int(k)] = dict(n_below=len(below), n_inplane=len(inplane), ok=len(below) == 3 and all(abs(x - NN) < 0.02 for x in below))
    return out


# ---------------- Kinks ----------------
for name, side in (("Kink-edge1", "edge1"), ("Kink-edge2", "edge2")):
    s = base_slab(8, 3)
    l2 = layer_idx(s, 1); p = s.get_positions()
    u, v, nx, ny = frac_rows(s, l2)
    assert nx == 8 and ny == 3
    strip_rows = [0, 1, 2, 3]
    strip_xy = [p[k, :2] for k, uu in zip(l2, u) if uu in strip_rows]
    kink_row = (strip_rows[0] - 1) % nx if side == "edge1" else (strip_rows[-1] + 1) % nx
    kink_xy = [p[k, :2] for k, uu, vv in zip(l2, u, v) if uu == kink_row and vv == 0]
    assert len(strip_xy) == 12 and len(kink_xy) == 1
    at = add_atoms(s, strip_xy + kink_xy)
    added = list(range(96, 109))
    reg = registry_ok(at, added, 12.2)
    kink = 108
    r = dict(n_atoms=len(at), min_AuAu_A=min_dist(at), strip_rows=strip_rows, kink_row=int(kink_row), kink_side=side,
             kink_atom=dict(n_below=reg[kink]["n_below"], n_inplane_strip=reg[kink]["n_inplane"], CN=reg[kink]["n_below"] + reg[kink]["n_inplane"]),
             edge_atoms_inplane_CN=sorted(set(reg[k]["n_inplane"] for k in added[:-1])),
             all_added_in_hollows=all(x["ok"] for x in reg.values()),
             edge_contour="straight edge of 3 atoms per period with one atom protruding one row outward at v=0: two corners per period on this edge",
             cell_A=[round(float(np.linalg.norm(at.cell[0])), 2), round(float(np.linalg.norm(at.cell[1])), 2)])
    r["validation"] = "PASS" if (r["n_atoms"] == 109 and r["min_AuAu_A"] > 2.9 and r["all_added_in_hollows"] and r["kink_atom"]["CN"] == 5) else "FAIL"
    rec[name] = r; print(name, json.dumps(r))
    if r["validation"] == "PASS": write(f"{LIB}/{name}.poscar", at, format="vasp", direct=False, sort=True)

# ---------------- Island-7-elongated / Pit-7-trench (6x6) ----------------
s66 = base_slab(6, 6)
l2 = layer_idx(s66, 1); l4 = layer_idx(s66, 3); p = s66.get_positions()
u2, v2, nx, ny = frac_rows(s66, l2)


def zigzag_footprint(idx, u, v, rows=(2, 3), vs=(1, 2, 3, 4)):
    """4 sites in row rows[0] at v in vs, plus the 3 sites of row rows[1] having >= 2 NN among them."""
    four = [k for k, uu, vv in zip(idx, u, v) if uu == rows[0] and vv in vs]
    cand = [k for k, uu in zip(idx, u) if uu == rows[1]]
    a = s66.copy(); a.set_pbc((True, True, False))
    P = a.get_positions(); C = a.get_cell().array
    def d(k, m):
        dd = P[m] - P[k]; f = np.linalg.solve(C[:2, :2].T, dd[:2]); f -= np.round(f)
        return min(np.linalg.norm((f + [i, j]) @ C[:2, :2]) for i in (-1, 0, 1) for j in (-1, 0, 1))
    three = [k for k in cand if sum(abs(d(k, m) - NN) < 0.02 for m in four) >= 2]
    assert len(four) == 4 and len(three) == 3, (len(four), len(three))
    return four + three


isl_idx = zigzag_footprint(l2, u2, v2)
isl = add_atoms(s66, [p[k, :2] for k in isl_idx])
added = list(range(144, 151)); reg = registry_ok(isl, added, 12.2)
xy = isl.get_positions()[added][:, :2]; cell = isl.get_cell().array
iso = periodic_image_isolation(xy, (cell[0][:2], cell[1][:2]))
r = dict(n_atoms=len(isl), min_AuAu_A=min_dist(isl), island_inplane_CN=[reg[k]["n_inplane"] for k in added],
         all_in_hollows=all(x["ok"] for x in reg.values()), connected=all(reg[k]["n_inplane"] >= 2 for k in added),
         periodic_image_separation_A=iso["defect_image_min_distance"], shape="3+4 two-row zigzag chain")
r["validation"] = "PASS" if (r["n_atoms"] == 151 and r["min_AuAu_A"] > 2.9 and r["all_in_hollows"] and r["connected"] and r["periodic_image_separation_A"] >= 8.0) else "FAIL"
rec["Island-7-elongated"] = r; print("Island-7-elongated", json.dumps(r))
if r["validation"] == "PASS": write(f"{LIB}/Island-7-elongated.poscar", isl, format="vasp", direct=False, sort=True)

u4, v4, _, _ = frac_rows(s66, l4)
pit_idx = zigzag_footprint(l4, u4, v4)
keep = [k for k in range(len(s66)) if k not in pit_idx]
pit = s66[keep]
vac_xy = p[pit_idx][:, :2]
iso = periodic_image_isolation(vac_xy, (cell[0][:2], cell[1][:2]))
# footprint connectivity: each vacancy has >= 2 NN vacancies
def nn_count(xy_all):
    out = []
    for a_ in xy_all:
        n = 0
        for b_ in xy_all:
            dd = b_ - a_; f = np.linalg.solve(cell[:2, :2].T, dd); f -= np.round(f)
            dist = min(np.linalg.norm((f + [i, j]) @ cell[:2, :2]) for i in (-1, 0, 1) for j in (-1, 0, 1))
            if abs(dist - NN) < 0.02: n += 1
        out.append(n)
    return out
r = dict(n_atoms=len(pit), n_removed=len(pit_idx), min_AuAu_A=min_dist(pit), vacancy_nn_counts=nn_count(vac_xy),
         connected=all(n >= 2 for n in nn_count(vac_xy)), periodic_image_separation_A=iso["defect_image_min_distance"], shape="3+4 two-row zigzag trench")
r["validation"] = "PASS" if (r["n_atoms"] == 137 and r["connected"] and r["periodic_image_separation_A"] >= 8.0) else "FAIL"
rec["Pit-7-trench"] = r; print("Pit-7-trench", json.dumps(r))
if r["validation"] == "PASS": write(f"{LIB}/Pit-7-trench.poscar", pit, format="vasp", direct=False, sort=True)


def hexagon(idx, u, v, uc, vc, atoms):
    """centre site (uc, vc) and its 6 in-plane NN among idx."""
    P = atoms.get_positions(); C = atoms.get_cell().array
    centre = [k for k, uu, vv in zip(idx, u, v) if uu == uc and vv == vc]; assert len(centre) == 1; c = centre[0]
    def d(k):
        dd = P[k] - P[c]; f = np.linalg.solve(C[:2, :2].T, dd[:2]); f -= np.round(f)
        return min(np.linalg.norm((f + [i, j]) @ C[:2, :2]) for i in (-1, 0, 1) for j in (-1, 0, 1))
    ring = [k for k in idx if k != c and abs(d(k) - NN) < 0.02]; assert len(ring) == 6, len(ring)
    return [c] + ring


# ---------------- C1: Step-8x4 + compact island on the lower terrace ----------------
s84 = read(f"{LIB}/Step-8x4.poscar")
p = s84.get_positions(); cell = s84.get_cell().array
strip = np.where(p[:, 2] > 13.5)[0]; assert len(strip) == 16
l2 = layer_idx(s84, 1); u2, v2, nx, ny = frac_rows(s84, l2); assert nx == 8 and ny == 4
f = s84.get_scaled_positions(wrap=True)
strip_rows = sorted(set(np.round(f[strip, 0] * nx - f[l2, 0].min() * nx).astype(int) % nx))
edge2_row = max(strip_rows); foot_row = (edge2_row + 1) % nx
hex_idx = hexagon(l2, u2, v2, (foot_row + 1) % nx, 1, s84)          # centre one row beyond the foot row -> hexagon rows foot..foot+2
c1 = add_atoms(s84, [p[k, :2] for k in hex_idx])
added = list(range(144, 151)); reg = registry_ok(c1, added, 12.2)
P = c1.get_positions(); i_, j_, d_ = nl(c1, 3.1)
strip_new = [k for k in range(len(c1)) if k < 144 and P[k, 2] > 13.5]
touch = sorted(set(int(jj) for ii, jj, dd in zip(i_, j_, d_) if ii in added and jj in strip_new and abs(dd - NN) < 0.02))
xy = P[added][:, :2]; iso = periodic_image_isolation(xy, (cell[0][:2], cell[1][:2]))
r = dict(n_atoms=len(c1), min_AuAu_A=min_dist(c1), strip_rows=[int(x) for x in strip_rows], island_rows=[int(foot_row), int((foot_row + 1) % nx), int((foot_row + 2) % nx)],
         all_in_hollows=all(x["ok"] for x in reg.values()), island_strip_contacts=len(touch),
         classification="step-attached (island foot row = edge2 foot row)" if touch else "step-detached",
         periodic_image_separation_A=iso["defect_image_min_distance"],
         note="10 A lower terrace has 4 rows, a hexagon spans 3: attachment to one foot row is unavoidable in Step-8x4; placed on the edge2 side")
r["validation"] = "PASS" if (r["n_atoms"] == 151 and r["min_AuAu_A"] > 2.9 and r["all_in_hollows"]) else "FAIL"
rec["C1-island-near-step"] = r; print("C1-island-near-step", json.dumps(r))
if r["validation"] == "PASS": write(f"{LIB}/C1-island-near-step.poscar", c1, format="vasp", direct=False, sort=True)

# ---------------- C2: island + pit, 8x8 cell ----------------
s88 = base_slab(8, 8); p = s88.get_positions(); cell = s88.get_cell().array
l2 = layer_idx(s88, 1); l4 = layer_idx(s88, 3)
u2, v2, nx, ny = frac_rows(s88, l2); u4, v4, _, _ = frac_rows(s88, l4)
isl_idx = hexagon(l2, u2, v2, 1, 1, s88)              # island rows 0,1,2
pit_idx = hexagon(l4, u4, v4, 5, 5, s88)              # pit rows 4,5,6: one vacant row on both sides of the 8-row cell; a2 offset 4 sites
keep = [k for k in range(len(s88)) if k not in pit_idx]
c2 = add_atoms(s88[keep], [p[k, :2] for k in isl_idx])
added = list(range(len(keep), len(keep) + 7)); reg = registry_ok(c2, added, 12.2)
P = c2.get_positions()
def mind(A_, B_):
    C = cell[:2, :2]; best = 9e9
    for a_ in A_:
        for b_ in B_:
            dd = b_ - a_; ff = np.linalg.solve(C.T, dd); ff -= np.round(ff)
            best = min(best, min(np.linalg.norm((ff + [i, j]) @ C) for i in (-1, 0, 1) for j in (-1, 0, 1)))
    return float(best)
sep_rim = mind(P[added][:, :2], p[pit_idx][:, :2])
iso_i = periodic_image_isolation(P[added][:, :2], (cell[0][:2], cell[1][:2])); iso_p = periodic_image_isolation(p[pit_idx][:, :2], (cell[0][:2], cell[1][:2]))
r = dict(n_atoms=len(c2), min_AuAu_A=min_dist(c2), island_rows=[0, 1, 2], pit_rows=[4, 5, 6], all_in_hollows=all(x["ok"] for x in reg.values()),
         island_rim_to_pit_rim_min_distance_A=sep_rim, island_image_separation_A=iso_i["defect_image_min_distance"], pit_image_separation_A=iso_p["defect_image_min_distance"],
         note="plan adjustment: 6x6 cannot hold two 3-row features with a 2-row rim gap (3+2+3 > 6); built in 8x8 (256 atoms net) with ONE vacant row between the rims on both sides (3+1+3+1 = 8) and a 4-site a2 offset; rim-to-rim distance recorded")
r["validation"] = "PASS" if (r["n_atoms"] == 256 and r["min_AuAu_A"] > 2.9 and r["all_in_hollows"] and sep_rim > 1.7 * NN) else "FAIL"
rec["C2-island+pit"] = r; print("C2-island+pit", json.dumps(r))
if r["validation"] == "PASS": write(f"{LIB}/C2-island+pit.poscar", c2, format="vasp", direct=False, sort=True)

json.dump(rec, open("03_pilot/report_assets/batch1/batchB_build_validation.json", "w"), indent=1)
print("summary:", {k: v["validation"] for k, v in rec.items()})
