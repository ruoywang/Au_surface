#!/usr/bin/env python3
"""Six multi-layer parent surfaces (two per class) for the size/complexity test, on the same 32 x 32 x 4 Au(111) base
as rough_sampling_v1 (a0 4.158, bottom two layers fixed). Every upper-level site comes from ONE seven-layer fcc(111)
build of the same cell, so level L sits in the true fcc hollows of level L-1 (ABC registration by construction), and
every added atom is checked to have its three supporting atoms below (no 2-D shape copied to a height without support).
Height is added above the terrace; the four-layer base is never thinned and no pit is deeper than one layer.

  M1  tiered islands      a: two tiers (wide lower island, smaller upper island)   b: three tiers
  M2  step bunching       two level-4/level-5 step edges that approach and, over part of their length, coincide as a
                          double step; elsewhere a 1-3 row terrace separates them (a, b: different widths/offsets)
  M3  open valley         two ridges (levels 4+5, b: one ridge gets a level-6 cap) separated by a valley whose floor
                          is the complete original terrace; a: narrow (~3 rows, 7.6 A), b: wide (~6 rows, 15 A)

Output: size_complex_test_v1/parents/<id>.poscar/.extxyz (arrays fixed, layer), parents_manifest.json,
        parents_summary.md, _parents_overview.png. The MD that follows uses rough_sampling_v1/md_driver.py with
        --parents-dir/--md-dir pointing here (same potential, protocol and frame rule; nothing re-validated).
Usage (from Au_Cl/):  rough_sampling_v1/pyrun_rs.sh size_complex_test_v1/build_complex_parents.py
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from ase import Atoms  # noqa: E402
from ase.build import fcc111  # noqa: E402
from ase.io import write  # noqa: E402
from matplotlib.patches import Circle  # noqa: E402

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
sys.path.insert(0, f"{ROOT}/rough_sampling_v1")
from build_parents import A0, NN, N1, N2, NLAYERS, FIXED_LAYERS, VACUUM, neighbours_ij, grow_patches, check  # noqa: E402

OUT = f"{ROOT}/size_complex_test_v1/parents"
D111 = A0 / np.sqrt(3)
N_UPPER = 3                       # levels 4, 5, 6 above the four-layer base
COL = {0: "#4a4a4a", 1: "#7a7a7a", 2: "#b9b9b9", 3: "#f2c14e", 4: "#e0603a", 5: "#8b1a1a", 6: "#4b0082"}


# ------------------------------------------------------------------------------------------------- lattice
def base_and_levels():
    """Base slab (layers 0-3) and, for each upper level L = 4..6, its 32 x 32 site grid (xy, z) from one 7-layer build."""
    seven = fcc111("Au", size=(N1, N2, NLAYERS + N_UPPER), a=A0, vacuum=VACUUM, orthogonal=False, periodic=True)
    P = seven.get_positions(); lay = np.round((P[:, 2] - P[:, 2].min()) / D111).astype(int); cell = seven.get_cell().array
    base = seven[lay < NLAYERS]; base_lay = lay[lay < NLAYERS]
    C = np.array([cell[0][:2], cell[1][:2]])
    levels = {}
    for L in range(NLAYERS, NLAYERS + N_UPPER):
        xy = P[lay == L][:, :2]; z = float(P[lay == L][0, 2])
        f = (xy - xy[0]) @ np.linalg.inv(C); ij = np.round(f * np.array([N1, N2])).astype(int) % np.array([N1, N2])
        assert len({tuple(x) for x in ij}) == N1 * N2
        grid_xy = np.zeros((N1, N2, 2)); grid_xy[ij[:, 0], ij[:, 1]] = xy
        levels[L] = dict(xy=grid_xy, z=z)
    # support map: site (i, j) of level L rests on the three level-(L-1) sites nearest in-plane (NN/sqrt(3) away)
    support = {}
    for L in range(NLAYERS + 1, NLAYERS + N_UPPER):
        lo = levels[L - 1]["xy"].reshape(-1, 2); hi = levels[L]["xy"].reshape(-1, 2)
        sup = np.zeros((N1 * N2, 3, 2), int)
        for k, p in enumerate(hi):
            best = np.full(len(lo), np.inf)
            for si in (-1, 0, 1):
                for sj in (-1, 0, 1):
                    best = np.minimum(best, np.linalg.norm(lo + si * cell[0][:2] + sj * cell[1][:2] - p, axis=1))
            idx = np.argsort(best)[:3]; assert np.all(np.abs(best[idx] - NN / np.sqrt(3)) < 0.05), "upper site not on an fcc hollow"
            sup[k] = np.array([divmod(int(q), N2) for q in idx])
        support[L] = sup.reshape(N1, N2, 3, 2)
    # the first upper level rests on the base top layer, which is complete: always supported
    return base, base_lay, levels, support, cell


def enforce_support(occ_by_level, support):
    """Remove upper-level sites whose three supports are not all occupied; iterate until stable. Returns the count removed."""
    removed = 0
    for L in sorted(occ_by_level):
        if L == NLAYERS: continue
        changed = True
        while changed:
            changed = False
            for i in range(N1):
                for j in range(N2):
                    if occ_by_level[L][i, j] and not all(occ_by_level[L - 1][a, b] for a, b in support[L][i, j]):
                        occ_by_level[L][i, j] = False; removed += 1; changed = True
    return removed


def grow_inside(rng, allowed, target, n_seeds):
    """grow_patches restricted to `allowed` sites (upper tiers grow only where they are supported)."""
    occ = np.zeros((N1, N2), bool); cand = list(zip(*np.where(allowed)))
    if not cand: return occ
    for _ in range(n_seeds):
        i, j = cand[rng.integers(0, len(cand))]; occ[i, j] = True
    while occ.sum() < min(target, allowed.sum()):
        front = {}
        for i in range(N1):
            for j in range(N2):
                if not occ[i, j]: continue
                for a, b in neighbours_ij(i, j):
                    if not occ[a, b] and allowed[a, b]: front[(a, b)] = front.get((a, b), 0) + 1
        if not front: break
        keys = list(front); w = np.array([front[k] ** 1.6 for k in keys], float); w /= w.sum()
        a, b = keys[rng.choice(len(keys), p=w)]; occ[a, b] = True
    return occ


def supported_sites(occ_below, support_L):
    s = np.zeros((N1, N2), bool)
    for i in range(N1):
        for j in range(N2):
            s[i, j] = all(occ_below[a, b] for a, b in support_L[i, j])
    return s


def strip(i0, i1, roughen=None, rng=None):
    occ = np.zeros((N1, N2), bool)
    for i in range(i0, i1): occ[i % N1, :] = True
    if roughen and rng is not None:
        for _ in range(roughen):
            for e_out, e_in in (((i0 - 1) % N1, i0 % N1), (i1 % N1, (i1 - 1) % N1)):
                j = rng.integers(0, N2); run = rng.integers(1, 4)
                for k in range(run):
                    jj = (j + k) % N2
                    if rng.random() < 0.5: occ[e_out, jj] = True
                    else: occ[e_in, jj] = False
    return occ


# ------------------------------------------------------------------------------------------------- classes
def M1(rng, variant, support):
    occ = {4: grow_patches(rng, int(0.24 * N1 * N2), 1)}
    occ[5] = grow_inside(rng, supported_sites(occ[4], support[5]), int((0.09 if variant == "a" else 0.11) * N1 * N2), 1)
    if variant == "b":
        occ[6] = grow_inside(rng, supported_sites(occ[5], support[6]), int(0.03 * N1 * N2), 1)
    return occ, dict(description="tiered island: wide lower tier, smaller upper tier(s)", tiers=2 if variant == "a" else 3)


def M2(rng, variant, support):
    """Lower strip rows [4, 20); upper strip on top with its low-side edge offset d(j) rows from the lower edge:
    d = 0 over part of the length (double step), 1-3 elsewhere (narrow intermediate terrace)."""
    lower = strip(4, 20, roughen=2, rng=rng)
    offsets = [0, 0, 1, 2, 3, 2, 1, 0] if variant == "a" else [0, 2, 4, 4, 2, 0, 0, 1]
    upper = np.zeros((N1, N2), bool); hi = 14 if variant == "a" else 12
    for j in range(N2):
        d = offsets[(j * len(offsets)) // N2]
        for i in range(4 + d, hi): upper[i, j] = True
    occ = {4: lower, 5: upper}
    return occ, dict(description="step bunching: two step edges locally merged into a double step", offsets_rows=offsets)


def M3(rng, variant, support):
    gap = 3 if variant == "a" else 6
    r1 = strip(2, 12, roughen=1, rng=rng); r2 = strip(12 + gap, 22 + gap, roughen=1, rng=rng)
    occ = {4: r1 | r2, 5: strip(4, 10) | strip(14 + gap, 20 + gap)}
    if variant == "b": occ[6] = strip(6, 8)
    return occ, dict(description="open valley between two multi-level ridges; floor = complete original terrace", valley_rows=gap, valley_width_A=round(gap * NN * np.sqrt(3) / 2, 1))


def assemble(base, base_lay, levels, occ):
    at = base.copy(); lay = base_lay.copy(); added = {}
    for L in sorted(occ):
        sel = np.argwhere(occ[L]); xy = levels[L]["xy"][sel[:, 0], sel[:, 1]]
        extra = Atoms("Au" * len(sel), positions=np.c_[xy, np.full(len(sel), levels[L]["z"])], cell=at.get_cell(), pbc=at.pbc)
        at = at + extra; lay = np.r_[lay, np.full(len(sel), L)]; added[L] = int(len(sel))
    fixed = lay < FIXED_LAYERS
    at.set_array("fixed", fixed.astype(int)); at.set_array("layer", lay.astype(int))
    return at, added


def overview(manifest, atoms_by_id):
    fig, axes = plt.subplots(2, 3, figsize=(27, 16)); plt.rcParams.update({"font.size": 14})
    for ax, m in zip(axes.ravel(), manifest):
        at = atoms_by_id[m["parent_id"]]; P = at.get_positions(); lay = at.get_array("layer"); c = at.get_cell().array
        for i in np.argsort(P[:, 2]):
            if lay[i] >= 2: ax.add_patch(Circle(P[i, :2], 0.48 * NN, fc=COL.get(int(lay[i]), "#000"), ec="k", lw=0.2))
        ax.set_xlim(-2, c[0][0] + c[1][0] + 2); ax.set_ylim(-2, c[1][1] + 2); ax.set_aspect("equal"); ax.set_axis_off()
        ax.set_title(f"{m['parent_id']}: {m['construction']['description']}\n{m['n_atoms']} Au; levels " +
                     ", ".join(f"L{L}: {n}" for L, n in m["added_per_level"].items()), fontsize=14)
    handles = [plt.Line2D([], [], marker="o", ls="", ms=12, mfc=COL[k], mec="k", label=f"layer {k}") for k in range(2, 7)]
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=14, frameon=False); fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(f"{OUT}/_parents_overview.png", dpi=90); plt.close(fig)


def main():
    os.makedirs(OUT, exist_ok=True)
    base, base_lay, levels, support, cell = base_and_levels()
    manifest = []; atoms_by_id = {}
    for cls, fn in (("M1", M1), ("M2", M2), ("M3", M3)):
        for variant in ("a", "b"):
            rng = np.random.default_rng(7000 + int(cls[1]) * 10 + ord(variant))
            occ, info = fn(rng, variant, support)
            n_unsupported = enforce_support(occ, support)
            at, added = assemble(base, base_lay, levels, occ)
            dmin = check(at); pid = f"P{cls}{variant}_s11"
            zt = at.get_positions()[:, 2].max() - at.get_positions()[:, 2].min()
            row = dict(parent_id=pid, cls=cls, level=variant, seed=11, parameters=dict(variant=variant), construction=dict(info, unsupported_sites_removed=n_unsupported),
                       n_atoms=len(at), n_base=len(base), n_fixed=int((at.get_array("layer") < FIXED_LAYERS).sum()), added_per_level={str(k): v for k, v in added.items()},
                       height_above_bottom_A=round(float(zt), 3), height_in_dft_cell_A=round(5.0 + float(zt), 3), min_AuAu_A=dmin, cell_A=cell.tolist(),
                       status="PASS" if dmin >= 2.6 else "FAIL")
            manifest.append(row); atoms_by_id[pid] = at
            write(f"{OUT}/{pid}.poscar", at, format="vasp", direct=False, sort=False); write(f"{OUT}/{pid}.extxyz", at, format="extxyz")
            print(f"  {pid}: {len(at)} Au, levels {added}, unsupported removed {n_unsupported}, top at {5.0 + zt:.2f} A in the DFT cell, min d {dmin:.3f} {row['status']}")
    json.dump(manifest, open(f"{OUT}/parents_manifest.json", "w"), indent=1)
    L = ["# Multi-layer parents for size_complex_test_v1", "", f"{len(manifest)} parents on the 32 x 32 x 4 base; every upper level from one 7-layer fcc build; "
         "every upper atom supported by its three fcc hollows (unsupported sites removed and counted).", "",
         "| id | class | description | atoms | L4 | L5 | L6 | top in DFT cell (A) | min Au-Au (A) | status |", "|---|---|---|---|---|---|---|---|---|---|"]
    for r in manifest:
        a = r["added_per_level"]; L.append(f"| {r['parent_id']} | {r['cls']} | {r['construction']['description']} | {r['n_atoms']} | {a.get('4',0)} | {a.get('5',0)} | {a.get('6',0)} | "
                                           f"{r['height_in_dft_cell_A']} | {r['min_AuAu_A']:.3f} | {r['status']} |")
    L += ["", "The production solvent window requires the highest Au below SOL_Z1 - 15 = 19.603 A (metal bottom at 5.0 A): three-level parents top at "
          f"{max(r['height_in_dft_cell_A'] for r in manifest):.2f} A before thermal motion, so cells cut from them will need the vertical-box check (step 4 of the README)."]
    open(f"{OUT}/parents_summary.md", "w").write("\n".join(L) + "\n"); overview(manifest, atoms_by_id)
    print(f"{len(manifest)} parents -> {OUT}")


if __name__ == "__main__":
    main()
