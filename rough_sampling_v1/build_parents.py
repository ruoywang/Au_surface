#!/usr/bin/env python3
"""Build the 32 large parent surfaces for rough_sampling_v1: four classes x two parameter levels x four seeds.

Common base: a continuous Au(111) slab, 32 x 32 in-plane, FOUR layers (4096 Au), a0 = 4.158 A, bottom two layers
fixed. Every upper-layer site that any construction uses comes from ONE five-layer fcc build of the same cell,
so islands and strips sit in the correct fcc hollows of the layer below by construction; no two slabs of
different thickness are ever stitched together.

Classes (parameters are coverage targets for candidate generation, not an experimental equilibrium):
  A  curved steps      strips of a fifth layer whose edges are locally advanced / retreated to make kinks,
                       bends and locally narrow terraces; mean terrace ~10 A or ~20 A
  B  island growth     several islands grown site by site on fifth-layer sites to 0.15 or 0.35 ML, so that at
                       the higher coverage islands touch, neck and merge
  C  pits / retreat    connected patches of the TOP layer removed, 0.10 or 0.25 ML, never deeper than one layer
  D  mixed             a stepped surface from which fifth-layer atoms are transferred onto the lower terrace next
                       to the step: islands in contact with steps, pits behind them, atom count conserved

Each parent is written as POSCAR (selective dynamics: bottom two layers F) and .extxyz, with a manifest row
recording class, level, seed, counts, coverage and every site added or removed. The final morphology is meant
to be enriched by the MD that follows; nothing here is claimed to be an equilibrium surface.

Usage (from Au_Cl/):  scripts/pyrun.sh rough_sampling_v1/build_parents.py
Output: rough_sampling_v1/parents/<id>.poscar, <id>.extxyz, parents_manifest.json, parents_summary.md
"""
import json
import os

import numpy as np
from ase import Atoms
from ase.build import fcc111
from ase.io import write

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
OUT = f"{ROOT}/rough_sampling_v1/parents"
A0 = 4.158
NN = A0 / np.sqrt(2)                   # 2.940 A
N1 = N2 = 32
NLAYERS = 4
FIXED_LAYERS = 2
VACUUM = 16.0                          # A each side; the production template sets the real cell later
LEVELS = {"A": {"1": dict(terrace_A=10.0), "2": dict(terrace_A=20.0)},
          "B": {"1": dict(coverage=0.15), "2": dict(coverage=0.35)},
          "C": {"1": dict(coverage=0.10), "2": dict(coverage=0.25)},
          "D": {"1": dict(transfer=0.08), "2": dict(transfer=0.16)}}
SEEDS = [11, 23, 37, 51]


# ------------------------------------------------------------------------------------------------- lattice
def base_and_sites():
    """The 4-layer slab plus the fifth-layer site lattice, from ONE 5-layer fcc(111) build."""
    five = fcc111("Au", size=(N1, N2, NLAYERS + 1), a=A0, vacuum=VACUUM, orthogonal=False, periodic=True)
    z = five.get_positions()[:, 2]
    layers = np.round((z - z.min()) / (A0 / np.sqrt(3)))
    top_sites = five.get_positions()[layers == NLAYERS]
    base = five[layers < NLAYERS]
    # layer index per atom, 0 = bottom
    zb = base.get_positions()[:, 2]
    lay = np.round((zb - zb.min()) / (A0 / np.sqrt(3))).astype(int)
    return base, lay, top_sites[:, :2], float(top_sites[0, 2]), five.get_cell().array


def frac_grid(xy, cell, origin_xy):
    """Integer (i, j) lattice indices of in-plane positions on the 32 x 32 site grid, measured from the first
    atom of the SAME layer. Each fcc(111) layer is offset from the next by a third of a primitive cell, which
    is not a multiple of 1/32 of the supercell, so rounding absolute fractions would split one layer's atoms
    between two integer indices; measuring from an atom of that layer makes every index exact."""
    C = np.array([cell[0][:2], cell[1][:2]])
    f = (xy - origin_xy) @ np.linalg.inv(C)
    return np.round(f * np.array([N1, N2])).astype(int) % np.array([N1, N2])


def neighbours_ij(i, j):
    """The six in-plane neighbours of site (i, j) on a triangular lattice indexed along a1, a2."""
    return [((i + 1) % N1, j), ((i - 1) % N1, j), (i, (j + 1) % N2), (i, (j - 1) % N2),
            ((i + 1) % N1, (j - 1) % N2), ((i - 1) % N1, (j + 1) % N2)]


# ------------------------------------------------------------------------------------------------- classes
def class_A(rng, terrace_A, cell):
    """Strips of width ~terrace along a1 (period 2 x terrace), edges roughened by local advance / retreat."""
    L1 = np.linalg.norm(cell[0][:2]); row = L1 / N1
    period_rows = max(4, int(round(2 * terrace_A / row)))
    n_strips = max(1, int(round(N1 / period_rows))); period_rows = N1 // n_strips
    half = period_rows // 2
    occ = np.zeros((N1, N2), bool)
    for s in range(n_strips):
        i0 = s * period_rows
        for i in range(i0, i0 + half): occ[i % N1, :] = True
    # roughen each edge: random runs of 1-4 sites pushed out or pulled in by one row, several passes
    for _ in range(3):
        for s in range(n_strips):
            i0 = s * period_rows; edges = [(i0 - 1) % N1, (i0 + half) % N1]   # row just outside each edge
            inner = [i0 % N1, (i0 + half - 1) % N1]                             # row just inside each edge
            for e_out, e_in in zip(edges, inner):
                j = rng.integers(0, N2); run = rng.integers(1, 5)
                for k in range(run):
                    jj = (j + k) % N2
                    if rng.random() < 0.5: occ[e_out, jj] = True        # advance
                    else: occ[e_in, jj] = False                          # retreat
    return occ, dict(n_strips=n_strips, period_rows=period_rows)


def grow_patches(rng, target_sites, n_seeds, attach_bias=1.6):
    """Site-by-site growth of several patches to a total of target_sites, with attachment probability rising
    with the number of occupied neighbours, which makes compact islands that touch and merge when crowded."""
    occ = np.zeros((N1, N2), bool)
    for _ in range(n_seeds):
        while True:
            i, j = rng.integers(0, N1), rng.integers(0, N2)
            if not occ[i, j]: occ[i, j] = True; break
    while occ.sum() < target_sites:
        front = {}
        for i in range(N1):
            for j in range(N2):
                if not occ[i, j]: continue
                for a, b in neighbours_ij(i, j):
                    if not occ[a, b]: front[(a, b)] = front.get((a, b), 0) + 1
        keys = list(front); w = np.array([front[k] ** attach_bias for k in keys], float); w /= w.sum()
        a, b = keys[rng.choice(len(keys), p=w)]
        occ[a, b] = True
    return occ


def class_B(rng, coverage, cell):
    n_seeds = 6 if coverage < 0.25 else 9
    occ = grow_patches(rng, int(round(coverage * N1 * N2)), n_seeds)
    return occ, dict(n_seeds=n_seeds)


def class_C(rng, coverage, cell):
    n_seeds = 5 if coverage < 0.2 else 8
    occ = grow_patches(rng, int(round(coverage * N1 * N2)), n_seeds)
    return occ, dict(n_seeds=n_seeds)              # occ marks sites REMOVED from the top layer


def class_D(rng, transfer, cell):
    """A two-strip stepped surface; `transfer` ML of fifth-layer atoms are taken from the strips (leaving pits
    and retreated edges) and placed on the lower terrace as islands touching the step foot. Atom count conserved."""
    occ, info = class_A(rng, 20.0, cell)
    n_move = int(round(transfer * N1 * N2))
    moved = 0; removed, added = [], []
    strip_sites = list(zip(*np.where(occ)))
    while moved < n_move and strip_sites:
        # remove a compact patch from inside a strip
        k = strip_sites[rng.integers(0, len(strip_sites))]
        patch = [k]
        for a, b in neighbours_ij(*k):
            if occ[a, b] and len(patch) < 4 and rng.random() < 0.7: patch.append((a, b))
        for a, b in patch:
            if occ[a, b]: occ[a, b] = False; removed.append((int(a), int(b))); moved += 1
        strip_sites = list(zip(*np.where(occ)))
        # place them on the lower terrace adjacent to a step foot (an empty site with an occupied neighbour)
        foot = [(i, j) for i in range(N1) for j in range(N2) if not occ[i, j]
                and any(occ[a, b] for a, b in neighbours_ij(i, j))]
        for _ in range(len(patch)):
            if not foot: break
            a, b = foot[rng.integers(0, len(foot))]
            if not occ[a, b]: occ[a, b] = True; added.append((int(a), int(b)))
            foot = [(i, j) for i in range(N1) for j in range(N2) if not occ[i, j]
                    and any(occ[x, y] for x, y in neighbours_ij(i, j))]
    info.update(n_transferred=moved, removed_sites=removed, added_sites=added)
    return occ, info


# ------------------------------------------------------------------------------------------------- assembly
def assemble(base, lay, site_xy, z_top, cell, cls, occ):
    """Return the Atoms object for one parent and the bookkeeping of what was added / removed."""
    at = base.copy()
    ij = frac_grid(site_xy, cell, site_xy[0])
    added, removed = [], []
    assert len({tuple(x) for x in ij}) == N1 * N2, "fifth-layer sites do not map 1:1 onto the site grid"
    if cls in ("A", "B", "D"):
        sel = [k for k in range(len(site_xy)) if occ[ij[k, 0], ij[k, 1]]]
        extra = Atoms("Au" * len(sel), positions=np.c_[site_xy[sel], np.full(len(sel), z_top)],
                      cell=at.get_cell(), pbc=at.pbc)
        at = at + extra
        added = [tuple(int(x) for x in ij[k]) for k in sel]
        lay = np.r_[lay, np.full(len(sel), NLAYERS)]
    else:                                                   # C: remove top-layer atoms at the marked sites
        top = np.flatnonzero(lay == NLAYERS - 1)
        txy = at.get_positions()[top][:, :2]
        tij = frac_grid(txy, cell, txy[0])
        assert len({tuple(x) for x in tij}) == N1 * N2, "top-layer atoms do not map 1:1 onto the site grid"
        kill = [top[k] for k in range(len(top)) if occ[tij[k, 0], tij[k, 1]]]
        removed = [tuple(int(x) for x in tij[k]) for k in range(len(top)) if occ[tij[k, 0], tij[k, 1]]]
        keep = np.ones(len(at), bool); keep[kill] = False
        at = at[keep]; lay = lay[keep]
    fixed = lay < FIXED_LAYERS
    at.set_array("fixed", fixed.astype(int))
    at.set_array("layer", lay.astype(int))
    return at, added, removed, fixed


def check(at):
    """Minimum Au-Au distance with in-plane periodicity; must be >= 2.6 A for an unrelaxed fcc construction."""
    p = at.get_positions(); cell = at.get_cell().array
    best = np.inf
    for si in (-1, 0, 1):
        for sj in (-1, 0, 1):
            sh = si * cell[0] + sj * cell[1]
            # chunked to keep memory flat for 4000+ atoms
            for a in range(0, len(p), 512):
                d = np.linalg.norm(p[a:a + 512, None, :] - (p[None, :, :] + sh), axis=-1)
                d = np.where(d > 0.1, d, np.inf); best = min(best, float(d.min()))
    return best


def main():
    os.makedirs(OUT, exist_ok=True)
    base, lay, site_xy, z_top, cell = base_and_sites()
    manifest = []
    for cls in "ABCD":
        for lvl, par in LEVELS[cls].items():
            for seed in SEEDS:
                rng = np.random.default_rng(seed * 1000 + ord(cls) * 10 + int(lvl))
                fn = {"A": class_A, "B": class_B, "C": class_C, "D": class_D}[cls]
                occ, info = fn(rng, list(par.values())[0], cell)
                at, added, removed, fixed = assemble(base, lay, site_xy, z_top, cell, cls, occ)
                dmin = check(at)
                pid = f"P{cls}{lvl}_s{seed}"
                top_n = int(((at.get_array("layer")) == NLAYERS).sum()) if cls != "C" else 0
                row = dict(parent_id=pid, cls=cls, level=lvl, seed=seed, parameters=par, construction=info,
                           n_atoms=len(at), n_base=len(base), n_fixed=int(fixed.sum()),
                           n_fifth_layer=top_n, n_top_removed=len(removed),
                           coverage_fifth_layer_ML=top_n / (N1 * N2), removed_fraction_top_ML=len(removed) / (N1 * N2),
                           min_AuAu_A=dmin, cell_A=cell.tolist(), added_sites=added, removed_sites=removed,
                           status="PASS" if dmin >= 2.6 else "FAIL")
                manifest.append(row)
                write(f"{OUT}/{pid}.poscar", at, format="vasp", direct=False, sort=False)
                write(f"{OUT}/{pid}.extxyz", at, format="extxyz")
                print(f"  {pid}: {len(at)} atoms, fifth layer {top_n} ({top_n/(N1*N2):.3f} ML), "
                      f"removed {len(removed)}, min d {dmin:.3f} A  {row['status']}")
    json.dump(manifest, open(f"{OUT}/parents_manifest.json", "w"), indent=1)
    L = ["# Parent surfaces for rough_sampling_v1", "",
         f"{len(manifest)} parents: 32 x 32 Au(111), {NLAYERS} layers ({len(base)} base Au), a0 = {A0} A, "
         f"bottom {FIXED_LAYERS} layers fixed. Fifth-layer sites come from one 5-layer fcc build of the same cell.", "",
         "| id | class | level | seed | atoms | fifth-layer ML | top removed ML | min Au-Au (A) | status |",
         "|---|---|---|---|---|---|---|---|---|"]
    for r in manifest:
        L.append(f"| {r['parent_id']} | {r['cls']} | {r['level']} | {r['seed']} | {r['n_atoms']} | "
                 f"{r['coverage_fifth_layer_ML']:.3f} | {r['removed_fraction_top_ML']:.3f} | {r['min_AuAu_A']:.3f} | {r['status']} |")
    open(f"{OUT}/parents_summary.md", "w").write("\n".join(L) + "\n")
    print(f"{len(manifest)} parents -> {OUT}/  ({sum(r['status']=='PASS' for r in manifest)} PASS)")


if __name__ == "__main__":
    main()
