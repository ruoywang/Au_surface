#!/usr/bin/env python3
"""Drawing regression for the structure gallery: does each figure agree with the coordinates it came from?

Five assertions, for every structure. Each one is written so that the bug it guards against FAILS it -- a check
that a helper agrees with itself proves nothing, so every test recomputes the quantity independently of the
function under test.

  1. CONNECTIVITY. The number of connected components of the disc-union outline equals the number of connected
     components of the feature atoms' bond graph (cutoff 1.15 x the shortest Au-Au distance in the slab,
     computed independently of the spacing function under test). Guards against a disc
     radius built from a cross-layer projected distance (1.6975 A) instead of the bond length (2.9401 A), which
     drew A3's trimer as three circles and Pit-7-compact's floor as twelve.

  2. PERIODIC REPLICATION. Every tile that tiles() yields carries the SAME array as the base cell, and the
     level-set coverage fraction is identical for every copy. Guards against re-evaluating the field at a
     shifted query, which reached only +-1 image and left Step-8x2's third copy 0.40% covered against 53.18%.

  3. A-A' CORRESPONDENCE. The section profile is compared against the field re-sampled independently at
     A + s*t_hat, with A reconstructed from cut_frame. Any shift applied to one and not the other fails here.

  4. BAND WIDTH, at the half-width the renderer actually uses. For each atom the side view draws, the
     perpendicular distance is recomputed from scratch by an explicit search over the (-1,0,1) images, and the
     along-line coordinate is recomputed from the SAME image that search picked. Comparing cut_band's mask
     against cut_band's own distance would pass trivially and would not see a wrong image for s.

  5. PERIODIC-IMAGE INVARIANCE. Adding or subtracting an in-plane lattice vector must leave (mask, s, d_perp)
     unchanged. Guards against folding the two coordinates independently: in a sheared cell the across vector
     has a component along the line, so with a = (10,0) and b = (5,8.6603) the same atom written as r and as
     r + b came out at s = 9.5 and s = 4.5.

Exits NON-ZERO when anything fails, so a publish step that only checks the exit status cannot ship a bad figure.

Usage (from Au_Cl/):  scripts/pyrun.sh scripts/check_gallery_drawing.py
"""
import json
import sys

import numpy as np
from scipy import ndimage
from ase.io import read

sys.path.insert(0, "/anvil/scratch/x-rywang/Au_Cl/scripts")
from render_gallery import (flatten_cell, surface_model, feature_sets, layer_spacing, cut_profile,  # noqa: E402
                            cut_band, cut_frame, parent_of, tiles, wrap_pad, _line, _runs,
                            BAND_FAMILIES, BAND_HALF_FRAC, TILE_TARGET)

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
GAL = f"{ROOT}/analysis/gallery/gallery.json"


def bond_cutoff(at):
    """1.15 x the shortest Au-Au distance in the slab, computed here and NOT from layer_spacing().

    The connectivity test must not take its bond length from the function it is testing: with a radius built
    from the broken projected spacing, an atom graph cut at 1.15 x that same broken value would split in the
    same places and the test would pass on two matching errors. The shortest interatomic distance is an
    independent property of the coordinates."""
    p = at.get_positions(); cell = at.get_cell().array
    best = np.inf
    for si in (-1, 0, 1):
        for sj in (-1, 0, 1):
            sh = (si * cell[0] + sj * cell[1])
            d = np.linalg.norm(p[:, None, :] - (p[None, :, :] + sh[None, None, :]), axis=-1)
            d = np.where(d > 0.1, d, np.inf)
            best = min(best, float(d.min()))
    return 1.15 * best


def atom_components(pts, cell, cutoff):
    """Connected components of the feature atoms' bond graph, minimum image."""
    n = len(pts)
    if n == 0: return 0
    par = list(range(n))

    def f(a):
        while par[a] != a: par[a] = par[par[a]]; a = par[a]
        return a
    C = np.array([cell[0][:2], cell[1][:2]]); Ci = np.linalg.inv(C)
    fr = (pts[:, None, :] - pts[None, :, :]) @ Ci; fr -= np.round(fr)
    d = np.linalg.norm(fr @ C, axis=-1)
    for i in range(n):
        for j in range(i + 1, n):
            if d[i, j] < cutoff: par[f(i)] = f(j)
    return len({f(i) for i in range(n)})


def disc_components(D, R):
    """Connected components of the disc union on the PERIODIC grid."""
    lab, k = ndimage.label(D < R)
    if k <= 1: return int(k)
    par = list(range(k + 1))

    def f(a):
        while par[a] != a: par[a] = par[par[a]]; a = par[a]
        return a
    for A, B in ((lab[0, :], lab[-1, :]), (lab[:, 0], lab[:, -1])):
        for u, v in zip(A, B):
            if u and v: par[f(u)] = f(v)
    return len({f(i) for i in range(1, k + 1)})


def perp_and_along(pos_xy, cell, ci):
    """(|d_perp|, s) recomputed from scratch: pick the image that minimises the perpendicular distance, then
    read BOTH coordinates off that one image. Deliberately not written in terms of cut_band."""
    rA, t, n, L, _ = cut_frame(cell, ci)
    across = cell[1 - ci["axis"]][:2]; along = cell[ci["axis"]][:2]
    best_d = np.full(len(pos_xy), np.inf); best_s = np.zeros(len(pos_xy))
    for m in (-2, -1, 0, 1, 2):
        for q in (-1, 0, 1):
            r = pos_xy - rA - m * across[None, :] - q * along[None, :]
            d = np.abs(r @ n)
            better = d < best_d
            best_d[better] = d[better]; best_s[better] = (r @ t)[better] % L
    return best_d, best_s


def main():
    gal = json.load(open(GAL))
    bad, rows = [], []
    for sid in sorted(gal):
        m = gal[sid]
        at = flatten_cell(read(f"{m['dir']}/{m['geom']}")); cell = at.get_cell().array
        pa = parent_of(sid, at)
        sm = surface_model(at, prefer_xy=(pa["added"] if pa else None))
        a0 = layer_spacing(at)
        bcut = bond_cutoff(at)
        hi, lo, _ = feature_sets(at)
        note = []

        # --- 1. connectivity, per feature
        for pts, lay in zip((hi, lo), sm["layers"] or [None, None]):
            if lay is None or not len(pts): continue
            na = atom_components(pts, cell, bcut)
            nd = disc_components(lay["D"], sm["R"])
            if na != nd:
                bad.append(f"{sid}: a feature of {len(pts)} atoms is {na} bonded cluster(s) but {nd} drawn blob(s)")
            note.append(f"{len(pts)}at/{na}cl")

        # --- 2. periodic replication, through the SAME generator the renderer uses
        t1 = max(1, min(3, int(round(TILE_TARGET / np.linalg.norm(cell[0][:2])))))
        t2 = max(1, min(3, int(round(TILE_TARGET / np.linalg.norm(cell[1][:2])))))
        field = sm["Z"] if sm["mode"] == "blob" else sm["H"]
        px, py, pF = wrap_pad(sm["gx"], sm["gy"], field, cell)
        lvl = 0.5 * 2.4 if sm["mode"] == "blob" else float(np.median(pF))
        base = None; ncopy = 0
        for x, y, A in tiles(px, py, pF, cell, range(-1, t1 + 1), range(-1, t2 + 1)):
            ncopy += 1
            if A is not pF:
                bad.append(f"{sid}: a tiled copy is not the base array translated"); break
            frac = float((A > lvl).mean())
            if base is None: base = frac
            elif abs(frac - base) > 1e-12:
                bad.append(f"{sid}: tiled copy coverage {100*frac:.2f}% differs from the first copy's "
                           f"{100*base:.2f}%")
        note.append(f"{ncopy}tiles")

        # --- 3, 4, 5. per section line
        prof_desc = ""
        half = BAND_HALF_FRAC * a0 if m.get("family") in BAND_FAMILIES else None
        for ci in sm["cuts"]:
            s, pr, L = cut_profile(field, ci, cell)
            rA, t, nv, LL, _ = cut_frame(cell, ci)
            C = np.array([cell[0][:2], cell[1][:2]]); Ci = np.linalg.inv(C)
            fr = ((rA[None, :] + s[:, None] * t[None, :]) @ Ci) % 1.0
            ii = np.clip((fr[:, 0] * field.shape[0]).astype(int), 0, field.shape[0] - 1)
            jj = np.clip((fr[:, 1] * field.shape[1]).astype(int), 0, field.shape[1] - 1)
            resid = float(np.abs(field[ii, jj] - pr).max())
            if resid > 1e-9:
                bad.append(f"{sid}: the profile does not match the marked {ci['label']}-{ci['label']}' line "
                           f"(max {resid:.3g})")
            lab = np.where(pr > 0.5, "H", np.where(pr < -0.5, "L", "T"))
            runs = []
            for c in lab:
                if runs and runs[-1][0] == c: runs[-1][1] += 1
                else: runs.append([c, 1])
            prof_desc += f" {ci['label']}:{'-'.join(f'{c}{n}' for c, n in runs)}"

            pos = at.get_positions()
            hw = half if half is not None else 1e9
            mk, ss, dp = cut_band(pos, cell, ci, hw)
            ref_d, ref_s = perp_and_along(pos[:, :2], cell, ci)
            if half is not None:
                if mk.any() and float(ref_d[mk].max()) > half + 1e-9:
                    bad.append(f"{sid}: an atom the {ci['label']} band draws is "
                               f"{float(ref_d[mk].max()):.3f} A out, past its stated {half:.3f} A")
                if (ref_d < half - 1e-9).sum() != mk.sum():
                    bad.append(f"{sid}: the {ci['label']} band holds {int(mk.sum())} atoms, independent "
                               f"recount gives {int((ref_d < half - 1e-9).sum())}")
            ds = np.abs((ss - ref_s + 0.5 * LL) % LL - 0.5 * LL)
            if mk.any() and float(ds[mk].max()) > 1e-6:
                bad.append(f"{sid}: the {ci['label']} side view places an atom {float(ds[mk].max()):.3f} A from "
                           f"where its own periodic image puts it")

            # --- 5. image invariance
            for shift in (cell[0][:2], -cell[0][:2], cell[1][:2], -cell[1][:2], cell[0][:2] + cell[1][:2]):
                q = pos.copy(); q[:, :2] += shift
                mk2, ss2, dp2 = cut_band(q, cell, ci, hw)
                d_s = np.abs((ss - ss2 + 0.5 * LL) % LL - 0.5 * LL)
                if (mk != mk2).any() or float(d_s.max()) > 1e-6 or float(np.abs(dp - dp2).max()) > 1e-6:
                    bad.append(f"{sid}: shifting every atom by an in-plane lattice vector changes the "
                               f"{ci['label']} side view (max ds {float(d_s.max()):.3f} A)")
                    break

        rows.append(f"{sid:38s} {sm['mode']:5s} cuts={len(sm['cuts'])} a0={a0:.4f} R={sm['R']:.3f} "
                    f"[{' '.join(note)}]{prof_desc}")

    print("\n".join(rows))
    print()
    if bad:
        print("FAILURES:")
        print("\n".join("  " + b for b in bad))
        raise SystemExit(1)
    print(f"all {len(gal)} structures pass: disc blobs = bonded clusters; every tiled copy is the base array "
          f"translated; profile = marked line; band width and side-view coordinates match an independent "
          f"recomputation; and the view is invariant under a lattice translation")


if __name__ == "__main__":
    main()
