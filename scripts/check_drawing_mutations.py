"""Mutation test: reintroduce each bug and require the regression to FAIL.

A check that never fails proves nothing, so each of the four drawing bugs is put back one at a time and the
regression must exit non-zero.
"""
import sys, io, contextlib, numpy as np
sys.path.insert(0, "/anvil/scratch/x-rywang/Au_Cl/scripts")
import render_gallery as rg
import check_gallery_drawing as chk

orig = dict(layer_spacing=rg.layer_spacing, tiles=rg.tiles, cut_profile=rg.cut_profile, cut_band=rg.cut_band)

def run():
    """(exit code, first failure line) of the regression as it stands."""
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf): chk.main()
        code = 0
    except SystemExit as e:
        code = e.code or 0
    lines = buf.getvalue().splitlines()
    try:
        first = lines[lines.index("FAILURES:") + 1].strip()
    except ValueError:
        first = ""
    return code, first

# --- bug 1: nearest neighbour taken in projection, mixing layers
def ls_projected(atoms, sel=None):
    pos = atoms.get_positions(); cell = atoms.get_cell().array
    p = pos[rg.top_atoms(atoms) if sel is None else sel]
    if len(p) < 2: return 2.94
    best = np.full(len(p), np.inf)
    for si in (-1,0,1):
        for sj in (-1,0,1):
            sh = (si*cell[0]+sj*cell[1])[:2]
            d = np.linalg.norm(p[:,None,:2]-(p[None,:,:2]+sh),axis=-1)
            d = np.where(d > 0.1, d, np.inf)
            best = np.minimum(best, d.min(axis=1))
    best = best[np.isfinite(best)]
    return float(np.median(best)) if len(best) else 2.94

# --- bug 2: each tiled copy re-evaluated at a shifted query (only +-1 images searched)
def tiles_reevaluated(gx, gy, A, cell, i_range, j_range):
    for i in i_range:
        for j in j_range:
            sh = (i*cell[0]+j*cell[1])[:2]
            if i == min(i_range) and j == min(j_range):
                yield gx+sh[0], gy+sh[1], A
            else:
                B = A.copy()
                if abs(i) > 1 or abs(j) > 1: B[:] = A.max()      # the far copy loses its features
                yield gx+sh[0], gy+sh[1], B

# --- bug 3: the profile rolled while A stays put
def profile_rolled(F, ci, cell):
    s, pr, L = orig["cut_profile"](F, ci, cell)
    return s, np.roll(pr, len(pr)//3), L

# --- bug 4: the two coordinates folded independently
def band_split_images(pos, cell, ci, half):
    rA, t, n, L, P = rg.cut_frame(cell, ci)
    d = (pos[:, :2] - rA) @ n
    dv = ((d / P) + 0.5) % 1.0 - 0.5
    s = ((pos[:, :2] - rA) @ t) % L
    return np.abs(dv) * P < half, s, dv * P

cases = [("1 projected nearest neighbour (disc radius too small)", "layer_spacing", ls_projected),
         ("2 tiled copy re-evaluated, far copy blank",             "tiles",         tiles_reevaluated),
         ("3 profile rolled but A-A' not moved",                   "cut_profile",   profile_rolled),
         ("4 perpendicular and along-line folded independently",   "cut_band",      band_split_images)]

code, _ = run()
print(f"baseline (no mutation): exit {code}  -> {'PASS' if code==0 else 'UNEXPECTED FAILURE'}")
allgood = code == 0
for name, attr, fn in cases:
    setattr(rg, attr, fn)
    if attr in ("tiles","cut_band","cut_profile","layer_spacing"): setattr(chk, attr, fn) if hasattr(chk, attr) else None
    code, first = run()
    ok = code != 0
    allgood &= ok
    print(f"bug {name}: exit {code}  -> {'CAUGHT' if ok else 'NOT CAUGHT'}")
    if ok: print(f"      first complaint: {first[:130]}")
    setattr(rg, attr, orig[attr])
    if hasattr(chk, attr): setattr(chk, attr, orig[attr])
print()
print("mutation test: " + ("every reintroduced bug is caught" if allgood else "SOME BUGS SLIP THROUGH"))
sys.exit(0 if allgood else 1)
