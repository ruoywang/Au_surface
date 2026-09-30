#!/usr/bin/env python3
"""Large schematic renders of every structure in the dataset, coloured by coordination number.

One PNG per structure: a tiled top view (so the periodic motif is visible) beside a side view. Atoms carry the
SAME coordination classes the spatial analysis uses to cut regions (CN <= 6 kink/adatom, 7-8 edge/rim, 9 terrace,
>= 10 sub-surface), so the gallery and the anion maps can be read against each other directly.

Geometry source: the POSCAR that was actually computed (05_production/<structure>/ideal__mu-4.9071), so the picture
is the calculation, not a re-generated idealisation. Structures whose production entry is a relaxation use the
accepted CONTCAR.

Usage (from Au_Cl/):  scripts/pyrun.sh scripts/render_gallery.py [--only NAME]
Outputs: analysis/gallery/<structure_id>.png  and  analysis/gallery/gallery.json
"""
import argparse
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from ase.io import read

sys.path.insert(0, "/anvil/scratch/x-rywang/Au_Cl/scripts")
from analysis_spatial import coordination, CN_CUT  # noqa: E402

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
OUT = f"{ROOT}/analysis/gallery"
R_AU = 1.44
CN_COLOR = {"kink": "#c0392b", "edge": "#e08a2e", "terrace": "#d8b34a", "bulk": "#9aa3ad"}
# The figures are labelled in English, so the CJK fallback face is gone and the default sans is used. It has
# U+00C5 and U+2212 directly; the mathtext Angstrom is kept only so existing call sites need no edit.
CJK = matplotlib.font_manager.FontProperties(family="DejaVu Sans")
AA = r"$\mathrm{\AA}$"


FAMILY_ZH = {k: k for k in ("flat Au(111)", "point defect", "reconstruction-related", "strip step",
                            "vicinal step face", "kink / edge rearrangement", "single-layer island",
                            "single-layer pit", "composite")}
CN_ZH = {"kink": "kink", "edge": "edge", "terrace": "terrace", "bulk": "sub-surface"}


def cn_class(c):
    return "kink" if c <= 6 else "edge" if c <= 8 else "terrace" if c == 9 else "bulk"


SURFACE_BAND = 3.0      # A below the highest atom; the SAME surface set analysis_spatial.column_labels uses


def exposed_atoms(atoms):
    """The surface set the picture emphasises: exactly the set the region analysis assigns columns to.

    It used to be "within 3 A of the highest atom", which on a flat (111) face also caught the second layer 2.4 A
    below. The region analysis was corrected to the un-buried rule; this follows it, so the atoms drawn in full
    colour are the atoms that actually carry a coordination label in the anion maps. Both terraces of a step and
    the over-coordinated step-foot row stay in the set, as they should."""
    return top_atoms(atoms)


REGISTRY_FAMILIES = {"flat Au(111)", "point defect", "reconstruction-related", "strip step",
                     "kink / edge rearrangement", "single-layer island", "single-layer pit", "composite"}
# NOT the vicinal faces: the test below uses the global z axis and the (111) interlayer spacing, which only means
# something when the macroscopic normal IS (111). On Au(211)/(221)/(332)/(554) it fired on a couple of atoms and
# those marks could not be read as hcp stacking. Restoring them there needs a local terrace frame first.


def registry_offset(atoms):
    """Lateral distance from each atom to the nearest atom TWO layers below, minimum image (inf if none).

    0 means the atom sits directly over the one two layers down = hcp termination; the fcc value is a0/sqrt(3)
    ~ 1.70 A. Reporting the distance rather than a yes/no lets a domain wall show its transition atoms instead of
    being forced into one of two bins."""
    pos = atoms.get_positions(); cell = atoms.get_cell().array
    best = np.full(len(pos), np.inf)
    for si in (-1, 0, 1):
        for sj in (-1, 0, 1):
            sh = (si * cell[0] + sj * cell[1])[:2]
            d = np.linalg.norm(pos[:, None, :2] - (pos[None, :, :2] + sh), axis=-1)
            dz = pos[:, None, 2] - pos[None, :, 2]
            m = np.abs(dz - 2 * D111) < 0.7
            dd = np.where(m, d, np.inf)
            best = np.minimum(best, dd.min(axis=1))
    return best


def registry_class(atoms, family=None):
    """0 = fcc-like, 1 = transition / domain wall, 2 = hcp-like, -1 = not applicable."""
    if family not in REGISTRY_FAMILIES: return np.full(len(atoms), -1)
    off = registry_offset(atoms)
    ref = layer_spacing(atoms) / np.sqrt(3.0)       # the fcc lateral offset
    cls = np.full(len(atoms), 0)
    cls[off > 1e6] = 0
    cls[(off < 0.35 * ref)] = 2
    cls[(off >= 0.35 * ref) & (off <= 0.75 * ref)] = 1
    return cls


def registry(atoms, family=None):
    return registry_class(atoms, family) == 2


BAND_FAMILIES = {"single-layer pit", "single-layer island", "point defect", "composite"}
BAND_HALF_FRAC = 0.62   # x the same-layer spacing. One constant, so the regression tests the width the
                        # renderer actually uses rather than a number copied by hand into the test.
ADD_FAMILIES = {"point defect", "single-layer island", "composite"}     # families whose defining atoms sit ON the terrace
MISS_FAMILIES = {"point defect"}                                       # families defined by removing terrace atoms


def added_and_missing(atoms, family=None):
    """(atoms standing above the terrace, empty lattice sites inside the terrace), gated by the structure family.

    A purely geometric detector is not reliable enough to stand alone: on a vicinal face the step-edge row sits
    above the median and reads as an "adatom", on a strip step the whole upper terrace does, and just outside an
    island edge there are sites that look like interior holes. The family comes from the frozen plan, i.e. from how
    the structure was actually built, so it is used to decide WHICH marking is meaningful for this structure, and
    the geometry only decides WHERE."""
    pos = atoms.get_positions(); cell = atoms.get_cell().array
    top = top_atoms(atoms)
    zs = pos[top, 2]
    modal = float(np.median(zs))
    added = np.flatnonzero(top & (pos[:, 2] > modal + 1.0))
    if family not in ADD_FAMILIES or len(added) > 0.30 * top.sum(): added = np.array([], int)
    if family not in MISS_FAMILIES: return added, np.zeros((0, 2))
    terr = pos[top & (np.abs(pos[:, 2] - modal) < 0.8)]
    if len(terr) < 3: return added, np.zeros((0, 2))
    d = np.linalg.norm(terr[:, None, :2] - terr[None, :, :2], axis=-1)
    nn = d[(d > 0.1) & (d < 4.0)]
    if not len(nn): return added, np.zeros((0, 2))
    a0 = float(np.median(nn[nn < np.percentile(nn, 30)]))
    vecs = [np.array([np.cos(t), np.sin(t)]) * a0 for t in np.arange(6) * np.pi / 3]
    # candidate empty sites, then keep only those ringed by terrace atoms (interior holes, not island exteriors)
    cand = np.vstack([terr[:, :2] + v for v in vecs])
    keep = []
    for c in cand:
        dd = np.min([np.linalg.norm(terr[:, :2] + (si * cell[0] + sj * cell[1])[:2] - c, axis=-1).min()
                     for si in (-1, 0, 1) for sj in (-1, 0, 1)])
        if dd < 0.6 * a0: continue                       # occupied
        n_ring = sum(int((np.linalg.norm(terr[:, :2] + (si * cell[0] + sj * cell[1])[:2] - c, axis=-1)
                          < 1.25 * a0).sum()) for si in (-1, 0, 1) for sj in (-1, 0, 1))
        if n_ring >= 4: keep.append(c)
    if keep:
        keep = np.array(keep); uniq = [keep[0]]
        for c in keep[1:]:
            if min(np.linalg.norm(np.array(uniq) - c, axis=1)) > 0.5 * a0: uniq.append(c)
        keep = np.array(uniq)
    else:
        keep = np.zeros((0, 2))
    return added, keep


def draw(ax, pts, colors, depth, emph, r):
    """Colour carries the coordination number. BRIGHTNESS and outline weight carry whether the atom is EXPOSED, so a
    grey atom at the surface (a step foot, over-coordinated because the terrace above leans on it) is never confused
    with a grey atom buried in the slab. Exposed atoms keep full colour and a dark outline; covered ones fade."""
    order = np.argsort(depth)
    dmax = depth.max() if len(depth) else 1.0
    for i in order:
        if emph[i]:
            f, ec, lw, z = 1.0, "#23272c", 0.85, 3000
        else:
            f = float(np.clip(1.0 - (dmax - depth[i]) / 6.0, 0.22, 0.85)) ** 1.15
            ec, lw, z = tuple(np.full(3, 0.96 - 0.28 * f)), 0.45, int(1500 * f)
        c = np.array(matplotlib.colors.to_rgb(colors[i]))
        ax.add_patch(Circle(pts[i], r, facecolor=tuple(c * f + (1 - f) * 0.975), edgecolor=ec,
                            linewidth=lw, zorder=z + int(60 * depth[i])))


TOP_DEPTH = 5.2         # A below the highest atom kept in the top view: two (111) terrace levels, no deep bulk
TILE_TARGET = 26.0      # A, the in-plane extent each tiled view aims for
D111 = 2.4              # A, (111) interlayer spacing -- the quantum the surface height comes in
R_FOOT = 1.75           # A, lateral radius an atom covers when the surface height map is built


R_COVER = 2.35   # A, lateral reach of a neighbour one layer up. In fcc(111) stacking an atom of the layer above sits
N_COVER = 3      # 1.70 A away laterally and a BURIED atom has three of them, while a step-foot atom has only one or
                 # two, so "three higher neighbours within R_COVER" separates buried from merely next to a step.
                 # 2.35 rather than ~1.8 because R2's compressed stripe layer pushes its third neighbour out to
                 # 2.35 A; measured, the choice is safe: under a normal layer the third neighbour is at 1.70 and an
                 # exposed step-foot atom's is at 3.39, and the buried count is identical for 2.05 through 2.50.


def top_atoms(atoms):
    """Atoms not buried under a complete layer: the set whose heights define the surface relief."""
    pos = atoms.get_positions(); cell = atoms.get_cell().array
    n_above = np.zeros(len(pos), int)
    for si in (-1, 0, 1):
        for sj in (-1, 0, 1):
            sh = (si * cell[0] + sj * cell[1])[:2]
            d = np.linalg.norm(pos[:, None, :2] - (pos[None, :, :2] + sh), axis=-1)
            n_above += ((d < R_COVER) & ((pos[None, :, 2] - pos[:, None, 2]) > 0.5)).sum(1)
    return n_above < N_COVER


def height_map(atoms, ng=150):
    """Surface height h(x,y) = the height of the nearest TOP atom (Voronoi over the un-buried atoms only).

    Three definitions were tried. "Highest atom within a fixed disc" let a step-edge atom's disc spill across a narrow
    vicinal terrace and swallow it. "Nearest atom in the top 3 A" put half of a FLAT terrace on the second layer, so
    every flat surface came out looking like a step. A sphere envelope fixed both but its within-atom bumps broke up
    a 4-row vicinal terrace into blobs. Voronoi over un-buried atoms has none of these: every exposed atom owns
    exactly its own patch, terraces keep their true width, and a flat surface is exactly flat."""
    cell = atoms.get_cell().array; pos = atoms.get_positions()[top_atoms(atoms)]
    n1 = max(40, min(ng, int(ng * np.linalg.norm(cell[0][:2]) / TILE_TARGET)))
    n2 = max(40, min(ng, int(ng * np.linalg.norm(cell[1][:2]) / TILE_TARGET)))
    f1, f2 = np.meshgrid((np.arange(n1) + .5) / n1, (np.arange(n2) + .5) / n2, indexing="ij")
    gx = f1 * cell[0][0] + f2 * cell[1][0]
    gy = f1 * cell[0][1] + f2 * cell[1][1]
    best = np.full(gx.shape, np.inf); H = np.zeros(gx.shape)
    for si in (-1, 0, 1):
        for sj in (-1, 0, 1):
            sh = (si * cell[0] + sj * cell[1])[:2]
            for p in pos:
                d2 = (gx - (p[0] + sh[0])) ** 2 + (gy - (p[1] + sh[1])) ** 2
                m = d2 < best
                best[m] = d2[m]; H[m] = p[2]
    return H, gx, gy, (n1, n2)


def flatten_cell(atoms):
    """Rotate in-plane so the longest cell vector lies along +x. Render-only: a vicinal cell is strongly sheared and
    tilted, and drawing it as stored wastes most of the panel on empty corners."""
    a = atoms.copy(); cell = a.get_cell().array.copy()
    iw = int(np.argmax([np.linalg.norm(cell[0][:2]), np.linalg.norm(cell[1][:2])]))
    w = cell[iw][:2]; th = -np.arctan2(w[1], w[0])
    R = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])
    a.set_cell(cell @ R.T, scale_atoms=False); a.set_positions(a.get_positions() @ R.T)
    if iw == 1:                                   # keep a1 as the long, now-horizontal vector
        c = a.get_cell().array.copy(); c[[0, 1]] = c[[1, 0]]
        c[1] *= -1                                # swapping two vectors flips the handedness; mirror a2 to restore it
        a.set_cell(c, scale_atoms=False)
    assert np.linalg.det(a.get_cell().array) > 0, "render cell must stay right-handed"
    a.wrap(pbc=(True, True, False), eps=1e-8)
    return a


MIN_LEVEL_AREA = 0.02   # levels thinner than this are grid noise on a boundary, not a feature of the surface


def levels_of(H):
    """Height quantised into (111) layers relative to the most common terrace level.

    Levels holding less than MIN_LEVEL_AREA of the cell are snapped to the nearest real level: a couple of stray
    pixels on a boundary would otherwise be reported as an extra terrace and turn a flat surface into a 'staircase'."""
    q = np.round((H - np.median(H)) / D111).astype(int)
    vals, cnt = np.unique(q, return_counts=True)
    q = q - vals[np.argmax(cnt)]
    vals, cnt = np.unique(q, return_counts=True)
    keep = vals[cnt / q.size >= MIN_LEVEL_AREA]
    if len(keep) and len(keep) < len(vals):
        q = keep[np.argmin(np.abs(q[..., None] - keep[None, None, :]), axis=-1)]
    return q


GAP_LAYER = 1.2     # A; a gap this wide in the height distribution means genuinely separated atomic levels


def height_gap(H):
    """Largest gap in the sorted height distribution, and the height that splits it.

    This is what separates a STEP from a RAMP. A vicinal face is not a staircase of 2.4 A terraces in this cell: its
    top-atom heights form an even ramp (Au554: ten rows, 0.26 A apart, 2.3 A end to end), because the terrace is a
    (111) micro-facet inclined to the cell's xy plane. Quantising such a surface into 2.4 A layers collapses it to
    'flat', which is why the layer picture had to go."""
    v = np.sort(H.ravel())
    d = np.diff(v)
    if not len(d): return 0.0, float(v[0]) if len(v) else 0.0
    i = int(np.argmax(d))
    return float(d[i]), float(0.5 * (v[i] + v[i + 1]))


def percolates(M):
    """True if the region spans the cell periodically rather than closing on itself.

    Area fraction cannot tell a pit from a step: Pit-19-8x8 occupies 36% of its cell and the old rule called it a
    'step'. What separates them is connectivity -- a strip step runs right across the cell, a pit or an island
    closes. A region that covers every index along one axis after projection spans the cell along the other."""
    return bool(M.any(axis=0).all() or M.any(axis=1).all())


def describe(H, family=None):
    """One Chinese line saying what the outline shows.

    The height field draws the outline; it does NOT name the structure. The family comes from the frozen plan and
    the open/closed distinction from connectivity, because an area fraction alone mislabels a broad shallow pit as
    a step."""
    rng = float(H.max() - H.min())
    if rng < 0.35:
        return "uniform height, no relief; it differs from the flat slab in stacking or in-plane registry, not in height"
    gap, cut = height_gap(H)
    if gap < GAP_LAYER:
        return (f"the local (111) terrace is inclined to the macroscopic face, relief {rng:.1f} A in the cell; "
                "the step sits on the periodic seam, so the height distribution shows no break")
    v = np.sort(H.ravel())
    cuts = [0.5 * (a + b) for a, b in zip(v[:-1], v[1:]) if b - a >= GAP_LAYER]
    base = float(np.median(H))
    hi_m = H > base + 0.5 * GAP_LAYER
    lo_m = H < base - 0.5 * GAP_LAYER
    part = []
    for m, raised in ((hi_m, True), (lo_m, False)):
        if not m.any(): continue
        if percolates(m):
            part.append(("upper terrace spanning the cell" if raised else "lower terrace spanning the cell")
                        + " (one side of a strip step)")
        else:
            part.append("closed protrusion (island or adatom cluster)" if raised
                        else "closed depression (pit or vacancy)")
    if len(cuts) >= 2:
        part.append(f"{len(cuts)+1} separated height levels")
    if not part:
        return f"relief {rng:.1f} A, no clear high or low region resolved"
    return "  ·  ".join(part) + f"  ·  level gap {gap:.1f} A"


PLAN_FILL = {2: "#dfb264", 1: "#e9c98f", 0: "#f2efe9", -1: "#bcd3dd", -2: "#9cbecd", -3: "#86adbf"}


LAYER_DZ = 0.6          # A; two atoms are in the same atomic layer if their heights differ by less than this


def layer_spacing(atoms, sel=None):
    """Median in-plane nearest-neighbour distance between atoms of the SAME atomic layer.

    This is the only lattice scale the drawing may use for connectivity: two surface atoms are bonded neighbours
    when they are this far apart, so discs of radius 0.62 x this overlap exactly for bonded pairs.

    It replaces an nn_spacing() that took the nearest neighbour in PROJECTION over all un-buried atoms and so
    mixed layers. An adatom sits in a hollow a0/sqrt(3) = 1.6975 A laterally from each of the three terrace atoms
    beneath it, and on A3, C1 and the compact pits more than half the un-buried atoms had such a cross-layer
    partner as their projected nearest, so the median came back 1.6975 A instead of 2.9401 A. Discs of radius
    0.62 x 1.6975 = 1.05 A cannot bridge a 2.94 A bond: A3's three-atom cluster was drawn as three separate
    circles (atom graph 1 component, disc graph 3) and Pit-7-compact's floor as twelve. Restricting the pair to
    |dz| < LAYER_DZ makes it the same-layer spacing it was always meant to be; R2's compressed stripe layer still
    reports its own, smaller value, and that value no longer leaks into any other structure's scale."""
    pos = atoms.get_positions(); cell = atoms.get_cell().array
    p = pos[top_atoms(atoms) if sel is None else sel]
    if len(p) < 2: return 2.94
    best = np.full(len(p), np.inf)
    for si in (-1, 0, 1):
        for sj in (-1, 0, 1):
            sh = (si * cell[0] + sj * cell[1])[:2]
            d = np.linalg.norm(p[:, None, :2] - (p[None, :, :2] + sh), axis=-1)
            d = np.where((np.abs(p[:, None, 2] - p[None, :, 2]) < LAYER_DZ) & (d > 0.1), d, np.inf)
            best = np.minimum(best, d.min(axis=1))
    best = best[np.isfinite(best)]
    return float(np.median(best)) if len(best) else 2.94


def blob_distance(pts_xy, gx, gy, cell):
    """Distance from every grid point to the nearest of pts_xy, minimum image.

    Call this ONCE, on the base-cell grid. The result is periodic by construction, so every tiled copy of the
    schematic must be drawn by translating this array -- never by re-evaluating the function at a shifted query
    position. The image search runs over -1,0,+1 only, so a query two cells out lost the atoms entirely: the
    third copy of Step-8x2 came out 0.40% covered against 53.18% for the first, i.e. its upper terrace vanished
    from a picture of a perfectly periodic structure."""
    best = np.full(gx.shape, np.inf)
    for si in (-1, 0, 1):
        for sj in (-1, 0, 1):
            sh = (si * cell[0] + sj * cell[1])[:2]
            for p in pts_xy:
                np.minimum(best, np.hypot(gx - (p[0] + sh[0]), gy - (p[1] + sh[1])), out=best)
    return best


def feature_sets(atoms):
    """(raised atoms, lowered atoms, vacant sites, modal terrace height) -- the ATOMS that define the feature.

    The outline is built from these, not from a height partition. The nearest-atom height partition has its own
    definition but it is NOT the outline of the atomic structure: a lower atom can win territory between two
    upper atoms and cut a connected cluster in two. Verified on the real coordinates -- A3's three adatoms and
    Island-7-elongated's seven island atoms are each ONE connected component at 2.94 A, yet the partition drew
    them as separate blobs. Contouring a distance field around the atoms themselves cannot do that."""
    pos = atoms.get_positions()
    top = top_atoms(atoms)
    modal = float(np.median(pos[top, 2]))
    hi = pos[top & (pos[:, 2] > modal + 1.0)][:, :2]
    lo = pos[top & (pos[:, 2] < modal - 1.0)][:, :2]
    return hi, lo, modal


def wrap_offset(L):
    """Always zero. Kept as a named function so the reason is recorded rather than rediscovered.

    An earlier version measured a net level change across the periodic seam and tiled the schematic with it, to draw
    a vicinal face as a descending staircase. That is wrong here: every cell in this set has a1_z = a2_z = 0, so the
    surface must return to the same height after one cell by construction, and there is no net descent to add back.
    The measured 'offset' was just the first and last tenth of the cell landing on different terraces, which also
    mislabelled every strip-step cell as a vicinal."""
    return 0, 0


def wrap_pad(gx, gy, A, cell):
    """Repeat the first row and column one period out, so tiled fills abut instead of leaving a hairline.

    The grid samples cell centres, from 0.5/n to (n-0.5)/n of the cell, so neighbouring tiles of a periodic
    field stopped a half cell short of each other and a white seam ran between every copy."""
    gx = np.concatenate([gx, gx[:1] + cell[0][0]], axis=0)
    gy = np.concatenate([gy, gy[:1] + cell[0][1]], axis=0)
    A = np.concatenate([A, A[:1]], axis=0)
    gx = np.concatenate([gx, gx[:, :1] + cell[1][0]], axis=1)
    gy = np.concatenate([gy, gy[:, :1] + cell[1][1]], axis=1)
    A = np.concatenate([A, A[:, :1]], axis=1)
    return gx, gy, A


def tiles(gx, gy, A, cell, i_range, j_range):
    """Yield (x, y, A) for each periodic copy of a base-cell field.

    The single place a tiled copy is produced, shared by the renderer and by check_gallery_drawing.py, so
    "every copy is a pure translation" is enforced rather than asserted in a comment. The ARRAY yielded is the
    same object every time; only the coordinates move. Re-evaluating the field at a shifted query instead is
    what blanked the third copy of Step-8x2 (53.18% / 51.11% / 0.40% coverage), because the image search behind
    it only reached +-1 cell."""
    for i in i_range:
        for j in j_range:
            sh = (i * cell[0] + j * cell[1])[:2]
            yield gx + sh[0], gy + sh[1], A


def smooth_periodic(A, sigma_cells=1.3):
    """Periodic Gaussian blur, in GRID cells. Applied to the distance field before it is both contoured and
    sliced, so rounding the per-atom scallops off the outline cannot make the outline and the section disagree."""
    n1, n2 = A.shape
    k1 = np.fft.fftfreq(n1)[:, None]; k2 = np.fft.fftfreq(n2)[None, :]
    return np.fft.ifft2(np.fft.fft2(A) * np.exp(-2 * (np.pi ** 2) * (sigma_cells ** 2)
                                                * (k1 ** 2 + k2 ** 2))).real


def schematic_surface(at, gx, gy):
    """The cartoon surface the plan view is drawn from: a flat terrace with each feature's disc union raised or
    sunk to that feature's MEASURED height.

    The point of returning a surface rather than only an outline is that the plan view and the section below it
    are then two views of ONE object. Before this, the outline was a union of discs while the section was sliced
    out of the nearest-atom height partition; the two have different definitions, so a cluster could be drawn
    connected above and show a notch between its atoms below."""
    cell = at.get_cell().array
    hi, lo, modal = feature_sets(at)
    pos = at.get_positions(); top = top_atoms(at)
    R = 0.62 * layer_spacing(at)
    Z = np.zeros(gx.shape); layers = []
    for pts, sgn in ((hi, +1.0), (lo, -1.0)):
        if not len(pts):
            layers.append(None); continue
        zz = pos[top & (sgn * (pos[:, 2] - modal) > 1.0), 2]
        h = float(np.median(zz) - modal) if len(zz) else sgn * D111
        D = smooth_periodic(blob_distance(pts, gx, gy, cell))   # ONCE; tiled copies translate this array
        Z[D < R] = h
        layers.append(dict(D=D, h=h, n=len(pts)))
    return Z, layers, R, modal


def _runs(mask):
    """(total samples set, longest CYCLIC run) -- a profile is a closed loop, so the run may wrap the seam."""
    n = len(mask); tot = int(mask.sum())
    if tot == 0 or tot == n: return tot, tot
    best = cur = 0
    for k in range(2 * n):
        if mask[k % n]:
            cur += 1; best = max(best, cur)
        else:
            cur = 0
    return tot, min(best, n)


def _line(F, ax_, idx):
    return F[:, idx] if ax_ == 0 else F[idx, :]


def _score(pr, hi_lvl, lo_lvl, target):
    """Rank a candidate line by the LONGEST unbroken crossing, then by total samples inside the feature.

    Ranking on the total alone rewards a line that shaves the scalloped rim of a blob: on C2 the best 'crosses
    both' line by that measure entered the pit three times for 7 samples each and the island twice for 4, which
    draws as several small dents rather than one pit and one island."""
    th, rh = _runs(pr > hi_lvl); tl, rl = _runs(pr < lo_lvl)
    if target == "hi": return (rh, th)
    if target == "lo": return (rl, tl)
    return (min(rh, rl), rh + rl)


CONTIG = 0.8            # a crossing counts as ONE passage through the feature if this much of it is unbroken
JOINT = 0.6             # and a joint line must give each feature this much of its own best line's crossing


def pick_cuts(F, hi_lvl, lo_lvl):
    """One section line, or TWO when no lattice-aligned line properly crosses both a raised and a sunken feature.

    C2 carries an island and a pit on different rows of the same cell. Nothing along a1 or a2 passes cleanly
    through both, and a line that merely clips their rims is worse than no line: it draws one pit as three.
    Two lines, A-A' through the raised feature and B-B' through the sunken one, is the honest answer and costs
    one panel. They are forced onto the same lattice direction so the two sections share an abscissa."""
    want_hi = bool((F > hi_lvl).any()); want_lo = bool((F < lo_lvl).any())
    if not (want_hi or want_lo): return []
    cand = [(a, i) for a in (0, 1) for i in range(F.shape[1 - a])]

    def best(target, axis=None):
        c = [x for x in cand if axis is None or x[0] == axis]
        return max(c, key=lambda x: _score(_line(F, *x), hi_lvl, lo_lvl, target))

    if want_hi and want_lo:
        bh = best("hi"); bl = best("lo")
        sh = _score(_line(F, *bh), hi_lvl, lo_lvl, "hi")[0]
        sl = _score(_line(F, *bl), hi_lvl, lo_lvl, "lo")[0]
        bb = best("both"); pr = _line(F, *bb)
        th, rh = _runs(pr > hi_lvl); tl, rl = _runs(pr < lo_lvl)
        if rh >= CONTIG * th and rl >= CONTIG * tl and rh >= JOINT * sh and rl >= JOINT * sl:
            return [dict(axis=bb[0], index=bb[1], label="A", target="both")]
        pick = None
        for a in (0, 1):
            ch = best("hi", a); cl = best("lo", a)
            s = (min(_score(_line(F, *ch), hi_lvl, lo_lvl, "hi")[0],
                     _score(_line(F, *cl), hi_lvl, lo_lvl, "lo")[0]), )
            if pick is None or s > pick[0]: pick = (s, ch, cl)
        _, ch, cl = pick
        return [dict(axis=ch[0], index=ch[1], label="A", target="hi"),
                dict(axis=cl[0], index=cl[1], label="B", target="lo")]
    t = "hi" if want_hi else "lo"
    b = best(t)
    if _score(_line(F, *b), hi_lvl, lo_lvl, t)[0] == 0: return []
    return [dict(axis=b[0], index=b[1], label="A", target=t)]


def centre_roll(pr, hi_lvl, lo_lvl, ramp=False):
    """How far to slide the display origin A along the line so the feature is not cut by the panel edge.

    This is a shift of the WHOLE display -- A, A', the profile and the side-view atom band all read the same
    offset out of cut_frame -- not a roll applied to one curve. Rolling only the profile is what put Pit-19's
    terrace in the middle of a picture whose plan view said the line starts on the terrace."""
    n = len(pr)
    if ramp:
        # put the riser at ~85% along, so one period reads terrace -> step -> the start of the next terrace
        # instead of losing the step under the panel's right edge
        d = np.diff(np.append(pr, pr[0]))
        return int(int(0.85 * n) - int(np.argmax(np.abs(d))))
    m = (pr > hi_lvl) | (pr < lo_lvl)
    idx = np.flatnonzero(m)
    if not len(idx) or len(idx) == n: return 0
    a = 2 * np.pi * idx / n
    c = int(round((np.arctan2(np.sin(a).mean(), np.cos(a).mean()) % (2 * np.pi)) * n / (2 * np.pi))) % n
    return int(n // 2 - c)


RAMP_MIN = 0.35         # A; below this the surface is flat and there is no downhill direction to section

# Miller indices of the vicinal slabs, taken from how they were built, not inferred from the picture.
VICINAL_MILLER = {"Au211": (2, 1, 1), "Au221": (2, 2, 1), "Au332": (3, 3, 2), "Au554": (5, 5, 4)}


def miller_angle_to_111(hkl):
    """Angle between the macroscopic (hkl) normal and (111), in degrees, for a cubic lattice.

    A crystallographic constant, arccos[(h+k+l)/(sqrt(3)sqrt(h^2+k^2+l^2))] -- NOT a measurement off this
    figure. The section used to print an angle computed as arctan(ptp(profile)/(L - grid step)), i.e. the
    apparent slope of the height-partition field along one lattice direction over nearly a whole period. That
    happens to land near the right value when the cut direction is the steepest descent and the terrace is
    wide, which is why Au554 came out 5.6 against 5.77 degrees, but it is not the same quantity and it is not
    generally close: Au211 printed 25.8 against a nominal 19.47."""
    h, k, l = hkl
    return float(np.degrees(np.arccos((h + k + l) / (np.sqrt(3.0) * np.sqrt(h * h + k * k + l * l)))))


def ramp_cut(H):
    """The section line for a surface with no separated levels: the lattice line of greatest relief.

    On a vicinal face this is the across-terrace direction, and the profile it returns is the staircase -- the
    inclined (111) micro-facet followed by the riser at the periodic seam. The plan view's colour gradient shows
    the height distribution but cannot show a step, which is why these four faces needed a real section."""
    best = None
    for ax_ in (0, 1):
        for idx in range(H.shape[1 - ax_]):
            pr = H[:, idx] if ax_ == 0 else H[idx, :]
            if best is None or float(np.ptp(pr)) > best[0]: best = (float(np.ptp(pr)), ax_, idx)
    if best is None or best[0] < RAMP_MIN: return []
    return [dict(axis=best[1], index=best[2], label="A", target="ramp")]


def cut_frame(cell, ci):
    """(rA, t_hat, n_hat, L, P) for a section line: where A is, the unit vector A->A', the in-plane normal, the
    line's length and the true PERPENDICULAR period across it.

    Everything downstream -- which atoms are in the band, how wide the band really is, and where each atom lands
    in the side view -- comes from this one frame. The band used to be selected from a fractional-coordinate
    difference times the other cell vector's length, which in a 60 degree cell overstates the perpendicular
    half-width by 1/sin(60): a band captioned +-3.23 A was in fact +-2.80 A. And the side view then projected on
    the GLOBAL x axis whatever direction the line ran, so for a line along (0.5, 0.866) the picture was not a
    view along A-A' at all."""
    ax_ = ci["axis"]
    along = cell[ax_][:2]; across = cell[1 - ax_][:2]
    L = float(np.linalg.norm(along)); t = along / L
    n = np.array([-t[1], t[0]])
    s0 = -(ci.get("roll", 0) / max(1, ci.get("n", 1))) * L        # the one place the display offset is applied
    rA = ci["frac"] * across + s0 * t
    # wrap A back into the base cell so its label lands inside the drawn panel. A lattice translation changes
    # neither s (taken mod L) nor the perpendicular offset (taken mod P), so nothing downstream shifts.
    C = np.array([cell[0][:2], cell[1][:2]])
    f = rA @ np.linalg.inv(C)
    rA = (f - np.floor(f)) @ C
    return rA, t, n, L, abs(float(across @ n))


def cut_band(pos, cell, ci, half):
    """(mask of atoms within `half` of the line, coordinate along it, signed offset across it) -- one frame.

    ONE periodic image per atom, chosen once, then BOTH coordinates read off it. The two used to be folded
    independently: the offset across the line was wrapped into the nearest band while the coordinate along the
    line kept the atom's original position. In a sheared cell the across vector has a component ALONG the line,
    so that is not allowed. With a = (10, 0) and b = (5, 8.6603), the same atom written as r and as r + b came
    out at s = 9.5 and s = 4.5 -- half a period apart, from nothing but a change of representation. Here the
    image index m is picked from the perpendicular distance and r - m*across is used for both.

    `half` is a TRUE perpendicular distance in angstrom, so the caption on the picture is the band that was
    actually taken. The signed offset doubles as the side view's depth cue."""
    rA, t, n, L, _ = cut_frame(cell, ci)
    across = cell[1 - ci["axis"]][:2]
    Pn = float(across @ n)                       # signed period across the line
    r = pos[:, :2] - rA
    m = np.round((r @ n) / Pn)
    r = r - m[:, None] * across[None, :]
    dv = r @ n
    s = (r @ t) % L
    return np.abs(dv) < half, s, dv


def cut_profile(F, ci, cell):
    """(s, profile) sampled along the line, in the order A -> A'.

    NOT rolled. The profile used to be circularly shifted to centre the feature while the A and A' labels on the
    plan view stayed put, so the two stopped corresponding: Pit-19-8x8 is sampled terrace(16) -> pit floor(75) ->
    terrace(8) along its marked A-A', and the roll turned that into pit(38) -> terrace(24) -> pit(37), which is
    why a pit read as a central bump. Centring is not worth breaking the correspondence the labels promise."""
    pr = np.roll(_line(F, ci["axis"], ci["index"]).astype(float), ci.get("roll", 0))
    L = float(np.linalg.norm(cell[ci["axis"]][:2]))
    return (np.arange(len(pr)) + 0.5) / len(pr) * L, pr, L


def prefer_row(F, ci, hi_lvl, lo_lvl, prefer_xy, cell):
    """Slide a chosen section line onto the row that also passes through `prefer_xy`, if that costs nothing.

    C1's seven island atoms sit at the same height as the step's upper terrace, so the height field holds no
    separate feature for them and the line was free to miss the island the structure is named for. The build
    knows where they are; if a line through them crosses the height feature nearly as well as the best line
    does, it is strictly the more informative cut."""
    if prefer_xy is None or not len(prefer_xy): return
    ax_ = ci["axis"]; n = F.shape[1 - ax_]
    C = np.array([cell[0][:2], cell[1][:2]])
    fr = (np.asarray(prefer_xy) @ np.linalg.inv(C)) % 1.0
    a = 2 * np.pi * fr[:, 1 - ax_]
    f = (np.arctan2(np.sin(a).mean(), np.cos(a).mean()) % (2 * np.pi)) / (2 * np.pi)
    want = int(np.clip(round(f * n - 0.5), 0, n - 1))
    # Passing through the island necessarily breaks the lower terrace's run in two -- that break IS the island --
    # so a longest-run test would reject exactly the line wanted. Accept while the feature is still well shown:
    # at least half its best total, and one unbroken passage worth a quarter of that total.
    base = _score(_line(F, ax_, ci["index"]), hi_lvl, lo_lvl, ci["target"])[1]
    alt_run, alt_tot = _score(_line(F, ax_, want), hi_lvl, lo_lvl, ci["target"])
    if alt_tot >= 0.5 * base and alt_run >= 0.25 * base:
        ci["index"] = want; ci["through_prefer"] = True


def surface_model(at, prefer_xy=None):
    """Everything the three panels share: the height map, the cartoon surface, and the section line(s).

    Computed once, before anything is drawn, because the schematic, the section and the side-view atom band must
    all be the same line in the same frame, and the figure size depends on the band."""
    cell = at.get_cell().array
    H, gx, gy, _ = height_map(at, ng=110)
    H = H - np.median(H)
    hi, lo, modal = feature_sets(at)
    gap, _ = height_gap(H)
    v = np.sort(H.ravel())
    levels = [float(0.5 * (a + b)) for a, b in zip(v[:-1], v[1:]) if b - a >= GAP_LAYER]
    Z = layers = None; R = 0.62 * layer_spacing(at)
    if levels and (len(hi) or len(lo)):
        mode = "blob"
        Z, layers, R, modal = schematic_surface(at, gx, gy)
        F = Z; hi_lvl, lo_lvl = 0.5 * D111, -0.5 * D111
    elif levels:
        mode = "level"; F = H; hi_lvl, lo_lvl = 0.5 * gap, -0.5 * gap
    else:
        mode = "ramp"; F = H; hi_lvl, lo_lvl = np.inf, -np.inf
    cuts = ramp_cut(H) if mode == "ramp" else pick_cuts(F, hi_lvl, lo_lvl)
    for c in cuts:
        if mode != "ramp" and len(cuts) == 1:
            prefer_row(F, c, hi_lvl, lo_lvl, prefer_xy, cell)
        c["frac"] = (c["index"] + 0.5) / F.shape[1 - c["axis"]]
        pr = _line(F, c["axis"], c["index"])
        c["n"] = int(len(pr))
        c["roll"] = centre_roll(pr, hi_lvl, lo_lvl, ramp=(mode == "ramp"))
    return dict(H=H, gx=gx, gy=gy, Z=Z, layers=layers, R=R, mode=mode, cuts=cuts, levels=levels,
                gap=gap, field=F, hi=hi, lo=lo, modal=modal)


def section_panel(axs, sm, ci, cell, reps, tag, hkl=None):
    """One section, drawn from the SAME field the plan view above it is drawn from and in the order A -> A'."""
    F = sm["Z"] if sm["mode"] == "blob" else sm["H"]
    s, pr, L = cut_profile(F, ci, cell)
    s = np.concatenate([s + k * L for k in range(reps)])
    p = np.tile(pr, reps)
    s = np.append(s, reps * L); p = np.append(p, p[0])
    body = 3.4
    floor = float(p.min()) - body
    axs.fill_between(s, p, floor, step="mid", color="#e6e2da", zorder=1)
    axs.step(s, p, where="mid", color="#2b3137", lw=2.1, zorder=3)
    axs.axhline(0.0, color="#b9b3a7", lw=0.9, ls=(0, (4, 3)), zorder=2)
    if ci.get("target") != "ramp":
        axs.text(reps * L, -0.55, "terrace baseline ", fontproperties=CJK, fontsize=7.6, color="#9aa1a8",
                 ha="right", va="top", zorder=5)
    axs.plot([0, reps * L], [floor] * 2, color="#c6c1b7", lw=1.0, zorder=3)
    top = float(p.max())
    # NO drawn "ion-accessible boundary": that surface is computed from SION and varies over the relief, which is
    # the whole point of the spatial analysis, so a fixed horizontal line would contradict it.
    axs.text(reps * L * 0.01, top + 0.5, "\u2191 electrolyte side", fontproperties=CJK, fontsize=8.6,
             color="#2e7d9a", va="bottom")
    a, b = tag
    axs.text(reps * L * 0.008, floor + 0.35, a, fontsize=9.5, color="#23272c", ha="left", va="bottom", zorder=5)
    axs.text(reps * L * 0.992, floor + 0.35, b, fontsize=9.5, color="#23272c", ha="right", va="bottom", zorder=5)
    extra = 6.6
    if ci.get("target") == "ramp":
        # The staircase, marked QUALITATIVELY: which part is the inclined (111) micro-facet and where the riser
        # is. No number is read off the height field here. The terrace width and the step height need a chosen
        # atom row and a stated measurement direction before they mean anything, and the "tilt" printed before
        # was the apparent slope of the height partition along one lattice vector, which is not the miscut
        # angle. The one number shown is the nominal (hkl)-to-(111) angle from the Miller indices.
        d = np.diff(np.append(pr, pr[0]))
        k = int(np.argmax(np.abs(d)))
        s_seam = float((k + 1.0) / len(pr) * L)
        for rep in range(reps):
            x = s_seam + rep * L
            axs.annotate("", xy=(x, float(pr.max())), xytext=(x, float(pr.min())), zorder=6,
                         arrowprops=dict(arrowstyle="<|-|>", lw=1.3, color="#b0543a", shrinkA=0, shrinkB=0))
        axs.text(s_seam, top + 1.7, "step", fontproperties=CJK, fontsize=8.4,
                 color="#b0543a", ha="center", va="bottom")
        m0, m1 = 0.06 * L, 0.78 * L
        axs.annotate("", xy=(m1, top + 3.2), xytext=(m0, top + 3.2), zorder=6,
                     arrowprops=dict(arrowstyle="<|-|>", lw=1.2, color="#1d4e8f", shrinkA=0, shrinkB=0))
        axs.text(0.5 * (m0 + m1), top + 3.5, "local (111) terrace", fontproperties=CJK, fontsize=8.2,
                 color="#1d4e8f", ha="center", va="bottom")
        extra = 9.4
        if hkl is not None:
            # left-aligned, in AXES fractions. Centred under the arrow it ran off the panel's left edge, and
            # placed in data units the two lines sat a few pixels apart on this short, non-aspect-locked panel.
            h, k_, l_ = hkl
            axs.text(0.0, 1.0, f"nominal ({h}{k_}{l_})\u2013(111) angle {miller_angle_to_111(hkl):.2f}\u00b0",
                     transform=axs.transAxes, fontproperties=CJK, fontsize=7.2, color="#1d4e8f",
                     ha="left", va="top")
            axs.text(0.0, 0.90, "(from the plane normals, not measured here)", transform=axs.transAxes,
                     fontproperties=CJK, fontsize=6.8, color="#9aa1a8", ha="left", va="top")
    axs.set_xlim(0, reps * L); axs.set_ylim(floor - 1.0, top + extra)
    return L


def simple_schematic(fig, cell_spec, at, sm, hkl=None):
    """The plain outline: no atoms at all. Top = plan view of the cartoon surface with its boundaries stroked,
    bottom = one section per marked line, cut through that SAME surface."""
    import matplotlib.gridspec as mgs
    cell = at.get_cell().array
    H, gx, gy, cuts = sm["H"], sm["gx"], sm["gy"], sm["cuts"]
    nsec = max(1, len(cuts))
    inner = mgs.GridSpecFromSubplotSpec(1 + nsec, 1, subplot_spec=cell_spec,
                                        height_ratios=[2.35] + [1.0] * nsec, hspace=0.42)

    # Smooth before contouring so a height-contoured outline reads as one clean curve instead of tracing the
    # atomic scallops. sigma is fixed in GRID cells, not angstrom: a vicinal face has terrace stripes barely 1 A
    # wide and a fixed 0.8 A kernel flattened them below the contour level, leaving a blank panel.
    n1g, n2g = H.shape
    k1 = np.fft.fftfreq(n1g)[:, None]; k2 = np.fft.fftfreq(n2g)[None, :]
    Hs = np.fft.ifft2(np.fft.fft2(H) * np.exp(-2 * (np.pi ** 2) * (1.5 ** 2) * (k1 ** 2 + k2 ** 2))).real

    t1 = max(1, min(3, int(round(TILE_TARGET / np.linalg.norm(cell[0][:2])))))
    t2 = max(1, min(3, int(round(TILE_TARGET / np.linalg.norm(cell[1][:2])))))
    # One tile of margin beyond the displayed window. A section whose A end sits mid-cell has its
    # A' end one period further on, which fell outside a window drawn to exactly t1 x t2 tiles, so
    # the far label vanished from the picture.
    d1, d2 = range(-1, t1 + 1), range(-1, t2 + 1)
    axp = fig.add_subplot(inner[0]); axp.set_aspect("equal"); axp.axis("off")
    cmap = plt.get_cmap("RdYlBu_r")
    if sm["mode"] == "blob":
        # OUTLINE FROM THE ATOMS: the level set of the distance to the nearest feature atom. The nearest-atom
        # height partition, which this replaced, let a lower atom win territory between two upper ones and so
        # split connected clusters. The distance field is computed ONCE on the base cell and TRANSLATED for each
        # tiled copy -- re-evaluating it at a shifted query lost the atoms beyond the +-1 image search and blanked
        # the third copy.
        axp.set_facecolor("none")
        for i in d1:
            for j in d2:
                sh = (i * cell[0] + j * cell[1])[:2]
                axp.fill([sh[0], sh[0] + cell[0][0], sh[0] + cell[0][0] + cell[1][0], sh[0] + cell[1][0]],
                         [sh[1], sh[1] + cell[0][1], sh[1] + cell[0][1] + cell[1][1], sh[1] + cell[1][1]],
                         color=PLAN_FILL[0], zorder=0)
        for lay, col in zip(sm["layers"], (PLAN_FILL[1], PLAN_FILL[-1])):
            if lay is None: continue
            px, py, pD = wrap_pad(gx, gy, lay["D"], cell)
            for x, y, D in tiles(px, py, pD, cell, d1, d2):
                axp.contourf(x, y, D, levels=[0, sm["R"]], colors=[col], zorder=1)
                axp.contour(x, y, D, levels=[sm["R"]], colors="#23272c",
                            linewidths=2.2, linestyles="solid", zorder=3)
    elif sm["mode"] == "level":
        cuts_h = sm["levels"]
        edges = [Hs.min() - 1] + cuts_h + [Hs.max() + 1]
        mids = [0.5 * (edges[k] + edges[k + 1]) for k in range(len(edges) - 1)]
        span = max(max(abs(m) for m in mids), 1.2)
        cols = [cmap(0.5 + 0.40 * m / span) for m in mids]
        px, py, pH = wrap_pad(gx, gy, Hs, cell)
        for x, y, Hh in tiles(px, py, pH, cell, d1, d2):
            axp.contourf(x, y, Hh, levels=edges, colors=cols, zorder=1)
            axp.contour(x, y, Hh, levels=cuts_h, colors="#23272c",
                        linewidths=2.2, linestyles="solid", zorder=3)
    else:
        rng = max(np.ptp(Hs), 1e-6)
        px, py, pH = wrap_pad(gx, gy, Hs, cell)
        for x, y, Hh in tiles(px, py, pH, cell, d1, d2):
            axp.pcolormesh(x, y, Hh, cmap=cmap, vmin=-0.6 * rng, vmax=0.6 * rng,
                           shading="gouraud", zorder=1)
        if np.ptp(Hs) > RAMP_MIN:                                  # mark which way the inclined terrace runs downhill
            g1 = float(np.mean(np.gradient(Hs, axis=0))); g2 = float(np.mean(np.gradient(Hs, axis=1)))
            d = -(g1 * cell[0][:2] / np.linalg.norm(cell[0][:2]) + g2 * cell[1][:2] / np.linalg.norm(cell[1][:2]))
            if np.linalg.norm(d) > 0:
                d = d / np.linalg.norm(d) * 0.22 * float(np.ptp(gx))
                c0 = np.array([gx.mean() * t1, gy.mean() * t2])
                axp.annotate("", xy=c0 + d, xytext=c0 - d, zorder=4,
                             arrowprops=dict(arrowstyle="-|>", lw=2.0, color="#23272c"))
                axp.text(*(c0 + 1.25 * d), "downhill", fontproperties=CJK, fontsize=9.0, color="#23272c",
                         ha="center", va="center", zorder=5)
    o = np.zeros(2)
    axp.plot(*zip(o, cell[0][:2], cell[0][:2] + cell[1][:2], cell[1][:2], o),
             color="#9aa1a8", lw=0.9, ls=(0, (4, 3)), zorder=4)
    px, py, _ = wrap_pad(gx, gy, H, cell)
    X = np.concatenate([(px + (i * cell[0] + j * cell[1])[0]).ravel() for i in range(t1) for j in range(t2)])
    Y = np.concatenate([(py + (i * cell[0] + j * cell[1])[1]).ravel() for i in range(t1) for j in range(t2)])
    xlo, xhi, ylo, yhi = X.min(), X.max(), Y.min(), Y.max()
    for ci in cuts:
        rA, _, _, _, _ = cut_frame(cell, ci)
        for q in (rA, rA + cell[ci["axis"]][:2]):
            xlo = min(xlo, q[0] - 1.2); xhi = max(xhi, q[0] + 1.2)
            ylo = min(ylo, q[1] - 1.2); yhi = max(yhi, q[1] + 1.2)
    axp.set_xlim(xlo, xhi); axp.set_ylim(ylo, yhi)
    # the window is extended to hold the A' marker, so a fixed "t1 x t2 cells" would no longer be the truth
    axp.set_title("outline \u00b7 plan view (dashes = one cell)", fontproperties=CJK,
                  fontsize=9.5, color="#4a5158", pad=3)

    # the section line(s) on the plan, so the reader can see WHERE each was taken
    for ci in cuts:
        # A and A' come out of cut_frame, the same call the profile and the atom band use, so the label can only
        # ever mark where the section actually starts.
        rA, tv, _, L, _ = cut_frame(cell, ci)
        along = cell[ci["axis"]][:2]
        a, b = (ci["label"], ci["label"] + "′")
        for i in d1:
            for j in d2:
                sh = (i * cell[0] + j * cell[1])[:2]
                p0 = rA + sh; p1 = p0 + along
                axp.plot([p0[0] - along[0], p1[0]], [p0[1] - along[1], p1[1]],
                         color="#23272c", lw=1.3, ls=(0, (6, 3)), zorder=6)
        axp.text(rA[0], rA[1], " " + a, fontsize=9, color="#23272c", va="center", ha="left", zorder=7)
        p1 = rA + along
        axp.text(p1[0], p1[1], b + " ", fontsize=9, color="#23272c", va="center", ha="right", zorder=7)

    if not cuts:
        axs = fig.add_subplot(inner[1]); axs.axis("off")
        use_a1 = float(np.ptp(H.max(axis=1))) >= float(np.ptp(H.max(axis=0)))
        prof = H.max(axis=1) if use_a1 else H.max(axis=0)
        L = float(np.linalg.norm(cell[0][:2] if use_a1 else cell[1][:2]))
        reps = max(1, min(2, int(round(TILE_TARGET / L))))
        s = np.append(np.concatenate([(np.arange(len(prof)) + .5) / len(prof) * L + k * L for k in range(reps)]),
                      reps * L)
        p = np.append(np.tile(prof, reps), prof[0])
        floor = float(p.min()) - 4.8
        axs.fill_between(s, p, floor, step="mid", color="#e6e2da", zorder=1)
        axs.step(s, p, where="mid", color="#2b3137", lw=2.1, zorder=3)
        axs.axhline(0.0, color="#b9b3a7", lw=0.9, ls=(0, (4, 3)), zorder=2)
        axs.set_xlim(0, reps * L); axs.set_ylim(floor - 0.4, float(p.max()) + 6.6)
        axs.set_title(f"projected envelope, not a section \u00b7 {reps}\u00d7", fontproperties=CJK,
                      fontsize=9.0, color="#4a5158", pad=2)
    else:
        for k, ci in enumerate(cuts):
            axs = fig.add_subplot(inner[1 + k]); axs.axis("off")
            L = float(np.linalg.norm(cell[ci["axis"]][:2]))
            reps = max(1, min(2, int(round(TILE_TARGET / L))))
            a, b = ci["label"], ci["label"] + "′"
            section_panel(axs, sm, ci, cell, reps, (a, b), hkl=hkl)
            src = "outline surface" if sm["mode"] == "blob" else "height field"
            cap = f"section {a}\u2013{b} \u00b7 cut from the {src} above \u00b7 {reps}\u00d7"
            if ci.get("target") == "ramp":
                cap += "\nheight in cell coordinates; no widths quoted"
            if k == len(cuts) - 1:
                cap += "\nsolid below = substrate only"
                if len(cuts) == 2: cap += "; island and pit are not on one lattice line"
            axs.set_title(cap, fontproperties=CJK, fontsize=8.8, color="#4a5158", pad=2)
    return describe(H)
LIB = f"{ROOT}/03_pilot/all_defect_structures"
# child -> the parent POSCAR it was built from. Taken from the build scripts, not guessed: the diff below IS the
# build-time atom mapping. Nothing geometric can recover it -- C1's seven island atoms sit at exactly the same
# height as the step's upper terrace, so a "higher than most atoms" rule cannot see them.
PARENT_FILE = {"C1-island-near-step": ("Step-8x4.poscar", "step edge of the parent Step-8x4"),
               "Step-8x2_edge-vacancy_plus_foot-adatom": ("Step-8x2.poscar", "straight edge of the parent Step-8x2")}
# the kinks' parent was built inline (an 8x3 slab plus a straight 12-atom strip) and never written to disk. It is
# recovered by the build's own rule: the strip is four full rows of three, plus one atom alone in a fifth row.
PARENT_PARTIAL_ROW = {"Kink-edge1", "Kink-edge2"}
PARENT_EDGE_COLOR = "#6b5b4a"


def partial_row_atoms(at):
    """Feature atoms sitting in a lattice row the rest of the feature leaves under-filled.

    Used only where the build added exactly such an atom. The test is strict: along one lattice direction every
    occupied row but one must hold the same number of atoms, and the odd row must hold fewer. Kink-edge1 gives
    rows {3,3,3,3,1} along a1 and {5,4,4} along a2, so a1 answers and a2 is rejected."""
    pos = at.get_positions(); cell = at.get_cell().array
    top = top_atoms(at); modal = float(np.median(pos[top, 2]))
    idx = np.flatnonzero(top & (pos[:, 2] > modal + 1.0))
    if len(idx) < 4: return np.array([], int)
    C = np.array([cell[0][:2], cell[1][:2]])
    fr = (pos[idx, :2] @ np.linalg.inv(C)) % 1.0
    a0 = layer_spacing(at)
    for ax_ in (0, 1):
        r = fr[:, ax_] * float(np.linalg.norm(cell[ax_][:2]))
        order = np.argsort(r); rows = [[order[0]]]
        for k in order[1:]:
            if r[k] - r[rows[-1][-1]] < 0.4 * a0: rows[-1].append(k)
            else: rows.append([k])
        cnt = np.array([len(g) for g in rows])
        full = int(np.median(cnt))
        odd = [g for g, c in zip(rows, cnt) if c < full]
        if len(odd) == 1 and full >= 2 and sum(1 for c in cnt if c == full) == len(cnt) - 1:
            return idx[np.array(odd[0])]
    return np.array([], int)


def parent_of(sid, at):
    """(parent Atoms in the SAME render frame, added xy, removed xy, caption), or None.

    Split from the overlay drawing because the section line wants to know where the added atoms are before the
    plotting grid exists."""
    cell = at.get_cell().array
    if sid in PARENT_FILE:
        fn, cap = PARENT_FILE[sid]
        par = flatten_cell(read(f"{LIB}/{fn}"))
        Pc = at.get_positions(); Pp = par.get_positions()
        C = np.array([cell[0][:2], cell[1][:2]]); Ci = np.linalg.inv(C)
        f = (Pc[:, None, :2] - Pp[None, :, :2]) @ Ci
        f -= np.round(f)
        d = np.linalg.norm(f @ C, axis=-1) + np.abs(Pc[:, None, 2] - Pp[None, :, 2])
        added = Pc[d.min(axis=1) > 0.8][:, :2]
        removed = Pp[d.min(axis=0) > 0.8][:, :2]
    elif sid in PARENT_PARTIAL_ROW:
        k = partial_row_atoms(at)
        if not len(k): return None
        par = at.copy(); del par[[int(x) for x in k]]
        added = at.get_positions()[k][:, :2]; removed = np.zeros((0, 2))
        cap = "straight edge of the parent strip, without the kink"
    else:
        return None
    return dict(parent=par, added=added, removed=removed, caption=cap)


def parent_overlay(pa, at, gx, gy):
    """The parent's own feature outline, on this structure's grid, to lay over the child as a dashed line."""
    if pa is None: return None
    cell = at.get_cell().array
    hi, lo, _ = feature_sets(pa["parent"])
    pts = hi if len(hi) else lo
    if not len(pts): return None
    D = smooth_periodic(blob_distance(pts, gx, gy, cell))
    px, py, pD = wrap_pad(gx, gy, D, cell)
    return dict(D=D, gx=px, gy=py, Dp=pD, R=0.62 * layer_spacing(pa["parent"]),
                added=pa["added"], removed=pa["removed"], caption=pa["caption"])


def off_canvas(fig, tol=1.0):
    """Labels that run off the canvas. English labels are wider than the CJK ones this layout was first tuned
    for, and a PNG truncates the overflow without complaint, so it has to be checked rather than eyeballed."""
    fig.canvas.draw(); r = fig.canvas.get_renderer()
    fw, fh = fig.canvas.get_width_height()
    out = []
    for t in fig.findobj(matplotlib.text.Text):
        if not t.get_text().strip(): continue
        try: bb = t.get_window_extent(renderer=r)
        except Exception: continue
        if bb.x0 < -tol or bb.y0 < -tol or bb.x1 > fw + tol or bb.y1 > fh + tol:
            out.append(t.get_text().replace("\n", " / ")[:60])
    return out


def render(sid, path, title, meta):
    at = flatten_cell(read(path))
    cell = at.get_cell().array
    pos = at.get_positions()
    cn = coordination(at)
    col = np.array([CN_COLOR[cn_class(c)] for c in cn])
    exp = exposed_atoms(at)
    zt = pos[:, 2].max()

    n1 = max(1, min(4, int(round(TILE_TARGET / np.linalg.norm(cell[0][:2])))))
    n2 = max(1, min(4, int(round(TILE_TARGET / np.linalg.norm(cell[1][:2])))))
    keep = pos[:, 2] > zt - TOP_DEPTH

    P, C, Z, E = [], [], [], []
    for i in range(n1):
        for j in range(n2):
            sh = (i * cell[0] + j * cell[1])[:2]
            P.append(pos[keep][:, :2] + sh); C.append(col[keep]); Z.append(pos[keep][:, 2]); E.append(exp[keep])
    P = np.vstack(P); C = np.concatenate(C); Z = np.concatenate(Z); E = np.concatenate(E)

    pa = parent_of(sid, at)
    sm = surface_model(at, prefer_xy=(pa["added"] if pa else None))
    cuts = sm["cuts"]
    # One side-view panel per section line, each projected in that line's OWN frame: s along A->A', z up, and
    # the perpendicular offset as the depth cue. Projecting on the global x axis whatever the line's direction
    # meant that for a cut along (0.5, 0.866) the "side view along A-A'" was nothing of the sort.
    bands = []
    for ci in (cuts if cuts else [None]):
        if ci is None:
            L = float(np.linalg.norm(cell[0][:2]))
            bands.append(dict(ci=None, mask=np.ones(len(pos), bool), s=pos[:, 0], depth=pos[:, 1], half=None, L=L))
            continue
        L = float(np.linalg.norm(cell[ci["axis"]][:2]))
        # BAND_HALF is a true perpendicular distance now, so it means what the caption says. 0.62 a0 keeps one
        # atom row: the next row along a (111) close-packed direction is a0*sqrt(3)/2 = 2.55 A away, and the old
        # nominal 1.1 a0 pulled in three rows, which stacked three different heights on top of one another.
        half = BAND_HALF_FRAC * layer_spacing(at) if meta.get("family") in BAND_FAMILIES else None
        m, s, dp = cut_band(pos, cell, ci, half if half is not None else 1e9)
        if half is not None and m.sum() < 8:
            half = None; m, s, dp = cut_band(pos, cell, ci, 1e9)
        bands.append(dict(ci=ci, mask=m, s=s, depth=dp, half=half, L=L))

    panels = []
    for bd in bands:
        ns = max(1, min(3, int(round(TILE_TARGET / bd["L"]))))
        m = bd["mask"]
        A = np.vstack([np.c_[bd["s"][m] + i * bd["L"], pos[m, 2]] for i in range(ns)])
        panels.append(dict(P=A, C=np.tile(col[m], ns), Z=np.tile(bd["depth"][m], ns), E=np.tile(exp[m], ns),
                           ns=ns, ci=bd["ci"], half=bd["half"], L=bd["L"]))

    # The marks are settled BEFORE the layout, because the legend's row count sets the bottom strip. It used to
    # be a fixed 0.44 inch, which a two-row legend overran -- the section's A label and its caption ended up
    # printed through the legend text.
    rcls = registry_class(at, meta.get("family"))
    added, missing = added_and_missing(at, meta.get("family"))
    ov = parent_overlay(pa, at, sm["gx"], sm["gy"])
    if ov is not None:                      # build-time provenance beats the height heuristic where both exist
        added = np.array([], int); missing = np.zeros((0, 2))
    hcp_top = np.flatnonzero((rcls == 2) & exp & keep)
    wall_top = np.flatnonzero((rcls == 1) & exp & keep)
    nsurf = int(np.sum(exp & keep))
    n_leg = 4 + (0 if ov is None else 1 + (len(ov["added"]) > 0) + (len(ov["removed"]) > 0)) \
        + (len(added) > 0) + (len(missing) > 0) + (len(hcp_top) > 0) + (len(wall_top) > 0)
    # Columns follow the LONGEST label, not a fixed 4: with wide English labels a four-column legend grew past
    # the 12-inch figure and its first column was pushed off the left edge.
    max_lab = max([24, len(f"hcp-like registry ({len(hcp_top)}/{nsurf})") if len(hcp_top) else 0,
                   len(f"transition / domain wall ({len(wall_top)}/{nsurf})") if len(wall_top) else 0,
                   len(ov["caption"]) if ov is not None else 0,
                   len("CN \u2265 10 (step/island foot, sub-surface)")])
    leg_ncol = min(n_leg, 4 if max_lab <= 34 else 3 if max_lab <= 46 else 2)
    leg_rows = int(np.ceil(n_leg / leg_ncol))

    def ext(A):
        return A[:, 0].max() - A[:, 0].min() + 3 * R_AU, A[:, 1].max() - A[:, 1].min() + 3 * R_AU
    w1, h1 = ext(P)
    ew = [ext(p["P"]) for p in panels]
    w2 = max(e[0] for e in ew); h2 = max(e[1] for e in ew) + 6.0
    nsp = len(panels)
    H = max(h1, nsp * h2)
    FIGW = 12.0
    panel_w = FIGW * 0.97
    scale = panel_w / (w1 + w2 + 1.2)                                # inches per angstrom
    LEG_H = 0.26 + 0.24 * leg_rows                                   # inches actually needed by the legend
    FIGH = H * scale + 1.10 + LEG_H                                  # + title strip and legend strip
    fig = plt.figure(figsize=(FIGW, FIGH), dpi=185)
    top_frac = 1 - 1.02 / FIGH; bot_frac = LEG_H / FIGH
    wsch = 0.52 * (w1 + w2)
    gs = fig.add_gridspec(1, 3, width_ratios=[wsch, w1, w2], wspace=1.2 / (w1 + w2) * 2,
                          left=0.015, right=0.985, top=top_frac, bottom=bot_frac)
    meta["schematic_width_frac"] = float(0.015 + (0.985 - 0.015) * wsch / (wsch + w1 + w2))
    meta["legend_top_frac"] = float(bot_frac)
    meta["schematic_note"] = simple_schematic(fig, gs[0, 0], at, sm, hkl=VICINAL_MILLER.get(sid))
    meta["n_sections"] = len(cuts)
    meta["section_axes"] = [f"a{c['axis'] + 1}" for c in cuts]

    ax = fig.add_subplot(gs[0, 1]); ax.set_aspect("equal"); ax.axis("off")
    draw(ax, P, C, Z, E, R_AU)
    # what was PUT ON and what was TAKEN OUT, and the stacking registry -- none of which colour or height show
    for i in range(n1):
        for j in range(n2):
            sh = (i * cell[0] + j * cell[1])[:2]
            if ov is not None:
                ax.contour(ov["gx"] + sh[0], ov["gy"] + sh[1], ov["Dp"], levels=[ov["R"]],
                           colors=PARENT_EDGE_COLOR, linewidths=1.6, linestyles="dashed", zorder=5900)
                # same array, translated: see tiles()
                for c in ov["added"]:
                    ax.add_patch(Circle(c + sh, R_AU * 1.32, facecolor="none", edgecolor="#1f6f3f",
                                        lw=1.9, zorder=6000))
                for c in ov["removed"]:
                    ax.add_patch(Circle(c + sh, R_AU * 1.05, facecolor="none", edgecolor="#8a3ffc",
                                        lw=1.9, ls=(0, (3, 2)), zorder=6000))
            for k in added:
                ax.add_patch(Circle(pos[k, :2] + sh, R_AU * 1.32, facecolor="none", edgecolor="#1f6f3f",
                                    lw=1.7, zorder=6000))
            for c in missing:
                ax.add_patch(Circle(c + sh, R_AU * 0.95, facecolor="none", edgecolor="#8a3ffc",
                                    lw=1.7, ls=(0, (3, 2)), zorder=6000))
            for k in hcp_top:
                ax.plot(*(pos[k, :2] + sh), marker="x", ms=4.4, mew=1.5, color="#1d4e8f", zorder=6100)
            for k in wall_top:
                ax.plot(*(pos[k, :2] + sh), marker="s", ms=4.0, mew=1.4, mfc="none", color="#b0543a", zorder=6100)
    for i in range(n1):
        for j in range(n2):
            o = i * cell[0][:2] + j * cell[1][:2]
            ax.plot(*zip(o, o + cell[0][:2], o + cell[0][:2] + cell[1][:2], o + cell[1][:2], o),
                    color="#8b9299", lw=0.7, ls=(0, (4, 3)), zorder=5000)
    ax.set_xlim(P[:, 0].min() - 1.5 * R_AU, P[:, 0].max() + 1.5 * R_AU)
    ax.set_ylim(P[:, 1].min() - 1.5 * R_AU, P[:, 1].max() + 1.5 * R_AU)
    ax.set_title(f"top view \u00b7 {n1}\u00d7{n2} cells \u00b7 top {TOP_DEPTH:.1f} " + AA,
                 fontproperties=CJK, fontsize=9.5, color="#4a5158", pad=3)

    import matplotlib.gridspec as mgs
    side = mgs.GridSpecFromSubplotSpec(nsp, 1, subplot_spec=gs[0, 2], hspace=0.30)
    for k, pn in enumerate(panels):
        ax2 = fig.add_subplot(side[k]); ax2.set_aspect("equal"); ax2.axis("off")
        draw(ax2, pn["P"], pn["C"], pn["Z"], pn["E"], R_AU)
        # see the note in section_panel: no fabricated accessibility line, only a side label
        ax2.text(pn["P"][:, 0].min() - R_AU, zt + 3.2, "\u2191 electrolyte side",
                 fontproperties=CJK, fontsize=8.6, color="#2e7d9a")
        if k == nsp - 1:
            ax2.text(pn["P"][:, 0].min() - R_AU, zt + 1.6, "(the fixed back layers are below)",
                     fontproperties=CJK, fontsize=7.4, color="#9aa1a8", va="bottom", ha="left")
        ax2.set_xlim(pn["P"][:, 0].min() - 1.5 * R_AU, pn["P"][:, 0].max() + 1.5 * R_AU)
        ax2.set_ylim(pn["P"][:, 1].min() - 1.5 * R_AU, zt + 7.0)
        if pn["ci"] is None:
            t_ = f"side view \u00b7 along a1 \u00b7 {pn['ns']} cells"
        else:
            lb = pn["ci"]["label"]; ab = f"{lb}–{lb}′"
            t_ = (f"side view \u00b7 \u00b1{pn['half']:.1f} " + AA + f" band about {ab}" if pn["half"]
                  else f"side view \u00b7 along {ab} \u00b7 {pn['ns']}\u00d7")
        ax2.set_title(t_, fontproperties=CJK, fontsize=9.3, color="#4a5158", pad=3)

    fig.text(0.5, 1 - 0.30 / FIGH, title, ha="center", va="top", fontsize=13.5, color="#14181c", weight="medium")
    sub = (f"{meta['n_atoms']} Au   \u00b7   area {meta['A_proj']:.0f} " + AA + "$^2$"
           f"   ·   {meta['family_zh']}   ·   {meta['cn_counts_zh']}")
    fig.text(0.5, 1 - 0.70 / FIGH, sub, ha="center", va="top", fontproperties=CJK, fontsize=9.2, color="#5a616a")
    h = [plt.Line2D([], [], marker="o", ls="", ms=7, mfc=CN_COLOR[k], mec="#2b2f36", mew=0.5, label=l)
         for k, l in [("kink", "CN \u2264 6 (kinks, adatoms)"), ("edge", "CN 7\u20138 (step, island, pit edges)"),
                      ("terrace", "CN 9 (flat terrace)"), ("bulk", "CN \u2265 10 (step/island foot, sub-surface)")]]
    if ov is not None:
        h.append(plt.Line2D([], [], ls=(0, (5, 3)), lw=1.6, color=PARENT_EDGE_COLOR, label=ov["caption"]))
        if len(ov["added"]): h.append(plt.Line2D([], [], marker="o", ls="", ms=9, mfc="none", mec="#1f6f3f",
                                                 mew=1.9, label=f"Au added vs the parent ({len(ov['added'])})"))
        if len(ov["removed"]): h.append(plt.Line2D([], [], marker="o", ls="", ms=8, mfc="none", mec="#8a3ffc",
                                                   mew=1.9, label=f"site removed from the parent ({len(ov['removed'])})"))
    if len(added): h.append(plt.Line2D([], [], marker="o", ls="", ms=9, mfc="none", mec="#1f6f3f", mew=1.7,
                                       label=f"added Au ({len(added)})"))
    if len(missing): h.append(plt.Line2D([], [], marker="o", ls="", ms=8, mfc="none", mec="#8a3ffc", mew=1.7,
                                         label=f"removed site ({len(missing)})"))
    if len(hcp_top): h.append(plt.Line2D([], [], marker="x", ls="", ms=6, mew=1.5, color="#1d4e8f",
                                         label=f"hcp-like registry ({len(hcp_top)}/{nsurf})"))
    if len(wall_top): h.append(plt.Line2D([], [], marker="s", ls="", ms=6, mew=1.4, mfc="none", color="#b0543a",
                                          label=f"transition / domain wall ({len(wall_top)}/{nsurf})"))
    assert len(h) == n_leg, (len(h), n_leg)        # the layout reserved space for exactly this many entries
    lg = fig.legend(handles=h, loc="lower center", ncol=leg_ncol, fontsize=8.2, frameon=False,
                    bbox_to_anchor=(0.5, 0.04 / FIGH), prop=CJK)
    for t in lg.get_texts(): t.set_fontproperties(CJK); t.set_fontsize(8.4)
    meta["registry_display_classes"] = dict(hcp_like=int(len(hcp_top)), transition=int(len(wall_top)),
                                            fcc_like=int(nsurf - len(hcp_top) - len(wall_top)), n_surface=nsurf)
    over = off_canvas(fig)
    fig.savefig(f"{OUT}/{sid}.png", facecolor="white")
    plt.close(fig)
    meta["off_canvas"] = over
    return {c: int((np.array([cn_class(x) for x in cn]) == c).sum()) for c in CN_COLOR}
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--only", default=None); a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    ch = json.load(open(f"{ROOT}/analysis/charging/charging.json"))["geometries"]
    S = json.load(open(f"{ROOT}/dataset_v1/states.json"))["states"]
    # one entry per structure: prefer its 'ideal' geometry, else the relaxed one
    best = {}
    for s in S.values():
        if s["campaign"] != "dataset_plan_v1_rev2": continue
        cfg = "relaxed" if s["config"] in ("relax", "relaxed") else s["config"]
        if cfg not in ("ideal", "relaxed"): continue
        cur = best.get(s["structure_id"])
        if cur is None or (cur["cfg"] != "ideal" and cfg == "ideal"):
            best[s["structure_id"]] = dict(cfg=cfg, dir=s["source_dir"], geom=s["geometry_file"],
                                           family=s["family"], tier=s["tier"], n_atoms=s["structure"]["n_atoms"],
                                           geometry_id=s["geometry_id"])
    out = {}
    for sid, m in sorted(best.items()):
        if a.only and a.only != sid: continue
        g = ch.get(m["geometry_id"], {})
        m["A_proj"] = g.get("A_proj", 0.0)
        at = read(f"{m['dir']}/{m['geom']}")
        cn = coordination(at)
        cc = collections_counter(cn)
        m["cn_counts"] = ", ".join(f"{k} {v}" for k, v in cc.items() if v)
        # The whole-slab CN histogram counts the BOTTOM surface too: a 4-layer flat slab has 32 CN-9 atoms, but
        # only 16 of them face the electrolyte. Report the upper surface separately so the caption cannot be read
        # as "the electrolyte touches 32 terrace atoms".
        at_r = flatten_cell(at)
        up = top_atoms(at_r); cn_up = coordination(at_r)[up]
        cu = {"kink": 0, "edge": 0, "terrace": 0, "bulk": 0}
        for x in cn_up: cu[cn_class(x)] += 1
        m["cn_counts_upper"] = cu
        m["n_upper_surface"] = int(up.sum())
        m["cn_counts_zh"] = ("upper surface " + str(int(up.sum())) + " atoms:  " +
                             "  ".join(f"{CN_ZH[k]} {v}" for k, v in cu.items() if v))
        m["family_zh"] = FAMILY_ZH.get(m["family"], m["family"])
        m["cn_counts_dict"] = cc
        render(sid, f"{m['dir']}/{m['geom']}", sid, m)
        pz = g.get("pzc", {})
        m.update(U_pzc_V=pz.get("U_pzc"), pzc_bracketed=pz.get("bracketed"),
                 C_median_uF_cm2=float(np.median([s["C_uF_per_cm2"] for s in g["secants"]])) if g.get("secants") else None,
                 dU_pzc_vs_cell_ref_mV=g.get("dU_pzc_vs_cell_ref_mV"), dC_vs_cell_ref_pct=g.get("dC_vs_cell_ref_pct"),
                 png=f"{sid}.png")
        out[sid] = m
        print(f"  {sid:38s} {m['n_atoms']:4d} atoms  {m['cn_counts']}")
    old = json.load(open(f"{OUT}/gallery.json")) if os.path.exists(f"{OUT}/gallery.json") else {}
    old.update(out); json.dump(old, open(f"{OUT}/gallery.json", "w"), indent=1)
    bad = {k: v["off_canvas"] for k, v in out.items() if v.get("off_canvas")}
    print(f"{len(out)} structures rendered -> {OUT}/")
    if bad:
        print(f"LABELS RUNNING OFF THE CANVAS in {len(bad)} figures:")
        for k, v in sorted(bad.items()):
            for t in v: print(f"   {k:34s} {t!r}")
        raise SystemExit(1)


def collections_counter(cn):
    c = {"kink": 0, "edge": 0, "terrace": 0, "bulk": 0}
    for x in cn: c[cn_class(x)] += 1
    return c


if __name__ == "__main__":
    main()
