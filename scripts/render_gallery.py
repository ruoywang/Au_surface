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
# Droid Sans Fallback is the only CJK face on this machine and matplotlib 3.5 has no per-glyph fallback, so the
# Angstrom sign must come from mathtext ($\AA$), which is typeset with the math font, not the text font.
CJK = matplotlib.font_manager.FontProperties(fname="/usr/share/fonts/google-droid/DroidSansFallback.ttf")
AA = r"$\mathrm{\AA}$"


FAMILY_ZH = {"flat Au(111)": "平整 Au(111)", "point defect": "点缺陷", "reconstruction-related": "重构相关",
             "strip step": "条带台阶", "vicinal step face": "邻晶面台阶", "kink / edge rearrangement": "拐角 / 边缘重排",
             "single-layer island": "单层岛", "single-layer pit": "单层坑", "composite": "复合形貌"}
CN_ZH = {"kink": "拐角", "edge": "边缘", "terrace": "平台", "bulk": "次表面"}


def cn_class(c):
    return "kink" if c <= 6 else "edge" if c <= 8 else "terrace" if c == 9 else "bulk"


SURFACE_BAND = 3.0      # A below the highest atom; the SAME surface set analysis_spatial.column_labels uses


def exposed_atoms(atoms):
    """The surface set: atoms within SURFACE_BAND of the highest one.

    This is exactly the set the region analysis assigns columns to, so the emphasised atoms in the picture are the
    atoms that carry a coordination label in the anion maps. It keeps BOTH terraces of a step (2.4 A apart) and the
    step-foot row, which is over-coordinated but still faces the electrolyte; an absolute depth fade would wrongly
    bury half of every stepped surface."""
    z = atoms.get_positions()[:, 2]
    return z > z.max() - SURFACE_BAND


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


def describe(H):
    """One Chinese line saying what the outline shows, generated from the height field itself."""
    rng = float(H.max() - H.min())
    if rng < 0.35:
        return "高度均匀，无起伏（若与平板有别，差别在层序或面内配准，不在高度）"
    gap, cut = height_gap(H)
    if gap < GAP_LAYER:
        return f"连续倾斜的平台，胞内高差 {rng:.1f} Å 且无断层——邻晶面的微斜切割，台阶并到胞边"
    v = np.sort(H.ravel())
    cuts = [0.5 * (a + b) for a, b in zip(v[:-1], v[1:]) if b - a >= GAP_LAYER]
    if len(cuts) >= 2:                                     # three or more separated levels
        edges = [-np.inf] + cuts + [np.inf]
        fr = [float(((H > lo) & (H <= hi)).mean()) for lo, hi in zip(edges[:-1], edges[1:])]
        return f"{len(fr)} 个分离的高度层，自上而下占 " + "、".join(f"{100*x:.0f}%" for x in reversed(fr))
    up = float((H > cut).mean()); dn = 1 - up
    if 0.30 <= up <= 0.70:
        return f"两层平台，上 {100*up:.0f}%、下 {100*dn:.0f}%（台阶），层间 {gap:.1f} Å"
    if up < dn:
        return f"高出一层的凸起，占面积 {100*up:.0f}%，高 {gap:.1f} Å"
    return f"低一层的凹陷，占面积 {100*dn:.0f}%，深 {gap:.1f} Å"


PLAN_FILL = {2: "#dfb264", 1: "#e9c98f", 0: "#f2efe9", -1: "#bcd3dd", -2: "#9cbecd", -3: "#86adbf"}


def wrap_offset(L):
    """Always zero. Kept as a named function so the reason is recorded rather than rediscovered.

    An earlier version measured a net level change across the periodic seam and tiled the schematic with it, to draw
    a vicinal face as a descending staircase. That is wrong here: every cell in this set has a1_z = a2_z = 0, so the
    surface must return to the same height after one cell by construction, and there is no net descent to add back.
    The measured 'offset' was just the first and last tenth of the cell landing on different terraces, which also
    mislabelled every strip-step cell as a vicinal."""
    return 0, 0


def simple_schematic(fig, cell_spec, at):
    """The plain outline: no atoms at all. Top = plan view of the surface levels with their boundaries stroked,
    bottom = the side silhouette. Both come from the real height map, so the cartoon cannot drift from the geometry."""
    import matplotlib.gridspec as mgs
    H, gx, gy, _ = height_map(at, ng=110)
    H = H - np.median(H)
    cell = at.get_cell().array
    inner = mgs.GridSpecFromSubplotSpec(2, 1, subplot_spec=cell_spec, height_ratios=[2.35, 1.0], hspace=0.34)

    # Smooth before contouring so the outline reads as one clean curve instead of tracing the atomic scallops of the
    # height map. Periodic Gaussian (via FFT), sigma ~1.2 A -- well under a nearest-neighbour spacing, so the feature
    # keeps its real size and shape; it only removes the per-atom ripple on the boundary.
    # sigma is fixed in GRID cells, not in angstrom: a vicinal face has terrace stripes barely 1 A wide, and a fixed
    # 0.8 A kernel flattened them below the +-0.5 contour level, leaving a blank panel. 1.5 cells always rounds the
    # per-atom scallops without ever erasing a feature the grid can resolve.
    n1g, n2g = H.shape
    k1 = np.fft.fftfreq(n1g)[:, None]; k2 = np.fft.fftfreq(n2g)[None, :]
    Hs = np.fft.ifft2(np.fft.fft2(H) * np.exp(-2 * (np.pi ** 2) * (1.5 ** 2) * (k1 ** 2 + k2 ** 2))).real

    # --- plan: one flat tone per level, one stroke on each boundary, the cell outline dashed. Nothing else.
    # A cell whose surface steps across the periodic seam (a vicinal terrace) is tiled with that offset added back,
    # so the staircase appears instead of one uniform block.
    t1 = max(1, min(3, int(round(TILE_TARGET / np.linalg.norm(cell[0][:2])))))
    t2 = max(1, min(3, int(round(TILE_TARGET / np.linalg.norm(cell[1][:2])))))
    axp = fig.add_subplot(inner[0]); axp.set_aspect("equal"); axp.axis("off")
    # A simple drawing needs FEW strokes. Cut the height field only where it is genuinely discontinuous -- one bold
    # outline per separated atomic level. A vicinal face has no such discontinuity inside the cell (its terrace is an
    # inclined plane), so it gets a smooth shade and a downhill arrow instead of a fan of meaningless contour bands.
    v = np.sort(H.ravel())
    cuts = [float(0.5 * (a + b)) for a, b in zip(v[:-1], v[1:]) if b - a >= GAP_LAYER]
    cmap = plt.get_cmap("RdYlBu_r")
    if cuts:
        edges = [Hs.min() - 1] + cuts + [Hs.max() + 1]
        mids = [0.5 * (edges[k] + edges[k + 1]) for k in range(len(edges) - 1)]
        span = max(max(abs(m) for m in mids), 1.2)
        cols = [cmap(0.5 + 0.40 * m / span) for m in mids]
        for i in range(t1):
            for j in range(t2):
                sh = (i * cell[0] + j * cell[1])[:2]
                axp.contourf(gx + sh[0], gy + sh[1], Hs, levels=edges, colors=cols, zorder=1)
                axp.contour(gx + sh[0], gy + sh[1], Hs, levels=cuts, colors="#23272c",
                            linewidths=2.2, linestyles="solid", zorder=3)
    else:
        rng = max(Hs.ptp(), 1e-6)
        for i in range(t1):
            for j in range(t2):
                sh = (i * cell[0] + j * cell[1])[:2]
                axp.pcolormesh(gx + sh[0], gy + sh[1], Hs, cmap=cmap, vmin=-0.6 * rng, vmax=0.6 * rng,
                               shading="gouraud", zorder=1)
        if Hs.ptp() > 0.35:                                        # mark which way the inclined terrace runs downhill
            g1 = float(np.mean(np.gradient(Hs, axis=0))); g2 = float(np.mean(np.gradient(Hs, axis=1)))
            d = -(g1 * cell[0][:2] / np.linalg.norm(cell[0][:2]) + g2 * cell[1][:2] / np.linalg.norm(cell[1][:2]))
            if np.linalg.norm(d) > 0:
                d = d / np.linalg.norm(d) * 0.22 * float(np.ptp(gx))
                c0 = np.array([gx.mean() * t1, gy.mean() * t2])
                axp.annotate("", xy=c0 + d, xytext=c0 - d, zorder=4,
                             arrowprops=dict(arrowstyle="-|>", lw=2.0, color="#23272c"))
                axp.text(*(c0 + 1.25 * d), "下坡", fontproperties=CJK, fontsize=8.5, color="#23272c",
                         ha="center", va="center", zorder=5)
    o = np.zeros(2)
    axp.plot(*zip(o, cell[0][:2], cell[0][:2] + cell[1][:2], cell[1][:2], o),
             color="#9aa1a8", lw=0.9, ls=(0, (4, 3)), zorder=4)
    X = np.concatenate([(gx + (i * cell[0] + j * cell[1])[0]).ravel() for i in range(t1) for j in range(t2)])
    Y = np.concatenate([(gy + (i * cell[0] + j * cell[1])[1]).ravel() for i in range(t1) for j in range(t2)])
    axp.set_xlim(X.min(), X.max()); axp.set_ylim(Y.min(), Y.max())
    axp.set_title(f"简笔示意 · 俯视轮廓（{t1}×{t2} 个胞）", fontproperties=CJK, fontsize=9.5, color="#4a5158", pad=3)

    # --- section: a real CUT through the middle of the feature, not a projection. Projecting the maximum would turn
    # a compact island into a plateau as wide as the island's whole footprint, i.e. make it look like a step.
    axs = fig.add_subplot(inner[1]); axs.axis("off")

    def circ_centre(mask_1d):
        """Index of the circular mean of a periodic boolean mask (a feature may straddle the cell boundary)."""
        n = len(mask_1d); idx = np.flatnonzero(mask_1d)
        if not len(idx): return n // 2
        a = 2 * np.pi * idx / n
        return int(round((np.arctan2(np.sin(a).mean(), np.cos(a).mean()) % (2 * np.pi)) * n / (2 * np.pi))) % n

    gap, cut = height_gap(H)
    M = H > cut if gap >= GAP_LAYER else np.zeros(H.shape, bool)
    if M.any() and not M.all():
        minority = M if M.mean() <= 0.5 else ~M
        i0 = circ_centre(minority.any(axis=1)); j0 = circ_centre(minority.any(axis=0))
        use_a1 = minority.any(axis=1).mean() <= minority.any(axis=0).mean()
        prof = H[:, j0].copy() if use_a1 else H[i0, :].copy()
        if prof.ptp() < 0.5 * gap:                                # the centre line misses it: use the silhouette
            prof = H.max(axis=1) if use_a1 else H.max(axis=0)
        c = circ_centre((prof > cut) if minority is M else (prof <= cut))
        prof = np.roll(prof, len(prof) // 2 - c)
    else:                                                          # flat, or a continuous ramp with no separated level
        use_a1 = H.max(axis=1).ptp() >= H.max(axis=0).ptp()
        prof = H.max(axis=1) if use_a1 else H.max(axis=0)
    length = np.linalg.norm(cell[0][:2] if use_a1 else cell[1][:2])
    reps = max(1, min(2, int(round(TILE_TARGET / length))))
    prof = np.tile(prof, reps)
    s = np.append(np.linspace(0, reps * length, len(prof), endpoint=False), reps * length)
    p = np.append(prof, prof[0]) - prof.min()
    body = 4.8
    axs.fill_between(s, p, -body, step="mid", color="#e6e2da", zorder=1)
    axs.step(s, p, where="mid", color="#2b3137", lw=2.1, zorder=3)
    axs.plot([0, reps * length], [-body] * 2, color="#c6c1b7", lw=1.0, zorder=3)
    top = p.max()
    axs.axhline(top + 3.2, color="#2e7d9a", lw=1.1, ls=(0, (5, 3)), zorder=2)
    axs.text(reps * length * 0.01, top + 3.7, "电解质", fontproperties=CJK, fontsize=8.2, color="#2e7d9a", va="bottom")
    axs.set_xlim(0, reps * length); axs.set_ylim(-body - 0.4, top + 6.2)
    axs.set_title(f"侧面剪影（{reps} 个周期，高度按真实比例）", fontproperties=CJK, fontsize=9.5, color="#4a5158", pad=2)
    return describe(H)


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

    ns = max(1, min(3, int(round(TILE_TARGET / np.linalg.norm(cell[0][:2])))))
    P2, C2, Z2, E2 = [], [], [], []
    for i in range(ns):
        P2.append(np.c_[pos[:, 0] + i * cell[0][0], pos[:, 2]]); C2.append(col); Z2.append(pos[:, 1])
        E2.append(exp)                                            # side view: same exposure flag, depth cue along y
    P2 = np.vstack(P2); C2 = np.concatenate(C2); Z2 = np.concatenate(Z2); E2 = np.concatenate(E2)

    def ext(A):
        return A[:, 0].max() - A[:, 0].min() + 3 * R_AU, A[:, 1].max() - A[:, 1].min() + 3 * R_AU
    w1, h1 = ext(P); w2, h2 = ext(P2); h2 += 6.0                     # room for the electrolyte marker
    H = max(h1, h2)
    FIGW = 12.0
    panel_w = FIGW * 0.97
    scale = panel_w / (w1 + w2 + 1.2)                                # inches per angstrom
    FIGH = H * scale + 1.55                                          # + title strip and legend strip
    fig = plt.figure(figsize=(FIGW, FIGH), dpi=185)
    top_frac = 1 - 1.02 / FIGH; bot_frac = 0.44 / FIGH
    gs = fig.add_gridspec(1, 3, width_ratios=[0.52 * (w1 + w2), w1, w2], wspace=1.2 / (w1 + w2) * 2,
                          left=0.015, right=0.985, top=top_frac, bottom=bot_frac)
    meta["schematic_note"] = simple_schematic(fig, gs[0, 0], at)

    ax = fig.add_subplot(gs[0, 1]); ax.set_aspect("equal"); ax.axis("off")
    draw(ax, P, C, Z, E, R_AU)
    for i in range(n1):
        for j in range(n2):
            o = i * cell[0][:2] + j * cell[1][:2]
            ax.plot(*zip(o, o + cell[0][:2], o + cell[0][:2] + cell[1][:2], o + cell[1][:2], o),
                    color="#8b9299", lw=0.7, ls=(0, (4, 3)), zorder=5000)
    ax.set_xlim(P[:, 0].min() - 1.5 * R_AU, P[:, 0].max() + 1.5 * R_AU)
    ax.set_ylim(P[:, 1].min() - 1.5 * R_AU, P[:, 1].max() + 1.5 * R_AU)
    ax.set_title(f"俯视图 · {n1}×{n2} 个胞 · 最外 {TOP_DEPTH:.1f} " + AA,
                 fontproperties=CJK, fontsize=9.5, color="#4a5158", pad=3)

    ax2 = fig.add_subplot(gs[0, 2]); ax2.set_aspect("equal"); ax2.axis("off")
    draw(ax2, P2, C2, Z2, E2, R_AU)
    ax2.axhline(zt + 4.2, color="#2e7d9a", lw=1.1, ls=(0, (5, 3)), zorder=6000)
    ax2.text(P2[:, 0].min() - R_AU, zt + 4.8, "此线以上为离子可达的电解质",
             fontproperties=CJK, fontsize=8.2, color="#2e7d9a")
    ax2.set_xlim(P2[:, 0].min() - 1.5 * R_AU, P2[:, 0].max() + 1.5 * R_AU)
    ax2.set_ylim(P2[:, 1].min() - 1.5 * R_AU, zt + 7.0)
    ax2.set_title(f"侧视图 · 沿 a1 方向 {ns} 个胞 · 全部原子层",
                  fontproperties=CJK, fontsize=9.5, color="#4a5158", pad=3)

    fig.text(0.5, 1 - 0.30 / FIGH, title, ha="center", va="top", fontsize=13.5, color="#14181c", weight="medium")
    sub = (f"{meta['n_atoms']} 个 Au   ·   投影面积 {meta['A_proj']:.0f} " + AA + "$^2$"
           f"   ·   {meta['family_zh']}   ·   {meta['cn_counts_zh']}")
    fig.text(0.5, 1 - 0.70 / FIGH, sub, ha="center", va="top", fontproperties=CJK, fontsize=9.2, color="#5a616a")
    h = [plt.Line2D([], [], marker="o", ls="", ms=7, mfc=CN_COLOR[k], mec="#2b2f36", mew=0.5, label=l)
         for k, l in [("kink", "CN≤6  拐角 / 吸附原子"), ("edge", "CN 7–8  台阶边、岛与坑的边缘"),
                      ("terrace", "CN 9  平整平台"), ("bulk", "CN≥10  台阶脚 / 岛脚 / 次表面")]]
    lg = fig.legend(handles=h, loc="lower center", ncol=4, fontsize=8.4, frameon=False,
                    bbox_to_anchor=(0.5, 0.02 / FIGH), prop=CJK)
    for t in lg.get_texts(): t.set_fontproperties(CJK); t.set_fontsize(8.4)
    fig.savefig(f"{OUT}/{sid}.png", facecolor="white")
    plt.close(fig)
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
        m["cn_counts_zh"] = "  ".join(f"{CN_ZH[k]} {v}" for k, v in cc.items() if v)
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
    print(f"{len(out)} structures rendered -> {OUT}/")


def collections_counter(cn):
    c = {"kink": 0, "edge": 0, "terrace": 0, "bulk": 0}
    for x in cn: c[cn_class(x)] += 1
    return c


if __name__ == "__main__":
    main()
