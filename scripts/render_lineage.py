#!/usr/bin/env python3
"""How one built geometry becomes several computed configurations -- the opening figure of the report.

It answers one question and no more: the 107 geometries are not 107 independently designed defects. 33 are
built by hand, and every other geometry is a relaxed, perturbed, collectively deformed or path image of one of
those 33. The structure families are not drawn here; they are in the gallery.

The footer names the FIVE potentials explicitly, because "constant potential" is otherwise easy to read as one
setting: every geometry is computed at a fixed electron chemical potential, three of them complete and two
still running, and the electron count is solved for rather than set.

Deliberately sparse, so the type can be large. The canvas is laid out in POINTS, so a label's size in points
is its size in layout units, and check_layout() asserts on the real rendered extents that nothing overlaps and
nothing runs off the canvas -- the two ways a figure like this silently goes wrong.

Usage (from Au_Cl/):  scripts/pyrun.sh scripts/render_lineage.py
Output: analysis/gallery/_lineage.png
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
OUT = f"{ROOT}/analysis/gallery/_lineage.png"
FONT = matplotlib.font_manager.FontProperties(family="DejaVu Sans")

INK, SUB, FAINT, LINE = "#14181c", "#5a616a", "#9aa1a8", "#8b9299"
GOLD, BLUE = "#c9a961", "#7fa8bd"
W, H = 1360.0, 782.0

MU0 = -4.9071
# (TARGETMU eV, U V, label, complete?)
POTENTIALS = [(-5.4071, +0.5, False), (-5.1071, +0.2, True), (-4.9071, 0.0, True),
              (-4.7071, -0.2, True), (-4.4071, -0.5, False)]


def rbox(ax, x, cy, w, h, fc, ec, lw=1.8):
    ax.add_patch(FancyBboxPatch((x, cy - h / 2), w, h, boxstyle="round,pad=0,rounding_size=9",
                                facecolor=fc, edgecolor=ec, linewidth=lw, zorder=2))


def tx(ax, x, y, s, size, color=INK, weight="normal", ha="left", store=None, style="normal"):
    t = ax.text(x, y, s, fontproperties=FONT, fontsize=size, color=color, fontweight=weight,
                ha=ha, va="center", zorder=5, style=style)
    if store is not None: store.append(t)
    return t


def node(ax, x, cy, w, h, title, count, note, fc, ec, T, ts=19.0, cs=15.0, ns=13.5):
    rbox(ax, x, cy, w, h, fc, ec)
    tx(ax, x + 18, cy - h / 2 + 24, title, ts, INK, "semibold", store=T)
    tx(ax, x + 18, cy - h / 2 + 50, count, cs, SUB, store=T)
    if note: tx(ax, x + 18, cy - h / 2 + 74, note, ns, FAINT, store=T)


def elbow(ax, x0, y0, x1, y1, color=LINE, lw=2.0):
    xm = x0 + (x1 - x0) * 0.45
    ax.plot([x0, xm, xm, x1], [y0, y0, y1, y1], color=color, lw=lw, zorder=1,
            solid_capstyle="round", solid_joinstyle="round")
    ax.annotate("", xy=(x1, y1), xytext=(x1 - 14, y1), zorder=3,
                arrowprops=dict(arrowstyle="-|>", lw=lw, color=color, shrinkA=0, shrinkB=0))


def check_layout(fig, texts, tol=1.0):
    """No label may overlap another or leave the canvas."""
    fig.canvas.draw(); r = fig.canvas.get_renderer()
    bb = [(t, t.get_window_extent(renderer=r)) for t in texts]
    fw, fh = fig.canvas.get_width_height()
    bad = []
    for t, a in bb:
        if a.x0 < -tol or a.y0 < -tol or a.x1 > fw + tol or a.y1 > fh + tol:
            bad.append((f"OFF-CANVAS {t.get_text()[:44]!r}", ""))
    for i in range(len(bb)):
        for j in range(i + 1, len(bb)):
            a, b = bb[i][1], bb[j][1]
            if min(a.x1, b.x1) - max(a.x0, b.x0) > tol and min(a.y1, b.y1) - max(a.y0, b.y0) > tol:
                bad.append((bb[i][0].get_text()[:34], bb[j][0].get_text()[:34]))
    return bad


def main():
    DPI = 150
    fig = plt.figure(figsize=(W / 72.0, H / 72.0), dpi=DPI)
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.axis("off")
    ax.add_patch(plt.Rectangle((0, 0), W, H, facecolor="white", edgecolor="none", zorder=0))
    T = []

    tx(ax, 24, 40, "How one built geometry becomes many computed states", 26, INK, "semibold", store=T)
    tx(ax, 24, 74, "33 structures are built by hand; every other geometry is a relaxed, perturbed or path "
                   "image of one of them.", 15, SUB, store=T)
    tx(ax, 24, 98, "107 geometries in all, not 107 separately designed defects.", 15, SUB, store=T)

    X0, W0 = 24.0, 268.0
    X1, W1 = 370.0, 316.0
    X2, W2 = 790.0, 540.0
    yA, yB = 206.0, 416.0
    yRoot = 0.5 * (yA + yB)
    y3 = [316.0, 408.0, 500.0]

    node(ax, X0, yRoot, W0, 96, "Hand-built geometry", "33", "one per structure; see the gallery",
         "#ffffff", INK, T, ts=19, cs=16)

    elbow(ax, X0 + W0, yRoot, X1 - 8, yA)
    node(ax, X1, yA, W1, 96, "Ideal  ·  ideal", "33 geometries", "static reference; never relaxed",
         "#f3efe7", GOLD, T)

    elbow(ax, X0 + W0, yRoot, X1 - 8, yB)
    node(ax, X1, yB, W1, 96, "Relaxed at the reference $\\mu_0$", "16 geometries · relax / relaxed",
         "local reference state for each structure", "#f3efe7", GOLD, T)

    kids = [("Random displacement  ·  pert05 / pert10", "16 + 16 geometries",
             "Gaussian 0.05 / 0.10 Å on the movable atoms"),
            ("Collective deformation  ·  coll", "14 geometries",
             "top-layer spacing −3%, in-plane strain +1%, edge bend"),
            ("Path image  ·  path", "12 geometries",
             "one atom moved along a hop; linear images, not NEB")]
    for (t_, c_, n_), y in zip(kids, y3):
        elbow(ax, X1 + W1, yB, X2 - 8, y)
        node(ax, X2, y, W2, 88, t_, c_, n_, "#eef2f5", BLUE, T, ts=17.5, cs=14.5, ns=13)

    # ---- the five potentials, spelled out ----
    TOP, BH = 566.0, 190.0
    rbox(ax, 24, TOP + BH / 2, W - 48, BH, "#f2f6f8", "#9dbccb", 1.6)
    tx(ax, 44, TOP + 24, "Every geometry is then computed at FIVE fixed potentials", 17, INK,
       "semibold", store=T)
    tx(ax, 44, TOP + 48, r"constant-$\mu$ DFT: the electron chemical potential is set, the electron count is "
                         r"solved for.   U = $\mu_0-\mu_e$,  $\mu_0$ = $-$4.9071 eV", 13.5, SUB, store=T)

    xL, xR = 150.0, W - 150.0
    ymark = TOP + 98
    ax.plot([xL - 46, xR + 46], [ymark, ymark], color="#9dbccb", lw=1.6, zorder=2)
    for k, (mu, u, done) in enumerate(POTENTIALS):
        x = xL + (xR - xL) * k / (len(POTENTIALS) - 1)
        ax.plot([x], [ymark], marker="o", ms=15 if done else 14,
                mfc=("#2e7d9a" if done else "#ffffff"), mec="#2e7d9a", mew=2.2, zorder=4)
        tx(ax, x, ymark - 24, f"U = {u:+.1f} V".replace("+0.0", " 0.0").replace("-", "\u2212"), 15,
           INK if done else SUB, "semibold", "center", store=T)
        tx(ax, x, ymark + 26, f"TARGETMU {mu:.4f}".replace("-", "−"), 12.5, SUB, ha="center", store=T)
    tx(ax, 44, TOP + 150, "filled = the three base potentials, complete (291 electronic states). "
                          "Every result on this page uses only these three.", 13, INK, store=T)
    tx(ax, 44, TOP + 172, "open = the ±0.5 V extension, still running. Drawn as a faint overlay on the σ(U) "
                          "chart; it enters no number here.", 13, SUB, store=T)

    bad = check_layout(fig, T)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, facecolor="white")
    plt.close(fig)
    if bad:
        print(f"LAYOUT PROBLEMS ({len(bad)}):")
        for a, b in bad[:10]: print(f"   {a!r} x {b!r}")
        raise SystemExit(1)
    print(f"{os.path.getsize(OUT)/1024:.0f} kB -> {OUT}  ({int(W/72*DPI)}x{int(H/72*DPI)} px, "
          f"canvas {W/72:.1f}x{H/72:.1f} in, {len(T)} labels, none overlapping, none off-canvas)")


if __name__ == "__main__":
    main()
