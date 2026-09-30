#!/usr/bin/env python3
"""One built geometry becomes several computed configurations -- the opening figure of the report.

It answers one question and no more: the 107 geometries are not 107 independently designed defects. 33 are
built by hand, and every other geometry is a relaxed, perturbed, collectively deformed or path image of one of
those 33. The structure families are not drawn here; they are in the gallery.

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
CJK = matplotlib.font_manager.FontProperties(fname="/usr/share/fonts/google-droid/DroidSansFallback.ttf")
AA = r"$\mathrm{\AA}$"          # Droid Sans Fallback has no U+00C5 / U+2212 / U+2080; those come from mathtext

INK, SUB, FAINT, LINE = "#14181c", "#5a616a", "#9aa1a8", "#8b9299"
W, H = 1120.0, 664.0


def rbox(ax, x, cy, w, h, fc, ec, lw=1.8):
    ax.add_patch(FancyBboxPatch((x, cy - h / 2), w, h, boxstyle="round,pad=0,rounding_size=9",
                                facecolor=fc, edgecolor=ec, linewidth=lw, zorder=2))


def tx(ax, x, y, s, size, color=INK, weight="normal", ha="left", store=None):
    t = ax.text(x, y, s, fontproperties=CJK, fontsize=size, color=color, fontweight=weight,
                ha=ha, va="center", zorder=5)
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

    tx(ax, 24, 40, "一个建构好的几何，如何变成多个被计算的构型", 26, INK, "semibold", store=T)
    tx(ax, 24, 76, "33 个结构是人工建构的；其余几何都是它们的弛豫、扰动与路径像。合计 107 个几何，"
                   "不是 107 个独立设计的缺陷。", 15, SUB, store=T)

    X0, W0 = 24.0, 254.0
    X1, W1 = 356.0, 306.0
    X2, W2 = 740.0, 356.0
    yA, yB = 186.0, 396.0
    yRoot = 0.5 * (yA + yB)
    y3 = [296.0, 388.0, 480.0]

    node(ax, X0, yRoot, W0, 96, "人工建构的初始几何", "33 个", "见下方结构图谱",
         "#ffffff", INK, T, ts=19, cs=16)

    elbow(ax, X0 + W0, yRoot, X1 - 8, yA)
    node(ax, X1, yA, W1, 96, "理想几何 · ideal", "33 个几何", "直接算静态参考态，不弛豫",
         "#f3efe7", "#c9a961", T)

    elbow(ax, X0 + W0, yRoot, X1 - 8, yB)
    node(ax, X1, yB, W1, 96, r"在公共参考电势 $\mu_0$ 下弛豫", "16 个几何 · relax / relaxed",
         "给出每个结构的局部参考态", "#f3efe7", "#c9a961", T)

    kids = [("随机位移 · pert05 / pert10", "16 + 16 个几何", "可动原子随机位移 0.05 / 0.10 " + AA),
            ("集体变形 · coll", "14 个几何", "顶层间距 -3%、面内应变 +1%、台阶边缘弯曲"),
            ("路径构型 · path", "12 个几何", "原子沿指定路径移动（过桥位、脱离、进出）")]
    for (t_, c_, n_), y in zip(kids, y3):
        elbow(ax, X1 + W1, yB, X2 - 8, y)
        node(ax, X2, y, W2, 88, t_, c_, n_, "#eef2f5", "#7fa8bd", T, ts=17.5, cs=14.5, ns=13)

    fy = 594.0
    rbox(ax, 24, fy, W - 48, 74, "#f2f6f8", "#9dbccb", 1.6)
    tx(ax, 44, fy - 16, "再乘上电势", 16, INK, "semibold", store=T)
    tx(ax, 44, fy + 12, r"同一个几何在若干 TARGETMU 下各算一个电子态，U = $\mu_0-\mu_e$。"
                        r"±0.2 V 窗口 291 个电子态已完成；±0.5 V 扩展进行中。并非每个几何都有五个电势点。",
       14.5, SUB, store=T)

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
