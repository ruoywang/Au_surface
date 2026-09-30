#!/usr/bin/env python3
"""The structure-generation lineage: where every structure in the library came from, and how one built geometry
becomes several computed configurations.

This answers a question none of the 33 structure figures answer: which comparisons are CONTROLLED. From the
gallery alone, Island-7-compact and Island-7-elongated look like two unrelated islands; the lineage says they
hold the same seven atoms in a different shape in the same cell. Every arrow carries the geometric OPERATION,
not just the product's name, and every operation is taken from the build script that actually wrote the POSCAR.

Two things it deliberately makes explicit:
  * Au(211)/(221)/(332)/(554) are NOT Au(111) plus a strip. They are cut from the crystal at those Miller
    indices by different generators (ase.build.fcc211 for (211); a cubic cell plus surface() for the other
    three), so they sit on their own branch.
  * The 107 geometries are not 107 independently designed defects. 33 are built by hand; the other 74 are the
    relaxed, perturbed, collectively deformed and path images of those 33.

Drawn with matplotlib rather than emitted as SVG so the result can be read back and inspected, and so the text
extents are real: check_no_overlap() asserts that no two labels collide, which is the failure this figure is
most exposed to.

Usage (from Au_Cl/):  scripts/pyrun.sh scripts/render_lineage.py
Output: analysis/gallery/_lineage.png
"""
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
OUT = f"{ROOT}/analysis/gallery/_lineage.png"
CJK = matplotlib.font_manager.FontProperties(fname="/usr/share/fonts/google-droid/DroidSansFallback.ttf")

INK, SUB, FAINT, LINE = "#14181c", "#5a616a", "#9aa1a8", "#8b9299"
AA = r"$\mathrm{\AA}$"   # Droid Sans Fallback has no U+00C5, U+2212 or U+2080; those come from mathtext
BOX_A, EDGE_A = "#f3efe7", "#c9a961"          # the Au(111) branch
BOX_B, EDGE_B = "#e9f0f4", "#7fa8bd"          # the high-index branch

# (family, n_structures, operation lines, member lines, one-line caveat)
ROWS_A = [
    ("平整 Au(111)", 3,
     ["切出四层 (111) 平板（层高 5.0 / 7.4 / 9.8 / 12.2 " + AA + "）；", "只改面内周期胞，表面本身不动"],
     ["T-4x4（4×4）　Flat-8x2（8×2）　Flat-16x1（16×1）"],
     "共同基底与尺寸／镜像间距参考，不是三种缺陷"),
    ("点缺陷 / 极小团簇", 6,
     ["顶层移走 1 / 2 / 3 个 Au；", "或在三重空位上加 1 或 3 个 Au"],
     ["V1 · V2 · V3（移走的位点）", "A1-fcc · A1-hcp（同一个加原子，落在不同空位）　A3（三原子团簇）"],
     "V 是移走，A 是加上；A1-fcc 与 A1-hcp 只差落位"),
    ("条带直台阶", 5,
     ["按 fcc 延续注册，在顶层加一条有限宽的单层条带", "（条带恒占胞长的一半）"],
     ["Step-8x1 · Step-16x1 · Step-24x1（上下台面各 4 / 8 / 12 行）",
      "Step-8x2 · Step-16x2（沿台阶方向周期 ×2）"],
     "一次同时产生上台面、下台面与两条不等价边缘"),
    ("拐角 / 边缘重排", 3,
     ["从直台阶出发：在一条边多放一个原子；", "或把一个边缘原子移到脚部空位"],
     ["Kink-edge1 · Kink-edge2（8×3 条带，edge1 / edge2 各加一个）",
      "Step-8x2_edge-vacancy_plus_foot-adatom（总 Au 数不变）"],
     "后者与 Step-8x2 原子数、原子序都相同，只是一颗 Au 换了位置"),
    ("单层岛", 4,
     ["在顶层之上加 7 或 19 个 Au；", "或把同样 7 个 Au 改排成 3+4 两排"],
     ["Island-7-compact · Island-19-8x8（改尺寸）　Island-7-elongated（改形状，原子数不变）",
      "Island-7-8x8（同一缺陷，扩大周围周期胞）"],
     "尺寸、形状、镜像间距是三组分开的对照"),
    ("单层坑", 4,
     ["只从顶层移走 7 或 19 个 Au；", "或把 7 个缺失位改成 3+4 沟槽"],
     ["Pit-7-compact · Pit-19-8x8（改宽度）　Pit-7-trench（改边缘形状）",
      "Pit-7-8x8（同一缺陷，扩大周围周期胞）"],
     "坑深恒为一层，比较的是宽度与边缘形状"),
    ("重构相关堆垛", 2,
     ["把整个顶层平移到 hcp 注册；", "或在 16 个位点的长度内放 17 个顶层 Au"],
     ["R1-hcp-terminated（整层改注册）", "R2-stripe-wall（多一个原子，形成注册过渡带）"],
     "R1 不是一个 hcp 位吸附原子；R2 不是完整鱼骨重构"),
    ("复合形貌", 2,
     ["在台阶脚部加七原子岛；", "或把坑里移走的 Au 就地用来堆岛"],
     ["C1-island-near-step（Step-8x4 + 7 个 Au，实际与台阶相连）",
      "C2-island+pit（8×8 胞内移走 7 个、加回 7 个，总 Au 数不变）"],
     "C1 的岛与台阶连通，不是孤立的岛"),
]
ROWS_B = [
    ("台面宽度系列", 3,
     ["由立方胞直接按 (hkl) 切割；", "面间距 a/(2√(h²+k²+l²)) = 0.693 / 0.443 / 0.256 " + AA],
     ["Au221 · Au332 · Au554", "与 (111) 的名义夹角 15.79° → 10.02° → 5.77°（台面渐宽）"],
     "同一条构建路线上的宽度系列"),
    ("另一类台阶环境", 1,
     ["由 ase.build.fcc211 直接切割，与上面三者不是同一个生成器"],
     ["Au211（与 (111) 的名义夹角 19.47°）"],
     "不在上面的宽度系列里"),
]

ROWH, Y0 = 84.0, 120.0
X_ROOT, W_ROOT = 18.0, 156.0
X_ROUTE, W_ROUTE = 214.0, 200.0
X_OP = 448.0
X_FAM, W_FAM = 916.0, 200.0
X_MEM = 1136.0
W = 1720.0
H = Y0 + ROWH * (len(ROWS_A) + len(ROWS_B)) + 452


def rbox(ax, x, y, w, h, fc, ec, lw=1.4):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=7",
                                facecolor=fc, edgecolor=ec, linewidth=lw, zorder=2))


def tx(ax, x, y, s, size=12.5, color=INK, weight="normal", ha="left", va="center", store=None):
    t = ax.text(x, y, s, fontproperties=CJK, fontsize=size, color=color, fontweight=weight,
                ha=ha, va=va, zorder=5)
    if store is not None: store.append(t)
    return t


def elbow(ax, x0, y0, x1, y1, color=LINE, lw=1.3, arrow=False):
    xm = x0 + (x1 - x0) * 0.42
    ax.plot([x0, xm, xm, x1], [y0, y0, y1, y1], color=color, lw=lw, zorder=1,
            solid_capstyle="round", solid_joinstyle="round")
    if arrow:
        ax.annotate("", xy=(x1, y1), xytext=(x1 - 16, y1), zorder=3,
                    arrowprops=dict(arrowstyle="-|>", lw=lw, color=color, shrinkA=0, shrinkB=0))


def check_layout(fig, texts, tol=1.0):
    """Every label must have the figure to itself AND stay inside the canvas.

    The layout is written in POINTS, while a text's extent comes back in device pixels, so this is also the
    check that the two are consistent: get the figure scale wrong and every line of every block collides.
    """
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    bb = [(t, t.get_window_extent(renderer=r)) for t in texts]
    bad = []
    fw, fh = fig.canvas.get_width_height()
    for t, a in bb:
        if a.x0 < -tol or a.y0 < -tol or a.x1 > fw + tol or a.y1 > fh + tol:
            bad.append((f"OFF-CANVAS {t.get_text()[:40]!r}", "", 0, 0))
    for i in range(len(bb)):
        for j in range(i + 1, len(bb)):
            a, b = bb[i][1], bb[j][1]
            ov_x = min(a.x1, b.x1) - max(a.x0, b.x0)
            ov_y = min(a.y1, b.y1) - max(a.y0, b.y0)
            if ov_x > tol and ov_y > tol:
                bad.append((bb[i][0].get_text()[:34], bb[j][0].get_text()[:34],
                            round(float(ov_x), 1), round(float(ov_y), 1)))
    return bad


def main():
    # The layout above is in POINTS, so the figure is W/72 inches wide and a label's size in points is its
    # size in layout units. DPI only sets how many pixels that is worth.
    DPI = 130
    fig = plt.figure(figsize=(W / 72.0, H / 72.0), dpi=DPI)
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.axis("off")
    ax.add_patch(plt.Rectangle((0, 0), W, H, facecolor="white", edgecolor="none", zorder=0))
    T = []

    tx(ax, 20, 36, "从 Au 晶体到缺陷构型：结构库的生成与采样", 21, INK, "semibold", store=T)
    tx(ax, 20, 66, "每条箭头写的是几何操作，不是结构名称；操作取自实际写出 POSCAR 的建构脚本。"
                   "上半部分说明结构种类怎么来，下半部分说明每个结构怎么变成多个计算构型。", 12, SUB, store=T)

    rows = [(r, "A") for r in ROWS_A] + [(r, "B") for r in ROWS_B]
    ys = [Y0 + ROWH * i + ROWH / 2 for i in range(len(rows))]
    yA = float(np.mean(ys[:len(ROWS_A)])); yB = float(np.mean(ys[len(ROWS_A):]))
    ymid = (yA + yB) / 2

    rbox(ax, X_ROOT, ymid - 32, W_ROOT, 64, "#ffffff", INK, 1.8)
    tx(ax, X_ROOT + W_ROOT / 2, ymid - 10, "fcc Au 晶体", 15.5, INK, "semibold", "center", store=T)
    tx(ax, X_ROOT + W_ROOT / 2, ymid + 13, r"建构用 $a_0$ = 4.158 " + AA, 11, SUB, ha="center", store=T)

    for y, fc, ec, lines in ((yA, BOX_A, EDGE_A, ["路线一", "切出 Au(111) 四层平板", "再做局部增删 / 重排"]),
                             (yB, BOX_B, EDGE_B, ["路线二", "按高指数晶面直接切割", "得到规则台阶面"])):
        elbow(ax, X_ROOT + W_ROOT, ymid, X_ROUTE, y)
        rbox(ax, X_ROUTE, y - 44, W_ROUTE, 88, fc, ec, 1.8)
        tx(ax, X_ROUTE + W_ROUTE / 2, y - 24, lines[0], 10.5, SUB, "semibold", "center", store=T)
        tx(ax, X_ROUTE + W_ROUTE / 2, y, lines[1], 12.5, INK, "semibold", "center", store=T)
        tx(ax, X_ROUTE + W_ROUTE / 2, y + 24, lines[2], 11, SUB, ha="center", store=T)

    for i, ((fam, n, ops, mem, note), route) in enumerate(rows):
        y = ys[i]
        src_y, ec, fc = (yA, EDGE_A, BOX_A) if route == "A" else (yB, EDGE_B, BOX_B)
        elbow(ax, X_ROUTE + W_ROUTE, src_y, X_OP - 12, y, ec, 1.2)
        for k, ln in enumerate(ops):
            tx(ax, X_OP, y + 17 * k - 8.5 * (len(ops) - 1), ln, 12, INK, store=T)
        ax.annotate("", xy=(X_FAM - 6, y), xytext=(X_FAM - 30, y), zorder=3,
                    arrowprops=dict(arrowstyle="-|>", lw=1.6, color=ec, shrinkA=0, shrinkB=0))
        rbox(ax, X_FAM, y - 26, W_FAM, 52, fc, ec)
        tx(ax, X_FAM + W_FAM / 2, y - 7, fam, 13, INK, "semibold", "center", store=T)
        tx(ax, X_FAM + W_FAM / 2, y + 13, f"{n} 个结构", 10.5, SUB, ha="center", store=T)
        top = y - 8.5 * (len(mem) - 1) - 6
        for k, ln in enumerate(mem):
            tx(ax, X_MEM, top + 17 * k, ln, 11.5, SUB, store=T)
        tx(ax, X_MEM, top + 17 * (len(mem) - 1) + 20, "· " + note, 10.5, FAINT, store=T)

    # ---------- lower panel ----------
    yb = Y0 + ROWH * len(rows) + 46
    ax.plot([20, W - 20], [yb - 16, yb - 16], color="#e2ddd4", lw=1.5, zorder=1)
    tx(ax, 20, yb + 12, "一个建构好的几何，如何变成多个被计算的构型", 16.5, INK, "semibold", store=T)
    tx(ax, 20, yb + 36, "上面 33 个结构是人工建构的；余下的几何都是它们的弛豫、扰动与路径像。"
                        "合计 107 个几何，不是 107 个独立设计的缺陷。", 12, SUB, store=T)

    ty = yb + 74
    rbox(ax, 20, ty, 232, 58, "#ffffff", INK, 1.8)
    tx(ax, 136, ty + 22, "人工建构的初始几何", 12.5, INK, "semibold", "center", store=T)
    tx(ax, 136, ty + 42, "33 个（见上半部分）", 10.5, SUB, ha="center", store=T)

    for k, (name, cnt, note) in enumerate(
            [("理想几何 · ideal", "33", "直接作为静态参考态计算，不弛豫"),
             (r"在公共参考电势 $\mu_0$ 下弛豫 · relax / relaxed", "16", "给出每个结构的局部参考态")]):
        y = ty + 6 + 78 * k
        elbow(ax, 252, ty + 29, 344, y + 23, LINE, 1.2)
        rbox(ax, 344, y, 300, 46, "#f6f4f0", "#c9c3b8")
        tx(ax, 358, y + 17, name, 12, INK, "semibold", store=T)
        tx(ax, 358, y + 34, f"{cnt} 个几何　·　{note}", 10.3, SUB, store=T)

    y_par = ty + 6 + 78
    for k, (name, cnt, note) in enumerate(
            [("随机位移 · pert05 / pert10", "16 + 16", "可动原子随机位移 0.05 / 0.10 " + AA),
             ("集体变形 · coll", "14", "顶层间距 -3%、面内应变 +1%、台阶边缘弯曲 0.15 " + AA),
             ("路径构型 · path", "12", "吸附原子过桥位、边缘原子脱离到脚部、拐角原子沿边移动、岛／坑边原子进出")]):
        y = ty + 96 + 56 * k
        elbow(ax, 644, y_par + 23, 700, y + 21, LINE, 1.2)
        rbox(ax, 700, y, 596, 42, "#f6f4f0", "#c9c3b8")
        tx(ax, 714, y + 16, name, 12, INK, "semibold", store=T)
        tx(ax, 714, y + 32, f"{cnt} 个几何　·　{note}", 10.3, SUB, store=T)

    px, py = 1330, ty + 96
    rbox(ax, px, py, 372, 154, "#f2f6f8", "#9dbccb")
    tx(ax, px + 16, py + 24, "再乘上电势", 12.5, INK, "semibold", store=T)
    tx(ax, px + 16, py + 47, "同一个几何在若干 TARGETMU 下各算一个电子态；", 11, SUB, store=T)
    tx(ax, px + 16, py + 65, r"U = $\mu_0-\mu_e$，$\mu_0$ = $-$4.9071 eV。", 11, SUB, store=T)
    tx(ax, px + 16, py + 91, "±0.2 V 窗口：291 个电子态，已完成。", 11, INK, store=T)
    tx(ax, px + 16, py + 111, "±0.5 V 扩展：进行中，本页结论不使用。", 11, INK, store=T)
    tx(ax, px + 16, py + 134, "并非每个几何都有五个电势点。", 10.3, FAINT, store=T)

    cy = ty + 272
    rbox(ax, 20, cy, 1276, 64, "#fdf6e8", "#e0c58a")
    tx(ax, 38, cy + 23, "这套采样是什么，不是什么", 12.5, INK, "semibold", store=T)
    tx(ax, 38, cy + 45, "结构类别由晶体学与几何操作构建；弛豫提供局部参考态；扰动与路径用于覆盖非平衡构型。"
                        "它们不是分子动力学中自然出现频率的统计，也不是 DFT 自动找出的全部稳定形貌。", 11.3, SUB,
       store=T)

    bad = check_layout(fig, T)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, facecolor="white")
    plt.close(fig)
    if bad:
        print(f"LAYOUT PROBLEMS ({len(bad)}):")
        for a, b, ox, oy in bad[:12]: print(f"   {a!r} x {b!r}  ({ox} x {oy} px)")
        raise SystemExit(1)
    print(f"{os.path.getsize(OUT)/1024:.0f} kB -> {OUT}  ({int(W/72*DPI)}x{int(H/72*DPI)} px, "
          f"{len(T)} labels, none overlapping, none off-canvas)")


if __name__ == "__main__":
    main()
