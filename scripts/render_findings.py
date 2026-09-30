#!/usr/bin/env python3
"""The two lead figures: what the charging data says that a heat map does not.

_finding_local.png -- the electrode can charge more while the liquid above the protruding site enriches LESS.
    Left: A1-hcp's region-averaged anion enrichment at U = +0.20 V, adatom region against terrace, beside the
    whole-cell charge of A1-hcp and of the flat member of the same cell.
    Right: the metal's response and the anion's response across a strip step, both normalised to their own
    cell mean, so the shapes can be compared without comparing units.

_finding_move.png -- moving ONE Au, with the atom count unchanged, shifts the charging curve sideways while
    leaving its slope nearly alone. Step-8x2 against Step-8x2_edge-vacancy_plus_foot-adatom.

Everything is read from the frozen +-0.2 V production analysis; nothing here uses the +-0.5 V extension.

Usage (from Au_Cl/):  scripts/pyrun.sh scripts/render_findings.py
Outputs: analysis/gallery/_finding_local.png, analysis/gallery/_finding_move.png
"""
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
from render_gallery import (flatten_cell, top_atoms, coordination, cn_class, CN_COLOR, draw,  # noqa: E402
                            parent_of, R_AU, CJK, AA)

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
GAL = f"{ROOT}/analysis/gallery"
PC = f"{ROOT}/analysis/spatial/percolumn"
INK, SUB, FAINT = "#14181c", "#5a616a", "#9aa1a8"
C_METAL, C_ION = "#b0543a", "#2e7d9a"


def charging():
    return json.load(open(f"{ROOT}/analysis/charging/charging.json"))["geometries"]


def geom(G, sid, cfg="ideal"):
    for g in G.values():
        if g["structure_id"] == sid and g["config"] == cfg: return g
    raise KeyError(f"{sid}/{cfg}")


def regions():
    return json.load(open(f"{ROOT}/analysis/spatial/regions.json"))["states"]


def smooth(y, L, w=1.0):
    n = len(y); k = max(1, int(round(w / (L / n))))
    ker = np.ones(2 * k + 1) / (2 * k + 1)
    return np.convolve(np.r_[y[-k:], y, y[:k]], ker, mode="valid")[:n]


def step_profiles(sid):
    """(s, metal response, anion response, upper-terrace mask, L), each normalised to its own cell mean.

    Per-column arrays are stored (ngy, ngx): averaging over y runs along the step edge, leaving a profile
    across the terrace. Both responses are the change from U = -0.2 V to U = +0.2 V."""
    lo = np.load(f"{PC}/{sid}__ideal__mu-4.7071.npz"); hi = np.load(f"{PC}/{sid}__ideal__mu-5.1071.npz")
    L = float(np.linalg.norm(hi["cell"][0][:2]))
    dq = -(hi["ne_col"].mean(0) - lo["ne_col"].mean(0))
    dg = hi["gamma_rel"].mean(0) - lo["gamma_rel"].mean(0)
    h = smooth(hi["top_height"].mean(0), L, 1.5)
    n = len(dq); s = (np.arange(n) + 0.5) / n * L
    return s, smooth(dq, L) / dq.mean(), smooth(dg, L) / dg.mean(), h > 0.5 * (h.max() + h.min()), L


def check(fig, texts, tol=1.0):
    """No label may overlap another or run off the canvas. The second half matters as much as the first: a
    caption that overflows the right edge is silently truncated in the PNG."""
    fig.canvas.draw(); r = fig.canvas.get_renderer()
    bb = [(t, t.get_window_extent(renderer=r)) for t in texts]
    bad = []
    fw, fh = fig.canvas.get_width_height()
    for t, a in bb:
        if a.x0 < -tol or a.y0 < -tol or a.x1 > fw + tol or a.y1 > fh + tol:
            bad.append((f"OFF-CANVAS {t.get_text()[:44]!r}", ""))
    for i in range(len(bb)):
        for j in range(i + 1, len(bb)):
            a, b = bb[i][1], bb[j][1]
            if min(a.x1, b.x1) - max(a.x0, b.x0) > tol and min(a.y1, b.y1) - max(a.y0, b.y0) > tol:
                bad.append((bb[i][0].get_text()[:36], bb[j][0].get_text()[:36]))
    return bad


# --------------------------------------------------------------------------------------- figure 1
def fig_local():
    G, R = charging(), regions()
    st = R["A1-hcp__ideal__mu-5.1071"]; U = st["U"]
    reg = st["regions"]
    rows = [("吸附原子所在的低配位区", reg["kink/adatom"]["K_rel"], reg["kink/adatom"]["area_fraction"]),
            ("周围台面区", reg["terrace"]["K_rel"], reg["terrace"]["area_fraction"]),
            ("次表面 / 脚部区", reg["sub-surface/foot"]["K_rel"], reg["sub-surface/foot"]["area_fraction"]),
            ("整胞平均", reg["all"]["K_rel"], 1.0)]
    g_ad = geom(G, "A1-hcp"); g_fl = geom(G, "T-4x4")
    p_ad = max(g_ad["points"], key=lambda p: p["U"]); p_fl = max(g_fl["points"], key=lambda p: p["U"])

    fig = plt.figure(figsize=(14.2, 7.0), dpi=170)
    # the left column is split so the sigma comparison has its own space instead of hanging off the axes
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.22],
                          left=0.235, right=0.985, top=0.80, bottom=0.235, wspace=0.14)
    T = []
    T.append(fig.text(0.020, 0.955, "整体充多少电，和液相在哪里富集，不是同一件事",
                      fontproperties=CJK, fontsize=20.6, color=INK, va="top", weight="semibold"))
    T.append(fig.text(0.020, 0.895, "左：同一个电势下，A1-hcp 的吸附原子区平均阴离子富集低于它周围的台面。"
                                    "右：条带台阶上，金属响应的峰在台阶边，阴离子响应的峰不在。",
                      fontproperties=CJK, fontsize=13.8, color=SUB, va="top"))

    ax = fig.add_subplot(gs[0, 0])
    y = np.arange(len(rows))[::-1]
    vals = [(k - 1) * 100 for _, k, _ in rows]
    cols = ["#c0392b", "#d8b34a", "#9aa3ad", "#5a616a"]
    ax.barh(y, vals, height=0.56, color=cols, edgecolor="#2b2f36", linewidth=0.6)
    for yy, (name, k, af), v in zip(y, rows, vals):
        T.append(ax.text(v + 0.35, yy, f"K = {k:.4f}", fontproperties=CJK, fontsize=13.1, va="center", color=INK))
        T.append(ax.text(-0.02, (yy + 0.75) / len(rows), f"{name}（占面积 {100*af:.0f}%）",
                         transform=ax.transAxes, fontproperties=CJK, fontsize=13.1,
                         va="center", ha="right", color=INK))
    ax.axvline(0, color="#2b3137", lw=1.2)
    ax.set_xlim(-9.0, 26.0); ax.set_ylim(-0.75, len(rows) - 0.25)
    ax.set_yticks([]); ax.tick_params(axis="x", labelsize=9.5, colors=SUB)
    ax.set_xlabel(r"区域平均阴离子富集相对体相的偏离 ($K-1$)，%", fontproperties=CJK, fontsize=13.1, color=SUB)
    for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_color("#c6c1b7")
    T.append(ax.set_title(f"A1-hcp（一个吸附原子），U = {U:+.4f} V", fontproperties=CJK, fontsize=15.0,
                          color=INK, pad=8))
    # a full-width caption rather than a box inside the left column, which the longest line overflowed
    box = (f"同一电势附近的整胞面电荷：A1-hcp σ = {p_ad['sigma_uC_per_cm2']:+.2f} μC/cm²"
           f"（U = {p_ad['U']:+.4f} V）；T-4x4 σ = {p_fl['sigma_uC_per_cm2']:+.2f} μC/cm²"
           f"（U = {p_fl['U']:+.4f} V）。\n"
           "加上吸附原子后整个电极更正，但它上方区域的平均富集反而低于台面。两者电势接近但不完全相等，只作量级比较。\n"
           "K > 1 说明吸附原子区仍然富集，只是不如台面；这不是“吸附原子排斥阴离子”。")
    T.append(fig.text(0.020, 0.135, box, fontproperties=CJK, fontsize=12.5, color=INK, va="top", ha="left",
                      linespacing=1.75))

    ax2 = fig.add_subplot(gs[0, 1])
    sid = "Step-16x1"
    s, dq, dg, up, L = step_profiles(sid)
    # draw the lower terrace shaded so the offset direction is readable
    ax2.fill_between(s, 0, 1, where=~up, transform=ax2.get_xaxis_transform(), color="#f0ece4",
                     lw=0, zorder=0)
    ax2.plot(s, dq, color=C_METAL, lw=2.2, zorder=3)
    ax2.plot(s, dg, color=C_ION, lw=2.2, zorder=3)
    ax2.axhline(1.0, color="#b9b3a7", lw=0.9, ls=(0, (4, 3)), zorder=1)
    n = len(s)
    bnd = [i for i in range(n) if up[i] != up[(i - 1) % n]]
    offs = []
    for i in bnd:
        into = -1 if up[i] else +1
        w = int(round(7.0 / (L / n)))
        jm = max([(i + d) % n for d in range(-3, 4)], key=lambda j: dq[j])
        jg = max([(i + into * k) % n for k in range(-2, w)], key=lambda j: dg[j])
        offs.append(((jg - jm + n // 2) % n - n // 2) * L / n * into)
        ax2.axvline(s[i], color="#8b9299", lw=1.0, ls=(0, (5, 3)), zorder=1)
        ax2.plot([s[jm]], [dq[jm]], "o", ms=6, mfc=C_METAL, mec="white", mew=1.2, zorder=4)
        ax2.plot([s[jg]], [dg[jg]], "o", ms=6, mfc=C_ION, mec="white", mew=1.2, zorder=4)
        ax2.annotate("", xy=(s[jg], dg[jg] + 0.06), xytext=(s[jm], dg[jg] + 0.06), zorder=5,
                     arrowprops=dict(arrowstyle="<|-|>", lw=1.1, color="#5a616a", shrinkA=0, shrinkB=0))
        T.append(ax2.text(0.5 * (s[jm] + s[jg]), dg[jg] + 0.10, f"{abs(offs[-1]):.1f} " + AA,
                          fontproperties=CJK, fontsize=11.9, color="#5a616a", ha="center", va="bottom"))
    ax2.set_xlim(0, L); ax2.set_ylim(0.55, max(dq.max(), dg.max()) + 0.42)
    ax2.set_xlabel("沿台阶法向的位置 (" + AA + ")", fontproperties=CJK, fontsize=13.1, color=SUB)
    ax2.set_ylabel("响应 / 该量自身的整胞平均", fontproperties=CJK, fontsize=13.1, color=SUB)
    ax2.tick_params(labelsize=9.5, colors=SUB)
    for sp in ("top", "right"): ax2.spines[sp].set_visible(False)
    for sp in ("bottom", "left"): ax2.spines[sp].set_color("#c6c1b7")
    T.append(ax2.set_title(f"{sid}：U 从 $-$0.2 V 到 +0.2 V 的响应（阴影为下台面）",
                           fontproperties=CJK, fontsize=15.0, color=INK, pad=8))
    T.append(ax2.text(0.015, 0.965, r"金属：正电荷增量 ($-\Delta n_e$)", transform=ax2.transAxes, fontproperties=CJK,
                      fontsize=13.1, color=C_METAL, va="top", weight="semibold"))
    T.append(ax2.text(0.015, 0.895, r"液相：阴离子过量增量 ($\Delta\Gamma_-$)", transform=ax2.transAxes, fontproperties=CJK,
                      fontsize=13.1, color=C_ION, va="top", weight="semibold"))
    T.append(ax2.text(0.985, 0.045, f"峰值：金属 {dq.max():.2f}×，阴离子 {dg.max():.2f}× 各自的整胞平均",
                      transform=ax2.transAxes, fontproperties=CJK, fontsize=12.2, color=INK,
                      ha="right", va="bottom"))

    bad = check(fig, T)
    fig.savefig(f"{GAL}/_finding_local.png", facecolor="white"); plt.close(fig)
    return bad, dict(U=U, K_adatom=reg["kink/adatom"]["K_rel"], K_terrace=reg["terrace"]["K_rel"],
                     K_all=reg["all"]["K_rel"], area_adatom=reg["kink/adatom"]["area_fraction"],
                     sigma_A1hcp=p_ad["sigma_uC_per_cm2"], U_A1hcp=p_ad["U"],
                     sigma_T4x4=p_fl["sigma_uC_per_cm2"], U_T4x4=p_fl["U"],
                     step=sid, metal_peak=float(dq.max()), ion_peak=float(dg.max()),
                     offsets_A=[float(o) for o in offs])


# --------------------------------------------------------------------------------------- figure 2
def fig_move():
    G = charging()
    A, B = "Step-8x2", "Step-8x2_edge-vacancy_plus_foot-adatom"
    gA, gB = geom(G, A), geom(G, B)
    fig = plt.figure(figsize=(13.6, 6.2), dpi=170)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.45], left=0.02, right=0.975, top=0.79,
                          bottom=0.145, wspace=0.14)
    T = []
    T.append(fig.text(0.02, 0.955, "原子数不变，只把一颗 Au 换个位置，充电曲线就整体平移",
                      fontproperties=CJK, fontsize=20.6, color=INK, va="top", weight="semibold"))
    T.append(fig.text(0.02, 0.895, "同一个胞、同样 72 个 Au：把台阶上边缘的一颗原子移到脚部空位。"
                                   "零电荷电势移动约 55 mV，而割线电容几乎不变。",
                      fontproperties=CJK, fontsize=13.8, color=SUB, va="top"))

    gal = json.load(open(f"{GAL}/gallery.json"))
    # one cell only: with the defect repeated over a tiling, six green circles read as six moved atoms
    atB = flatten_cell(read(f"{gal[B]['dir']}/{gal[B]['geom']}"))
    pa = parent_of(B, atB)                       # added = where the atom went, removed = where it came from
    for k, sid in enumerate([A, B]):
        m = gal[sid]
        at = flatten_cell(read(f"{m['dir']}/{m['geom']}")); cell = at.get_cell().array
        pos = at.get_positions(); cn = coordination(at)
        col = np.array([CN_COLOR[cn_class(c)] for c in cn]); exp = top_atoms(at)
        keep = pos[:, 2] > pos[:, 2].max() - 5.2
        ax = fig.add_subplot(gs[0, k]); ax.set_aspect("equal"); ax.axis("off")
        P, C, Z, E = [], [], [], []
        for j in range(2):                        # repeat only along the short vector, to show the periodicity
            sh = (j * cell[1])[:2]
            P.append(pos[keep][:, :2] + sh); C.append(col[keep]); Z.append(pos[keep][:, 2]); E.append(exp[keep])
        P = np.vstack(P); C = np.concatenate(C); Z = np.concatenate(Z); E = np.concatenate(E)
        draw(ax, P, C, Z, E, R_AU)
        if pa is not None:
            for j in range(2):
                sh = (j * cell[1])[:2]
                if k == 0:                        # the parent: mark the atom that is about to move
                    for c in pa["removed"]:
                        ax.add_patch(Circle(c + sh, R_AU * 1.5, facecolor="none", edgecolor="#8a3ffc",
                                            lw=2.2, zorder=6000))
                else:                             # the child: where it came from, and where it went
                    for c in pa["removed"]:
                        ax.add_patch(Circle(c + sh, R_AU * 1.2, facecolor="none", edgecolor="#8a3ffc",
                                            lw=2.0, ls=(0, (3, 2)), zorder=6000))
                    for c in pa["added"]:
                        ax.add_patch(Circle(c + sh, R_AU * 1.5, facecolor="none", edgecolor="#1f6f3f",
                                            lw=2.2, zorder=6000))
                    for c0, c1 in zip(pa["removed"], pa["added"]):
                        ax.annotate("", xy=tuple(c1 + sh), xytext=tuple(c0 + sh), zorder=6200,
                                    arrowprops=dict(arrowstyle="-|>", lw=1.8, color="#1f6f3f",
                                                    shrinkA=9, shrinkB=9))
        ax.set_xlim(P[:, 0].min() - 2, P[:, 0].max() + 2); ax.set_ylim(P[:, 1].min() - 2, P[:, 1].max() + 2)
        lab = "直台阶 Step-8x2（紫圈＝将要移动的那颗）" if k == 0 else "移位后（箭头＝它去了哪里）"
        T.append(ax.set_title(lab, fontproperties=CJK, fontsize=14.4, color=INK, pad=6))
    T.append(fig.text(0.02, 0.055, "紫色圈＝原来的边缘位点　绿色圈＝它的新位置（台阶脚部空位）　"
                                   "两图各画 1 个胞、沿短边重复 2 次",
                      fontproperties=CJK, fontsize=12.0, color=SUB, va="bottom"))

    ax = fig.add_subplot(gs[0, 2])
    out = {}
    for sid, g, c, mk in ((A, gA, "#2b6cb0", "o"), (B, gB, "#c0392b", "s")):
        pts = sorted(g["points"], key=lambda p: p["U"])
        U = np.array([p["U"] for p in pts]); S = np.array([p["sigma_uC_per_cm2"] for p in pts])
        Cm = float(np.median([x["C_uF_per_cm2"] for x in g["secants"]]))
        pz = g["pzc"]["U_pzc"]
        ax.plot(U, S, marker=mk, ms=7, lw=2.0, color=c, mec="white", mew=1.2, zorder=3,
                label=f"{sid}")
        ax.plot([pz], [0], marker="v", ms=9, color=c, mec="white", mew=1.0, zorder=4, clip_on=False)
        out[sid] = dict(U_pzc_mV=1000 * pz, C=Cm)
    ax.axhline(0, color="#2b3137", lw=1.1, zorder=1)
    ax.axvline(0, color="#d8d3c9", lw=0.9, ls=(0, (4, 3)), zorder=1)
    d = out[B]["U_pzc_mV"] - out[A]["U_pzc_mV"]
    ax.annotate("", xy=(out[B]["U_pzc_mV"] / 1000, 0.55), xytext=(out[A]["U_pzc_mV"] / 1000, 0.55), zorder=5,
                arrowprops=dict(arrowstyle="<|-|>", lw=1.4, color=INK, shrinkA=0, shrinkB=0))
    T.append(ax.text(0.5 * (out[A]["U_pzc_mV"] + out[B]["U_pzc_mV"]) / 1000, 0.68,
                     f"零电荷电势移动 {d:+.1f} mV", fontproperties=CJK, fontsize=13.8, color=INK,
                     ha="center", va="bottom", weight="semibold"))
    T.append(ax.text(0.03, 0.965,
                     f"{A}：U0 = {out[A]['U_pzc_mV']:+.1f} mV，割线电容 {out[A]['C']:.2f} μF/cm²\n"
                     f"移位后：U0 = {out[B]['U_pzc_mV']:+.1f} mV，割线电容 {out[B]['C']:.2f} μF/cm²\n"
                     f"电容差 {100*(out[B]['C']/out[A]['C']-1):+.1f}%：两条曲线斜率接近，横向位置不同",
                     transform=ax.transAxes, fontproperties=CJK, fontsize=12.2, color=INK, va="top",
                     linespacing=1.6,
                     bbox=dict(boxstyle="round,pad=0.5", facecolor="#f6f4f0", edgecolor="#c9c3b8", lw=1.0)))
    ax.set_xlabel(r"U = $\mu_0-\mu_e$  (V，内部参考，不对 RHE)", fontproperties=CJK, fontsize=13.1, color=SUB)
    ax.set_ylabel("面电荷 σ (μC/cm²)", fontproperties=CJK, fontsize=13.1, color=SUB)
    ax.tick_params(labelsize=9.5, colors=SUB)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("bottom", "left"): ax.spines[sp].set_color("#c6c1b7")
    lg = ax.legend(loc="lower right", frameon=False, fontsize=11.9, prop=CJK)
    for t in lg.get_texts(): t.set_fontproperties(CJK); t.set_fontsize(9.5)
    T.append(ax.text(0.97, 0.30, "▼ 标记＝线性插值得到的零电荷电势", transform=ax.transAxes,
                     fontproperties=CJK, fontsize=11.6, color=FAINT, ha="right", va="bottom"))
    T.append(ax.set_title("两条 σ(U)：三个基础电势点", fontproperties=CJK, fontsize=15.0, color=INK, pad=6))

    bad = check(fig, T)
    fig.savefig(f"{GAL}/_finding_move.png", facecolor="white"); plt.close(fig)
    out["dU_pzc_mV"] = d
    out["dC_pct"] = 100 * (out[B]["C"] / out[A]["C"] - 1)
    return bad, out


def relax_shift():
    """How far relaxing the SAME structure moves its zero-charge potential, over every structure that has both.

    The point is not any single row: it is that the structure's NAME does not fix its charging response, so a
    ranking built on one idealised POSCAR per defect would be reading the idealisation as much as the defect."""
    G = charging()
    by = {}
    for g in G.values():
        pz = g["pzc"]
        if g["config"] in ("ideal", "relaxed") and pz.get("U_pzc") is not None and pz.get("bracketed"):
            by.setdefault(g["structure_id"], {})[g["config"]] = (
                1000 * pz["U_pzc"], float(np.median([s["C_uF_per_cm2"] for s in g["secants"]])))
    rows = [(k, v["ideal"], v["relaxed"]) for k, v in by.items() if "ideal" in v and "relaxed" in v]
    d = sorted([(k, r[0] - i[0], 100 * (r[1] / i[1] - 1), i[0], r[0]) for k, i, r in rows],
               key=lambda x: -abs(x[1]))
    return dict(n=len(d), median_abs_dU_mV=float(np.median([abs(x[1]) for x in d])),
                max_abs_dU_mV=float(max(abs(x[1]) for x in d)),
                median_abs_dC_pct=float(np.median([abs(x[2]) for x in d])),
                top=[dict(structure=k, dU_mV=du, dC_pct=dc, ideal_mV=i, relaxed_mV=r)
                     for k, du, dc, i, r in d[:4]])


def main():
    os.makedirs(GAL, exist_ok=True)
    b1, s1 = fig_local()
    b2, s2 = fig_move()
    s3 = relax_shift()
    json.dump(dict(local=s1, move=s2, relax=s3), open(f"{GAL}/_findings.json", "w"),
              indent=1, ensure_ascii=False)
    for name, bad in (("_finding_local", b1), ("_finding_move", b2)):
        if bad:
            print(f"OVERLAPPING LABELS in {name}:")
            for a, b in bad[:10]: print(f"   {a!r} x {b!r}")
    if b1 or b2: raise SystemExit(1)
    print(f"2 figures -> {GAL}/_finding_local.png, {GAL}/_finding_move.png  (no label overlaps)")
    print(json.dumps(dict(local=s1, move=s2, relax=s3), ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
