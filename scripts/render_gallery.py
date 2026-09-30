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


def render(sid, path, title, meta):
    at = read(path)
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
    gs = fig.add_gridspec(1, 2, width_ratios=[w1, w2], wspace=1.2 / (w1 + w2) * 2,
                          left=0.015, right=0.985, top=top_frac, bottom=bot_frac)

    ax = fig.add_subplot(gs[0, 0]); ax.set_aspect("equal"); ax.axis("off")
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

    ax2 = fig.add_subplot(gs[0, 1]); ax2.set_aspect("equal"); ax2.axis("off")
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
