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


def cn_class(c):
    return "kink" if c <= 6 else "edge" if c <= 8 else "terrace" if c == 9 else "bulk"


def draw(ax, pts, colors, zs, r, lw=0.5):
    order = np.argsort(zs)
    zmin, zmax = (zs.min(), zs.max()) if len(zs) else (0, 1)
    for i in order:
        f = 0.45 + 0.55 * ((zs[i] - zmin) / (zmax - zmin) if zmax > zmin else 1.0)   # depth cue
        c = np.array(matplotlib.colors.to_rgb(colors[i]))
        ax.add_patch(Circle(pts[i], r, facecolor=tuple(c * f + (1 - f) * 0.92), edgecolor="#2b2f36",
                            linewidth=lw, zorder=int(1000 * f)))


TOP_DEPTH = 5.2         # A below the highest atom kept in the top view: two (111) terrace levels, no deep bulk
TILE_TARGET = 26.0      # A, the in-plane extent each tiled view aims for


def render(sid, path, title, meta):
    at = read(path)
    cell = at.get_cell().array
    pos = at.get_positions()
    cn = coordination(at)
    col = np.array([CN_COLOR[cn_class(c)] for c in cn])
    zt = pos[:, 2].max()

    n1 = max(1, min(4, int(round(TILE_TARGET / np.linalg.norm(cell[0][:2])))))
    n2 = max(1, min(4, int(round(TILE_TARGET / np.linalg.norm(cell[1][:2])))))
    keep = pos[:, 2] > zt - TOP_DEPTH

    P, C, Z = [], [], []
    for i in range(n1):
        for j in range(n2):
            sh = (i * cell[0] + j * cell[1])[:2]
            P.append(pos[keep][:, :2] + sh); C.append(col[keep]); Z.append(pos[keep][:, 2])
    P = np.vstack(P); C = np.concatenate(C); Z = np.concatenate(Z)

    ns = max(1, min(3, int(round(TILE_TARGET / np.linalg.norm(cell[0][:2])))))
    P2, C2, Z2 = [], [], []
    for i in range(ns):
        P2.append(np.c_[pos[:, 0] + i * cell[0][0], pos[:, 2]]); C2.append(col); Z2.append(pos[:, 1])
    P2 = np.vstack(P2); C2 = np.concatenate(C2); Z2 = np.concatenate(Z2)

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
    draw(ax, P, C, Z, R_AU)
    for i in range(n1):
        for j in range(n2):
            o = i * cell[0][:2] + j * cell[1][:2]
            ax.plot(*zip(o, o + cell[0][:2], o + cell[0][:2] + cell[1][:2], o + cell[1][:2], o),
                    color="#8b9299", lw=0.7, ls=(0, (4, 3)), zorder=5000)
    ax.set_xlim(P[:, 0].min() - 1.5 * R_AU, P[:, 0].max() + 1.5 * R_AU)
    ax.set_ylim(P[:, 1].min() - 1.5 * R_AU, P[:, 1].max() + 1.5 * R_AU)
    ax.set_title(f"top view · {n1}x{n2} cells · outermost {TOP_DEPTH:.1f} $\\AA$",
                 fontsize=8.5, color="#4a5158", pad=3)

    ax2 = fig.add_subplot(gs[0, 1]); ax2.set_aspect("equal"); ax2.axis("off")
    draw(ax2, P2, C2, Z2, R_AU)
    ax2.axhline(zt + 4.2, color="#2e7d9a", lw=1.1, ls=(0, (5, 3)), zorder=6000)
    ax2.text(P2[:, 0].min() - R_AU, zt + 4.8, "ion-accessible electrolyte above this line",
             fontsize=7.2, color="#2e7d9a")
    ax2.set_xlim(P2[:, 0].min() - 1.5 * R_AU, P2[:, 0].max() + 1.5 * R_AU)
    ax2.set_ylim(P2[:, 1].min() - 1.5 * R_AU, zt + 7.0)
    ax2.set_title(f"side view · {ns} cell(s) along a$_1$ · all layers", fontsize=8.5, color="#4a5158", pad=3)

    fig.text(0.5, 1 - 0.30 / FIGH, title, ha="center", va="top", fontsize=13.5, color="#14181c", weight="medium")
    sub = (f"{meta['n_atoms']} Au   ·   $A_{{proj}}$ = {meta['A_proj']:.0f} $\\AA^2$   ·   {meta['family']}"
           f"   ·   {meta['cn_counts']}")
    fig.text(0.5, 1 - 0.68 / FIGH, sub, ha="center", va="top", fontsize=8.6, color="#5a616a")
    h = [plt.Line2D([], [], marker="o", ls="", ms=7, mfc=CN_COLOR[k], mec="#2b2f36", mew=0.5, label=l)
         for k, l in [("kink", "CN$\\leq$6  kink / adatom"), ("edge", "CN 7-8  step edge, island & pit rim"),
                      ("terrace", "CN 9  terrace"), ("bulk", "CN$\\geq$10  sub-surface / step foot")]]
    fig.legend(handles=h, loc="lower center", ncol=4, fontsize=7.6, frameon=False,
               bbox_to_anchor=(0.5, 0.02 / FIGH))
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
