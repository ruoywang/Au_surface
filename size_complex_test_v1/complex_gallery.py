#!/usr/bin/env python3
"""Gallery of the multi-layer parents and their representative cuts: for every complex candidate that passed the 6x6
(or 8x8) screen, one row = parent MD frame with the cut frame and the protected set -> zoom on the cut with the
protected atoms ringed (the declared target: upper edge + foot, both steps + terrace, both walls + floor) and the
repaired seam band -> the reconstructed cell tiled 2 x 2 with a side view. Same drawing conventions as the lineage
figure (discs by layer, nothing drawn that is not in the coordinates); layer colours extended to levels 5 and 6.

Usage (from Au_Cl/):  rough_sampling_v1/pyrun_rs.sh size_complex_test_v1/complex_gallery.py
Output: size_complex_test_v1/complex_gallery/<candidate>.png and complex_gallery/index.html, plus _parents_after_md.png
"""
import csv
import html
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from ase.io import read  # noqa: E402
from matplotlib.patches import Circle, Polygon  # noqa: E402

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
R = f"{ROOT}/rough_sampling_v1"; T = f"{ROOT}/size_complex_test_v1"
sys.path.insert(0, R)
import render_rough_lineage as RL  # noqa: E402
import trajio  # noqa: E402
from extract_cells import lattice, boundary_offset, layers_from_z, R_CORE  # noqa: E402

RL.LAYER_COL.update({5: "#8b1a1a", 6: "#4b0082", 7: "#2b0045"})
OUT = f"{T}/complex_gallery"


def side_view(ax, at, title):
    c = at.get_cell().array; P = at.get_positions(); lay = layers_from_z(at); z0 = P[:, 2].min()
    e1 = c[0][:2] / np.linalg.norm(c[0][:2]); u = P[:, :2] @ e1; h = P[:, 2] - z0
    for i in np.argsort(-(P[:, :2] @ (c[1][:2] / np.linalg.norm(c[1][:2])))):
        ax.add_patch(Circle((u[i], h[i]), 1.0, fc=RL.LAYER_COL.get(int(lay[i]), "#000"), ec="k", lw=0.3))
    ax.set_xlim(u.min() - 2, u.max() + 2); ax.set_ylim(-1.5, h.max() + 1.8); ax.set_aspect("equal"); ax.set_yticks([]); ax.set_xticks([])
    ax.set_title(title, fontsize=RL.FS)


def main():
    os.makedirs(OUT, exist_ok=True); trajio.MD_DIR = f"{T}/md"
    rows = [r for r in csv.DictReader(open(f"{T}/size_screen.complex.csv")) if r["preferred_cell"] in ("6x6", "8x8")]
    cands = {json.loads(l)["candidate_id"]: json.loads(l) for l in open(f"{T}/complex_centres.jsonl")}
    cards = []
    for r in rows:
        c = cands[r["state_id"]]; pid = c["parent_id"]; step = c["step"]
        parent = trajio.frame_by_step(pid, step); p0 = read(f"{T}/parents/{pid}.extxyz"); parent.set_array("fixed", p0.get_array("fixed"))
        cell_dir = r["new_cell_dir"]; cell_at = read(f"{cell_dir}/cell.extxyz"); info = json.load(open(f"{cell_dir}/cut_info.json"))
        a1, a2, cellp = lattice(parent); off = boundary_offset(parent, a1, a2); n1, n2 = (int(x) for x in r["preferred_cell"].split("x")); oij = info["origin"]
        cut_origin = (oij[0] + off[0]) * a1 + (oij[1] + off[1]) * a2
        Cp = np.array([cellp[0][:2], cellp[1][:2]]); f = cut_origin @ np.linalg.inv(Cp); f -= np.floor(f); cut_origin = f @ Cp
        v1, v2 = n1 * a1, n2 * a2; centre = c["atom"]; prot = np.array(c["protected_atom_ids"], int)
        fig = plt.figure(figsize=(27, 9.2)); axes = [fig.add_axes([0.03, 0.08, 0.30, 0.84]), fig.add_axes([0.36, 0.08, 0.30, 0.84]), fig.add_axes([0.69, 0.36, 0.29, 0.56])]
        sx = fig.add_axes([0.69, 0.06, 0.29, 0.24])
        c_img = RL.centre_image(parent, cut_origin, v1, v2, centre)
        RL.panel_parent(axes[0], parent, cut_origin, v1, v2, c_img, f"{pid} ({c['cls']}), {c['time_ps']:.0f} ps ({len(parent)} Au)\nblue = cut frame {r['preferred_cell']}, green = centre 6 Å")
        RL.panel_cut(axes[1], parent, cut_origin, v1, v2, centre, prot, f"target '{c['target_environment']}': {len(prot)} protected atoms (ringed) kept to 1e-3 Å\n{info['repair'][:70]}")
        RL.panel_cell(axes[2], cell_at, info["centre_in_cell"], f"reconstructed {r['preferred_cell']} ({len(cell_at)} Au), tiled 2 x 2; seam-affected {info['stats']['n_seam_affected']}")
        side_view(sx, cell_at, "side view along a2")
        fig.legend(handles=[plt.Line2D([], [], marker="o", ls="", ms=12, mfc=RL.LAYER_COL[k], mec="k", label=f"layer {k}") for k in range(2, 7)] +
                   [plt.Line2D([], [], marker="o", ls="", ms=12, mfc="none", mec="#008800", mew=1.6, label="protected atom")], loc="lower center", ncol=7, fontsize=RL.FS, frameon=False)
        png = f"{OUT}/{r['state_id']}.png"; fig.savefig(png, dpi=100); plt.close(fig)
        cap = f"{r['state_id']} | {c['cls']} {c['target_environment']} | {r['preferred_cell']}, {r['new_n_atoms']} Au | protected {len(prot)} | 6x6: {r['six_status']} {r['six_reasons']}"
        cards.append(f'<figure><a href="{html.escape(os.path.basename(png))}"><img loading="lazy" src="{html.escape(os.path.basename(png))}"></a><figcaption>{html.escape(cap)}</figcaption></figure>')
    # parents after MD (last frame), top views
    M = json.load(open(f"{T}/parents/parents_manifest.json")); fig, axes = plt.subplots(2, 3, figsize=(27, 16))
    for ax, m in zip(axes.ravel(), M):
        try: at = trajio.frame_by_step(m["parent_id"], trajio.TOTAL_STEPS)
        except Exception: ax.set_axis_off(); ax.set_title(f"{m['parent_id']}: no final frame"); continue
        P = at.get_positions(); lay = layers_from_z(at); cc = at.get_cell().array
        for i in np.argsort(P[:, 2]):
            if lay[i] >= 2: ax.add_patch(Circle(P[i, :2], 0.48 * 2.94, fc=RL.LAYER_COL.get(int(lay[i]), "#000"), ec="k", lw=0.2))
        ax.set_xlim(-2, cc[0][0] + cc[1][0] + 2); ax.set_ylim(-2, cc[1][1] + 2); ax.set_aspect("equal"); ax.set_axis_off()
        ax.set_title(f"{m['parent_id']} after 200 ps: layers {np.bincount(lay).tolist()}", fontsize=14)
    fig.tight_layout(); fig.savefig(f"{OUT}/_parents_after_md.png", dpi=80); plt.close(fig)
    page = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><title>Complex gallery</title><style>
body{{font:14px "DejaVu Sans",sans-serif;margin:16px}} figure{{margin:0 0 18px}} img{{width:100%;max-width:1600px;border:1px solid #ccc}} figcaption{{color:#555;font-size:12px}}</style></head><body>
<h1>Multi-layer parents and representative cuts ({len(rows)} candidates that passed the screen)</h1>
<p>Each row: parent frame with the cut frame → zoom with the protected target atoms ringed → reconstructed cell (2×2 tiling) and side view. Click an image for full size.</p>
<figure><a href="_parents_after_md.png"><img src="_parents_after_md.png"></a><figcaption>the six parents after 200 ps (top views, layers ≥ 2)</figcaption></figure>
{''.join(cards)}</body></html>"""
    open(f"{OUT}/index.html", "w").write(page); print(f"wrote {OUT}/index.html with {len(rows)} cards")


if __name__ == "__main__":
    main()
