#!/usr/bin/env python3
"""One large figure per class showing how a rough state is made:
    parent surface (32 x 32, after MD)  ->  cut frame with core and periphery  ->  reconstructed periodic cell (tiled)
Atoms are discs coloured by layer, exactly as in the gallery; nothing is drawn that is not in the coordinates.

Usage (from Au_Cl/):  rough_sampling_v1/pyrun_rs.sh rough_sampling_v1/render_rough_lineage.py [--cells cells] [--manifest rough200/rough200_manifest.json] [--out figures]
Output: rough_sampling_v1/<out>/rough_lineage_<class>.png (one per class) and rough_lineage_all.png
"""
import argparse
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
R = f"{ROOT}/rough_sampling_v1"
sys.path.insert(0, R)
from extract_cells import lattice, boundary_offset, neighbours, layers_from_z, R_CORE  # noqa: E402
import trajio  # noqa: E402

A0 = 4.158; D111 = A0 / np.sqrt(3); RAD = 0.62 * A0 / np.sqrt(2) / 1.0
LAYER_COL = {0: "#4a4a4a", 1: "#7a7a7a", 2: "#b9b9b9", 3: "#f2c14e", 4: "#e0603a", 5: "#8b1a1a"}
FS = 15
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": FS})


def discs(ax, at, sel=None, alpha=1.0, edge="k", lw=0.3, shift=(0.0, 0.0)):
    P = at.get_positions(); lay = layers_from_z(at)
    idx = np.arange(len(P)) if sel is None else sel
    order = idx[np.argsort(P[idx, 2])]
    for i in order:
        ax.add_patch(Circle((P[i, 0] + shift[0], P[i, 1] + shift[1]), RAD, fc=LAYER_COL.get(int(lay[i]), "#000"), ec=edge, lw=lw, alpha=alpha))


def frame_polygon(origin, v1, v2):
    return np.array([origin, origin + v1, origin + v1 + v2, origin + v2])


def panel_parent(ax, parent, cut_origin, v1, v2, centre_xy, title):
    cell = parent.get_cell().array
    discs(ax, parent, lw=0.15)
    poly = frame_polygon(cut_origin, v1, v2)
    ax.add_patch(Polygon(poly, closed=True, fill=False, ec="#0033cc", lw=3.0))
    ax.add_patch(Circle(centre_xy, R_CORE, fill=False, ec="#008800", lw=2.5, ls="--"))
    ax.add_patch(Polygon(frame_polygon(np.zeros(2), cell[0][:2], cell[1][:2]), closed=True, fill=False, ec="k", lw=1.0, ls=":"))
    lo = np.min([np.zeros(2), cell[0][:2], cell[1][:2], cell[0][:2] + cell[1][:2]], axis=0) - 3
    hi = np.max([np.zeros(2), cell[0][:2], cell[1][:2], cell[0][:2] + cell[1][:2]], axis=0) + 3
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_aspect("equal"); ax.set_title(title, fontsize=FS + 1)
    ax.set_xlabel("x (Å)"); ax.set_ylabel("y (Å)")


def panel_cut(ax, parent, cut_origin, v1, v2, centre_idx, core_ids, title):
    """Zoom on the cut: atoms inside the frame drawn fully, outside faded; core ringed in green, periphery plain."""
    P = parent.get_positions()
    S = np.array([v1, v2]); g = (P[:, :2] - cut_origin) @ np.linalg.inv(S)
    inside = (g[:, 0] >= 0) & (g[:, 0] < 1) & (g[:, 1] >= 0) & (g[:, 1] < 1)
    margin = 6.0
    lo = cut_origin - margin; hi = cut_origin + v1 + v2 + margin
    near = (P[:, 0] > min(lo[0], cut_origin[0] + v2[0] - margin)) & (P[:, 0] < max(hi[0], cut_origin[0] + v1[0] + margin)) & (P[:, 1] > lo[1]) & (P[:, 1] < hi[1])
    discs(ax, parent, sel=np.flatnonzero(near & ~inside), alpha=0.25, edge="none")
    discs(ax, parent, sel=np.flatnonzero(near & inside), lw=0.4)
    for i in core_ids:
        ax.add_patch(Circle(P[i, :2], RAD * 0.55, fill=False, ec="#008800", lw=1.6))
    ax.add_patch(Circle(P[centre_idx, :2], RAD * 0.45, fc="#008800", ec="none"))
    ax.add_patch(Circle(P[centre_idx, :2], R_CORE, fill=False, ec="#008800", lw=2.0, ls="--"))
    ax.add_patch(Polygon(frame_polygon(cut_origin, v1, v2), closed=True, fill=False, ec="#0033cc", lw=3.0))
    ax.set_xlim(min(lo[0], cut_origin[0] + v2[0] - margin), max(hi[0], cut_origin[0] + v1[0] + margin)); ax.set_ylim(lo[1], hi[1])
    ax.set_aspect("equal"); ax.set_title(title, fontsize=FS + 1); ax.set_xlabel("x (Å)")


def panel_cell(ax, cell_at, centre_in_cell, title):
    """The reconstructed cell tiled 2 x 2 so the seams can be seen; the cell outline on the original copy."""
    c = cell_at.get_cell().array
    for i in (0, 1):
        for j in (0, 1):
            sh = i * c[0][:2] + j * c[1][:2]
            discs(ax, cell_at, alpha=1.0 if (i, j) == (0, 0) else 0.55, lw=0.4, shift=tuple(sh))
    P = cell_at.get_positions()
    ax.add_patch(Circle(P[centre_in_cell, :2], R_CORE, fill=False, ec="#008800", lw=2.5, ls="--"))
    ax.add_patch(Polygon(frame_polygon(np.zeros(2), c[0][:2], c[1][:2]), closed=True, fill=False, ec="#0033cc", lw=3.0))
    pts = np.array([[0, 0], c[0][:2], c[1][:2], c[0][:2] + c[1][:2]]) * 2
    ax.set_xlim(pts[:, 0].min() - 3, pts[:, 0].max() + 3); ax.set_ylim(pts[:, 1].min() - 3, pts[:, 1].max() + 3)
    ax.set_aspect("equal"); ax.set_title(title, fontsize=FS + 1); ax.set_xlabel("x (Å)")


def parent_frame(state):
    pid = state["parent_id"]
    if state["source"] == "md":
        return trajio.frame_by_step(pid, state["step"])
    return read(f"{R}/parents/{pid}.extxyz")


def draw_state(fig, axes, state, cell_row):
    parent = parent_frame(state); cell_at = read(f"{state['cell_dir']}/cell.extxyz")
    a1, a2, cellp = lattice(parent); off = boundary_offset(parent, a1, a2)
    n1, n2 = (int(x) for x in state["cell"].split("x")); oij = cell_row["origin"]
    cut_origin = (oij[0] + off[0]) * a1 + (oij[1] + off[1]) * a2
    # wrap the cut origin into the parent cell for drawing, and the centre with it
    Cp = np.array([cellp[0][:2], cellp[1][:2]]); f = cut_origin @ np.linalg.inv(Cp); f -= np.floor(f); cut_origin = f @ Cp
    v1, v2 = n1 * a1, n2 * a2
    P = parent.get_positions(); centre = state["centre_atom"] if "centre_atom" in state else cell_row["atom"]
    core_ids, _ = neighbours(parent, centre, R_CORE)
    # draw the parent copy whose image of the centre lies inside the frame
    S = np.array([v1, v2]); best = None
    for si in (-1, 0, 1):
        for sj in (-1, 0, 1):
            cxy = P[centre, :2] + si * cellp[0][:2] + sj * cellp[1][:2]; g = (cxy - cut_origin) @ np.linalg.inv(S)
            if 0 <= g[0] < 1 and 0 <= g[1] < 1: best = (si, sj)
    if best is None: best = (0, 0)
    shift = best[0] * cellp[0][:2] + best[1] * cellp[1][:2]
    par = parent.copy(); par.set_positions(P + np.r_[shift, 0.0])
    # wrap all parent atoms into the parent cell for the overview, but keep the shifted copy for the zoom
    panel_parent(axes[0], parent, cut_origin - shift, v1, v2, P[centre, :2],
                 f"{state['cls']}: parent {state['parent_id']}, step {state['step']} = {state.get('time_ps', 0):.0f} ps ({len(parent)} Au)\n32 x 32 cell; blue = cut frame {state['cell']}, green = 6 Å core")
    panel_cut(axes[1], par, cut_origin, v1, v2, centre, core_ids,
              f"cut: {len(core_ids) + 1} core atoms (ringed) kept to 1e-3 Å,\nperiphery = the rest of the frame; faded = outside")
    panel_cell(axes[2], cell_at, cell_row["centre_in_cell"],
               f"reconstructed cell {state['cell']} ({len(cell_at)} Au), tiled 2 x 2\n{cell_row['repair']}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", default="cells"); ap.add_argument("--manifest", default="rough200/rough200_manifest.json"); ap.add_argument("--out", default="figures")
    a = ap.parse_args()
    out = f"{R}/{a.out}"; os.makedirs(out, exist_ok=True)
    rows = {json.loads(l)["cell_id"]: json.loads(l) for l in open(f"{R}/{a.cells}/cells_manifest.jsonl")}
    man = json.load(open(f"{R}/{a.manifest}"))
    picked = {}
    for s in man["states"]:
        if s["cls"] not in picked and s["cell_id"] in rows: picked[s["cls"]] = s
    classes = [c for c in "ABCD" if c in picked]
    fig, axes = plt.subplots(len(classes), 3, figsize=(27, 8.6 * len(classes)), squeeze=False)
    for k, cls in enumerate(classes):
        s = picked[cls]; row = rows[s["cell_id"]]
        draw_state(fig, axes[k], s, row)
    handles = [plt.Line2D([], [], marker="o", ls="", ms=14, mfc=LAYER_COL[k], mec="k", label=f"layer {k}" + (" (fixed)" if k < 2 else " (terrace)" if k == 3 else " (adatom level)" if k == 4 else " (pit floor)" if k == 2 else "")) for k in range(5)]
    handles += [plt.Line2D([], [], marker="o", ls="", ms=14, mfc="none", mec="#008800", mew=1.6, label="core atom (kept exactly)"),
                plt.Line2D([], [], ls="--", color="#008800", lw=2, label="6 Å protection radius"), plt.Line2D([], [], ls="-", color="#0033cc", lw=3, label="cut frame = new periodic cell")]
    fig.legend(handles=handles, loc="lower center", ncol=8, fontsize=FS, frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle("Rough sampling: parent surface -> cut frame (core kept exactly, periphery repaired only if needed) -> reconstructed periodic electrode", fontsize=FS + 3, y=0.995)
    fig.tight_layout(rect=(0, 0.02, 1, 0.98))
    fig.savefig(f"{out}/rough_lineage_all.png", dpi=110); print(f"wrote {out}/rough_lineage_all.png ({len(classes)} classes)")
    for k, cls in enumerate(classes):
        f2, ax2 = plt.subplots(1, 3, figsize=(27, 8.6)); s = picked[cls]; draw_state(f2, ax2, s, rows[s["cell_id"]])
        f2.legend(handles=handles, loc="lower center", ncol=8, fontsize=FS, frameon=False); f2.tight_layout(rect=(0, 0.05, 1, 1))
        f2.savefig(f"{out}/rough_lineage_{cls}.png", dpi=110); plt.close(f2)


if __name__ == "__main__":
    main()
