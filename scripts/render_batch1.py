#!/usr/bin/env python3
"""Batch-1 diagrams, revision 2. Fixes vs rev 1:
- retracted the "edge1/edge2 verified NOT equivalent (12 vs 11)" claim --
  replaced with the re-checked, hedged finding from check_step_edges.py
- Au211 caption no longer lists Au(511) as a family member (wrong family)
- Island-7 caption no longer claims the 4x4 case was "geometrically
  impossible" -- that applies only to a 19-atom same-layer island, not 7
- thickness split into base_z_span / all_atoms_z_span, read from manifest
- Pit-7 gets a real thin-slice cross-section (narrow y-band through the
  pit center) instead of a full-width side view that visually fills in
  the gap with front/back rows
- Au211/Au221 side view tiled across 3 periods with a traced upper-surface
  envelope, so the terrace/step staircase is unambiguous
- "periodic-mirror view" renamed "periodic-repeat view" (lattice
  translation copies, not reflections)
- exports both a transparent and a white-background version
"""
import os
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
from ase.io import read
from ase.build import fcc111

SRC = "03_pilot/all_defect_structures"
OUT = "03_pilot/report_assets/batch1"
os.makedirs(OUT, exist_ok=True)
MANIFEST = json.load(open(f"{OUT}/manifest.json"))

ADD_COLOR = "#2f8f8a"
VAC_COLOR = "#c2483d"
EDGE1_COLOR = "#1f6a66"
EDGE2_COLOR = "#a2452e"
LAYER_SHADES = ["#7a5a1e", "#9c7526", "#c2952f", "#e3b545"]


def depth_colors(z, cmap_name="YlOrBr_r", lo=0.15, hi=0.85):
    zn = (z - z.min()) / max(z.max() - z.min(), 1e-6)
    cmap = cm.get_cmap(cmap_name)
    return [cmap(lo + (hi - lo) * v) for v in zn]


def cell_outline_segments(cell2d, i, j):
    a1, a2 = cell2d
    shift = i * a1 + j * a2
    c0, c1, c2, c3 = shift, shift + a1, shift + a1 + a2, shift + a2
    return np.array([c0, c1, c2, c3, c0])


def save_both(fig, tag):
    fig.patch.set_alpha(0)
    plt.savefig(f"{OUT}/{tag}.png", transparent=True, dpi=230)
    fig.patch.set_facecolor("white")
    for ax in fig.axes:
        ax.set_facecolor("white")
    plt.savefig(f"{OUT}/{tag}_whitebg.png", transparent=False, dpi=230, facecolor="white")
    plt.close(fig)
    print("rendered", tag, "(+ whitebg)")


def draw_repeat_panel(ax, pos, colors, sizes, cell2d, tiles):
    a1, a2 = cell2d
    order = np.argsort(pos[:, 2])
    for i, j in tiles:
        shift = i * a1 + j * a2
        ax.scatter(pos[order, 0] + shift[0], pos[order, 1] + shift[1], c=[colors[k] for k in order],
                   s=[sizes[k] for k in order], edgecolors="#00000040", linewidths=0.5, zorder=3)
        outline = cell_outline_segments(cell2d, i, j)
        ax.plot(outline[:, 0], outline[:, 1], color="#00000055", lw=1.0, ls=(0, (4, 3)), zorder=6)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("periodic-repeat view", fontsize=13, color="#6b6b6b", pad=10)


def base_panel_top(ax, pos, colors, sizes, cell2d):
    a1, a2 = cell2d
    for i in range(-1, 2):
        for j in range(-1, 2):
            if i == 0 and j == 0:
                continue
            shift = i * a1 + j * a2
            ax.scatter(pos[:, 0] + shift[0], pos[:, 1] + shift[1], c="none",
                       s=[s * 0.6 for s in sizes], edgecolors="#00000012", linewidths=0.6, zorder=1)
    order = np.argsort(pos[:, 2])
    ax.scatter(pos[order, 0], pos[order, 1], c=[colors[k] for k in order],
               s=[sizes[k] * 0.62 for k in order], edgecolors="#00000045", linewidths=0.6, zorder=3)
    outline = cell_outline_segments(cell2d, 0, 0)
    ax.plot(outline[:, 0], outline[:, 1], color="#00000050", lw=1.2, ls=(0, (4, 3)), zorder=6)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("top view", fontsize=13, color="#6b6b6b", pad=10)


def base_panel_side(ax, pos, colors, sizes):
    order = np.argsort(pos[:, 0])
    ax.scatter(pos[order, 0], pos[order, 2], c=[colors[k] for k in order],
               s=[sizes[k] * 0.75 for k in order], edgecolors="#00000045", linewidths=0.6, zorder=3)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("side view", fontsize=13, color="#6b6b6b", pad=10)


def layer_colors_sizes_4layer(pos, base_layers, top_extra_z=None, extra_color=ADD_COLOR, extra_size=210):
    colors, sizes = [], []
    for p in pos:
        if top_extra_z is not None and np.isclose(p[2], top_extra_z, atol=0.05):
            colors.append(extra_color); sizes.append(extra_size)
        else:
            li = int(np.argmin(np.abs(base_layers - p[2])))
            colors.append(LAYER_SHADES[min(li, 3)]); sizes.append(190 + li * 12)
    return colors, sizes


# ---------------------------------------------------------------- Step-16x4
def do_step16x4():
    m = MANIFEST["Step-16x4"]
    a = read(f"{SRC}/Step-16x4.poscar")
    pos = a.get_positions()
    cell = a.get_cell()
    a1, a2 = cell[0][:2], cell[1][:2]
    zvals = np.round(pos[:, 2], 2)
    z_levels = np.sort(np.unique(zvals))
    base_layers = z_levels[:4]
    strip_z = z_levels[4]
    colors, sizes = layer_colors_sizes_4layer(pos, base_layers, strip_z)

    fig = plt.figure(figsize=(16.5, 5.6), dpi=230)
    ax1 = fig.add_subplot(1, 3, 1); ax2 = fig.add_subplot(1, 3, 2); ax3 = fig.add_subplot(1, 3, 3)
    pad = 2.4

    base_panel_top(ax1, pos, colors, sizes, (a1, a2))
    for u, name, col in [(0.0, "edge 1", EDGE1_COLOR), (0.5, "edge 2", EDGE2_COLOR)]:
        p0 = u * a1; p1 = u * a1 + a2
        ax1.plot([p0[0], p1[0]], [p0[1], p1[1]], color=col, lw=2.4, ls=(0, (1, 1)), zorder=7)
        mid = p0 + 0.5 * (p1 - p0)
        ax1.annotate(name, mid, xytext=(mid[0], mid[1] - 4.5), fontsize=12, color=col, fontweight="bold",
                     ha="center", arrowprops=dict(arrowstyle="-", color=col, lw=1.1))
    ax1.set_xlim(pos[:, 0].min() - pad, pos[:, 0].max() + pad)
    ax1.set_ylim(pos[:, 1].min() - pad - 3, pos[:, 1].max() + pad)

    base_panel_side(ax2, pos, colors, sizes)
    ax2.set_xlim(pos[:, 0].min() - pad, pos[:, 0].max() + pad)
    ax2.set_ylim(pos[:, 2].min() - 2.5, pos[:, 2].max() + 4)

    draw_repeat_panel(ax3, pos, colors, sizes, (a1, a2), [(0, 0), (1, 0)])

    caption = (
        f"Step-16x4  |  decorated-(111)-slab, 4-layer base + half-layer strip  |  N={m['N_total']} "
        f"({m['N_base']} base + {m['N_added']} strip)  |  base_z_span {m['base_z_span']} A, all_atoms_z_span {m['all_atoms_z_span']} A\n"
        f"terrace width {m['terrace_width_each_side']} A/side (perp. to edge)  |  L_par {m['L_parallel']} A, L_perp {m['L_perpendicular']} A\n"
        "edge classification RETRACTED pending re-check -- see check_step_edges.py: edge1_top/edge2_top (actual\n"
        "upper-edge atoms) both have CN=7, identical shell composition (3 below+4 in-layer+0 above) at 3.2-3.6 A\n"
        "cutoffs -- suggestive of equivalence, not confirmed distinct, and not confirmed equivalent either\n"
        "role: constrained_model | relaxation_status: not_run | status: geometry_only (registry bug fixed)"
    )
    fig.text(0.5, 0.02, caption, ha="center", fontsize=9.5, color="#4a4030", family="monospace")
    plt.tight_layout(rect=[0, 0.14, 1, 1])
    save_both(fig, "Step-16x4")


# ---------------------------------------------------------------- Au211/221
def do_vicinal(tag, poscar, facet, family, microfacet, width_axis_index, n_periods):
    m = MANIFEST[tag]
    a = read(f"{SRC}/{poscar}")
    pos = a.get_positions()
    cell = a.get_cell()
    a1, a2 = cell[0][:2], cell[1][:2]
    colors = depth_colors(pos[:, 2])
    sizes = [230] * len(pos)

    fig = plt.figure(figsize=(16.5, 5.6), dpi=230)
    ax1 = fig.add_subplot(1, 3, 1); ax2 = fig.add_subplot(1, 3, 2); ax3 = fig.add_subplot(1, 3, 3)
    pad = 1.5

    base_panel_top(ax1, pos, colors, sizes, (a1, a2))
    ax1.set_xlim(pos[:, 0].min() - pad, pos[:, 0].max() + pad)
    ax1.set_ylim(pos[:, 1].min() - pad, pos[:, 1].max() + pad)

    # cross-section tiled across n_periods along the width (terrace-crossing) axis,
    # with the upper-surface envelope traced so the staircase is unambiguous
    width_vec = [a1, a2][width_axis_index]
    for k in range(n_periods):
        shift = k * width_vec[0]
        ax2.scatter(pos[:, 0] + shift, pos[:, 2], c=colors, s=[s * 0.75 for s in sizes],
                    edgecolors="#00000045", linewidths=0.6, zorder=3)
    # envelope: for each x-bin, the max z (topmost exposed atom)
    all_x = np.concatenate([pos[:, 0] + k * width_vec[0] for k in range(n_periods)])
    all_z = np.tile(pos[:, 2], n_periods)
    bins = np.linspace(all_x.min(), all_x.max(), 60)
    idx = np.digitize(all_x, bins)
    env_x, env_z = [], []
    for b in range(1, len(bins)):
        sel = idx == b
        if sel.any():
            env_x.append(all_x[sel].mean()); env_z.append(all_z[sel].max())
    ax2.plot(env_x, env_z, color="#a2452e", lw=1.6, zorder=5)
    period_len = np.linalg.norm(width_vec)
    ax2.set_title(f"cross-section x{n_periods} (envelope traced, period={period_len:.2f} A)",
                  fontsize=12, color="#6b6b6b", pad=10)
    ax2.set_aspect("equal"); ax2.axis("off")

    draw_repeat_panel(ax3, pos, colors, sizes, (a1, a2), [(i, j) for i in range(3) for j in range(3)])

    caption = (
        f"{tag}  |  true crystallographic vicinal slab, {facet}  |  step family: {family}  |  N={m['N_total']}\n"
        f"all_atoms_z_span {m['all_atoms_z_span']} A (macroscopic-normal height, real geometry not layer count)  |  "
        f"step microfacet: {microfacet}\n"
        "hkl normal aligned to z  |  back termination natural (unmodified)  |  role: constrained_model | "
        "relaxation_status: not_run | status: geometry_only (newly built)"
    )
    fig.text(0.5, 0.02, caption, ha="center", fontsize=9.5, color="#4a4030", family="monospace")
    plt.tight_layout(rect=[0, 0.12, 1, 1])
    save_both(fig, tag)


# ---------------------------------------------------------------- Island-7
def do_island7():
    m = MANIFEST["Island-7"]
    a = read(f"{SRC}/Island-7-6x6.poscar")
    pos = a.get_positions()
    cell = a.get_cell()
    a1, a2 = cell[0][:2], cell[1][:2]
    zvals = np.round(pos[:, 2], 2)
    z_levels = np.sort(np.unique(zvals))
    base_layers = z_levels[:4]
    top_z = z_levels[-1]
    colors, sizes = layer_colors_sizes_4layer(pos, base_layers,
                                               top_z if top_z > base_layers[-1] + 0.3 else None,
                                               extra_size=300)

    fig = plt.figure(figsize=(16.5, 5.6), dpi=230)
    ax1 = fig.add_subplot(1, 3, 1); ax2 = fig.add_subplot(1, 3, 2); ax3 = fig.add_subplot(1, 3, 3)
    pad = 2.4

    base_panel_top(ax1, pos, colors, sizes, (a1, a2))
    ax1.set_xlim(pos[:, 0].min() - pad, pos[:, 0].max() + pad)
    ax1.set_ylim(pos[:, 1].min() - pad, pos[:, 1].max() + pad)
    base_panel_side(ax2, pos, colors, sizes)
    ax2.set_xlim(pos[:, 0].min() - pad, pos[:, 0].max() + pad)
    ax2.set_ylim(pos[:, 2].min() - 2.5, pos[:, 2].max() + 4)
    draw_repeat_panel(ax3, pos, colors, sizes, (a1, a2), [(i, j) for i in range(2) for j in range(2)])

    caption = (
        f"Island-7 (6x6 base)  |  compact 7-atom monolayer island, fcc-continuation registry  |  "
        f"N={m['N_total']} ({m['N_base']} base + {m['N_added']} island)  |  base_z_span {m['base_z_span']} A, "
        f"all_atoms_z_span {m['all_atoms_z_span']} A\n"
        f"nearest lattice-translation copy of the island: {m['periodic_image_separation']:.1f} A away (geometric "
        "disconnection only -- not a claim about EDL/response convergence)\n"
        "note: the earlier 4x4 base was invalid for a 19-atom island (16 sites/layer), NOT for this 7-atom one "
        "(7/16 fits); 6x6 was chosen here to increase island-image separation, checked above\n"
        "role: relaxed_candidate (not yet relaxed) | relaxation_status: not_run | status: geometry_only (newly built)"
    )
    fig.text(0.5, 0.02, caption, ha="center", fontsize=9.5, color="#4a4030", family="monospace")
    plt.tight_layout(rect=[0, 0.14, 1, 1])
    save_both(fig, "Island-7-6x6")


# ---------------------------------------------------------------- Pit-7 (V7)
def do_pit7():
    m = MANIFEST["Pit-7"]
    a = read(f"{SRC}/V7.poscar")
    pos = a.get_positions()
    cell = a.get_cell()
    a1, a2 = cell[0][:2], cell[1][:2]

    ref = fcc111('Au', size=(6, 6, 4), a=4.158, vacuum=None, orthogonal=False, periodic=True)
    pref = ref.get_positions(); pref[:, 2] += (5.0 - pref[:, 2].min())
    zvals = np.round(pos[:, 2], 2)
    z_levels = np.sort(np.unique(zvals))
    top_z = z_levels[-1]
    ref_top_xy = pref[np.isclose(np.round(pref[:, 2], 2), top_z, atol=0.05)][:, :2]
    cur_top_xy = pos[np.isclose(zvals, top_z, atol=0.05)][:, :2]
    vac_xy = np.array([xy for xy in ref_top_xy if not np.any(np.all(np.isclose(cur_top_xy, xy, atol=0.15), axis=1))])
    pit_center = vac_xy.mean(axis=0)

    colors, sizes = [], []
    for p in pos:
        li = int(np.argmin(np.abs(z_levels - p[2])))
        colors.append(LAYER_SHADES[min(li, 3)]); sizes.append(190 + li * 12)

    fig = plt.figure(figsize=(16.5, 5.6), dpi=230)
    ax1 = fig.add_subplot(1, 3, 1); ax2 = fig.add_subplot(1, 3, 2); ax3 = fig.add_subplot(1, 3, 3)
    pad = 2.4

    base_panel_top(ax1, pos, colors, sizes, (a1, a2))
    for xy in vac_xy:
        ax1.scatter([xy[0]], [xy[1]], s=230, facecolors="none", edgecolors=VAC_COLOR, linewidths=2.6, zorder=5)
    ax1.annotate("O = removed reference site\n(no atom here)", pit_center, xytext=(pit_center[0] + 5, pit_center[1] - 5),
                 fontsize=9.5, color=VAC_COLOR, ha="left")
    ax1.set_xlim(pos[:, 0].min() - pad, pos[:, 0].max() + pad)
    ax1.set_ylim(pos[:, 1].min() - pad, pos[:, 1].max() + pad)

    # thin slice through the pit center: keep only atoms within +-2A of the
    # pit centroid's y, so front/back rows don't visually fill in the gap
    slice_mask = np.abs(pos[:, 1] - pit_center[1]) < 2.0
    ps = pos[slice_mask]; cs = [colors[i] for i in range(len(pos)) if slice_mask[i]]
    sz = [sizes[i] for i in range(len(pos)) if slice_mask[i]]
    order = np.argsort(ps[:, 0])
    ax2.scatter(ps[order, 0], ps[order, 2], c=[cs[k] for k in order], s=[sz[k] * 0.85 for k in order],
                edgecolors="#00000045", linewidths=0.6, zorder=3)
    # mark where the missing top atoms would have been, within this slice's y-band
    vac_in_slice = vac_xy[np.abs(vac_xy[:, 1] - pit_center[1]) < 2.0]
    for xy in vac_in_slice:
        ax2.scatter([xy[0]], [top_z], s=260, facecolors="none", edgecolors=VAC_COLOR, linewidths=2.2,
                    linestyles="dashed", zorder=6)
    ax2.set_aspect("equal"); ax2.axis("off")
    ax2.set_title("thin slice through pit center (|dy|<2A)", fontsize=12, color="#6b6b6b", pad=10)
    ax2.set_xlim(pos[:, 0].min() - pad, pos[:, 0].max() + pad)
    ax2.set_ylim(pos[:, 2].min() - 2.5, pos[:, 2].max() + 4)

    tiles = [(i, j) for i in range(2) for j in range(2)]
    order_all = np.argsort(pos[:, 2])
    for i, j in tiles:
        shift = i * a1 + j * a2
        ax3.scatter(pos[order_all, 0] + shift[0], pos[order_all, 1] + shift[1],
                    c=[colors[k] for k in order_all], s=[sizes[k] * 0.62 for k in order_all],
                    edgecolors="#00000040", linewidths=0.5, zorder=3)
        for xy in vac_xy:
            ax3.scatter([xy[0] + shift[0]], [xy[1] + shift[1]], s=230, facecolors="none",
                        edgecolors=VAC_COLOR, linewidths=2.2, zorder=5)
        outline = cell_outline_segments((a1, a2), i, j)
        ax3.plot(outline[:, 0], outline[:, 1], color="#00000055", lw=1.0, ls=(0, (4, 3)), zorder=6)
    ax3.set_aspect("equal"); ax3.axis("off")
    ax3.set_title("periodic-repeat view", fontsize=13, color="#6b6b6b", pad=10)

    caption = (
        f"Pit-7 (= V7, reclassified)  |  6x6 base, center + 6 first-shell neighbors removed, depth 1 layer  |  "
        f"N={m['N_total']} ({m['N_base']} base minus {m['N_removed']} removed)  |  base_z_span {m['base_z_span']} A\n"
        f"local metal thickness under the pit floor: ~{2*2.401:.2f} A (2 remaining layers), not the full "
        f"{m['base_z_span']} A base span\n"
        f"nearest lattice-translation copy of the pit footprint: {m['periodic_image_separation']:.1f} A away "
        "(geometric disconnection only)  |  keeps original 'seven-vacancy cluster' construction tag\n"
        "role: relaxed_candidate (not yet relaxed) | relaxation_status: not_run | status: geometry_only (reused, no rebuild)"
    )
    fig.text(0.5, 0.02, caption, ha="center", fontsize=9, color="#4a4030", family="monospace")
    plt.tight_layout(rect=[0, 0.14, 1, 1])
    save_both(fig, "Pit-7")


if __name__ == "__main__":
    do_step16x4()
    do_vicinal("Au211", "Au211.poscar", "Au(211) = 3(111)x(100)", "n(111)x(100)", "{100} -> A-type",
               width_axis_index=0, n_periods=3)
    do_vicinal("Au221", "Au221.poscar", "Au(221) = n(111)x(111)", "n(111)x(111), with Au(332)/Au(554)",
               "{111} -> B-type", width_axis_index=0, n_periods=3)
    do_island7()
    do_pit7()
