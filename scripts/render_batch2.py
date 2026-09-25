#!/usr/bin/env python3
"""Batch 2 diagrams: Au332/Au554 vicinal slabs, Island-7/19 and Pit-7/19 in
the shared 8x8 cell (size comparison), and the Step-24x4 render (cheap,
piggybacking on the same series as Step-8x4/16x4)."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ase.io import read
from ase.build import fcc111

from render_batch1 import (SRC, OUT, LAYER_SHADES, ADD_COLOR, VAC_COLOR, EDGE1_COLOR, EDGE2_COLOR,
                            save_both, base_panel_top, base_panel_side, draw_repeat_panel,
                            layer_colors_sizes_4layer, do_vicinal, cell_outline_segments)

MANIFEST = json.load(open(f"{OUT}/manifest.json"))


# ------------------------------------------------------------- Step-24x4
def do_step24x4():
    m = MANIFEST["Step-24x4"]
    a = read(f"{SRC}/Step-24x4.poscar")
    pos = a.get_positions()
    cell = a.get_cell()
    a1, a2 = cell[0][:2], cell[1][:2]
    zvals = np.round(pos[:, 2], 2)
    z_levels = np.sort(np.unique(zvals))
    base_layers = z_levels[:4]
    strip_z = z_levels[4]
    colors, sizes = layer_colors_sizes_4layer(pos, base_layers, strip_z)

    fig = plt.figure(figsize=(18, 5.6), dpi=230)
    ax1 = fig.add_subplot(1, 3, 1); ax2 = fig.add_subplot(1, 3, 2); ax3 = fig.add_subplot(1, 3, 3)
    pad = 2.4
    base_panel_top(ax1, pos, colors, sizes, (a1, a2))
    for u, name, col in [(0.0, "edge 1", EDGE1_COLOR), (0.5, "edge 2", EDGE2_COLOR)]:
        p0 = u * a1; p1 = u * a1 + a2
        ax1.plot([p0[0], p1[0]], [p0[1], p1[1]], color=col, lw=2.0, ls=(0, (1, 1)), zorder=7)
    ax1.set_xlim(pos[:, 0].min() - pad, pos[:, 0].max() + pad)
    ax1.set_ylim(pos[:, 1].min() - pad - 3, pos[:, 1].max() + pad)

    base_panel_side(ax2, pos, colors, sizes)
    ax2.set_xlim(pos[:, 0].min() - pad, pos[:, 0].max() + pad)
    ax2.set_ylim(pos[:, 2].min() - 2.5, pos[:, 2].max() + 4)

    draw_repeat_panel(ax3, pos, colors, sizes, (a1, a2), [(0, 0), (1, 0)])

    caption = (
        f"Step-24x4  |  decorated-(111)-slab, widest of the 3-point terrace series  |  N={m['N_total']} "
        f"({m['N_base']} base + {m['N_added']} strip)  |  base_z_span {m['base_z_span']} A, all_atoms_z_span {m['all_atoms_z_span']} A\n"
        f"terrace width {m['terrace_width_each_side']} A/side  |  L_par {m['L_parallel']} A, L_perp {m['L_perpendicular']} A  |  "
        "convergence-series candidate, not required for Batch 1/2 completion\n"
        "role: constrained_model | relaxation_status: not_run | status: geometry_only (built, now rendered)"
    )
    fig.text(0.5, 0.02, caption, ha="center", fontsize=9.5, color="#4a4030", family="monospace")
    plt.tight_layout(rect=[0, 0.1, 1, 1])
    save_both(fig, "Step-24x4")


# ------------------------------------------------------------- Island/Pit size comparison
def do_island_pit(tag, kind, ref_size):
    m = MANIFEST[tag]
    a = read(f"{SRC}/{tag}.poscar")
    pos = a.get_positions()
    cell = a.get_cell()
    a1, a2 = cell[0][:2], cell[1][:2]
    zvals = np.round(pos[:, 2], 2)
    z_levels = np.sort(np.unique(zvals))
    base_layers = z_levels[:4]
    top_z = z_levels[-1]

    vac_xy = None
    if kind == "add":
        colors, sizes = layer_colors_sizes_4layer(pos, base_layers,
                                                   top_z if top_z > base_layers[-1] + 0.3 else None, extra_size=260)
    else:
        colors, sizes = [], []
        for p in pos:
            li = int(np.argmin(np.abs(base_layers - p[2])))
            colors.append(LAYER_SHADES[min(li, 3)]); sizes.append(190 + li * 12)
        ref = fcc111('Au', size=ref_size, a=4.158, vacuum=None, orthogonal=False, periodic=True)
        pref = ref.get_positions(); pref[:, 2] += (5.0 - pref[:, 2].min())
        ref_top_xy = pref[np.isclose(np.round(pref[:, 2], 2), top_z, atol=0.05)][:, :2]
        cur_top_xy = pos[np.isclose(zvals, top_z, atol=0.05)][:, :2]
        vac_xy = np.array([xy for xy in ref_top_xy if not np.any(np.all(np.isclose(cur_top_xy, xy, atol=0.15), axis=1))])

    fig = plt.figure(figsize=(16.5, 5.6), dpi=230)
    ax1 = fig.add_subplot(1, 3, 1); ax2 = fig.add_subplot(1, 3, 2); ax3 = fig.add_subplot(1, 3, 3)
    pad = 2.4
    base_panel_top(ax1, pos, colors, sizes, (a1, a2))
    if vac_xy is not None:
        for xy in vac_xy:
            ax1.scatter([xy[0]], [xy[1]], s=200, facecolors="none", edgecolors=VAC_COLOR, linewidths=2.4, zorder=5)
    ax1.set_xlim(pos[:, 0].min() - pad, pos[:, 0].max() + pad)
    ax1.set_ylim(pos[:, 1].min() - pad, pos[:, 1].max() + pad)

    base_panel_side(ax2, pos, colors, sizes)
    ax2.set_xlim(pos[:, 0].min() - pad, pos[:, 0].max() + pad)
    ax2.set_ylim(pos[:, 2].min() - 2.5, pos[:, 2].max() + 4)

    draw_repeat_panel(ax3, pos, colors, sizes, (a1, a2), [(0, 0), (1, 0)])
    if vac_xy is not None:
        shift = a1
        for xy in vac_xy:
            ax3.scatter([xy[0], xy[0] + shift[0]], [xy[1], xy[1] + shift[1]], s=200, facecolors="none",
                        edgecolors=VAC_COLOR, linewidths=2.0, zorder=5)

    n_key = "N_added" if kind == "add" else "N_removed"
    caption = (
        f"{tag}  |  8x8 common cell (shared with the {19 if '7' in tag else 7}-atom counterpart for a true "
        f"size comparison)  |  N={m['N_total']} ({m['N_base']} base {'+ ' if kind=='add' else '- '}{m[n_key]})\n"
        f"nearest lattice-translation-copy distance: {m['periodic_image_separation']:.1f} A (geometric "
        "disconnection only, not an EDL/response convergence claim)\n"
        "role: relaxed_candidate (not yet relaxed) | relaxation_status: not_run | status: geometry_only (new, Batch 2)"
    )
    fig.text(0.5, 0.02, caption, ha="center", fontsize=9.5, color="#4a4030", family="monospace")
    plt.tight_layout(rect=[0, 0.1, 1, 1])
    save_both(fig, tag)


if __name__ == "__main__":
    do_vicinal("Au332", "Au332.poscar", "Au(332) = n(111)x(111)", "n(111)x(111), with Au(221)/Au(554)",
               "{111} -> B-type", width_axis_index=0, n_periods=3)
    do_vicinal("Au554", "Au554.poscar", "Au(554) = n(111)x(111)", "n(111)x(111), with Au(221)/Au(332)",
               "{111} -> B-type", width_axis_index=0, n_periods=3)
    do_step24x4()
    do_island_pit("Island-7-8x8", "add", (8, 8, 4))
    do_island_pit("Island-19-8x8", "add", (8, 8, 4))
    do_island_pit("Pit-7-8x8", "vac", (8, 8, 4))
    do_island_pit("Pit-19-8x8", "vac", (8, 8, 4))
