#!/usr/bin/env python3
"""Render all 8 defect-type structures, large, for the report."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ase.io import read

SRC = "03_pilot/all_defect_structures"
OUT = "03_pilot/report_assets"
os.makedirs(OUT, exist_ok=True)

LAYER_COLORS = ["#7a5a1e", "#9c7526", "#c2952f", "#e3b545"]
DEFECT_ADD = "#2f8f8a"
VACANCY_RING = "#c2483d"

TAGS = {
    "T": (None, "none"), "V1": (1, "vac"), "V2": (2, "vac"), "V3": (3, "vac"),
    "V7": (7, "vac"), "A1_fcc": (1, "add"), "A1_hcp": (1, "add"), "A3": (3, "add"),
}

# reference full (non-defective) top-layer site sets for each supercell size, to diff against
_ref_cache = {}


def ref_top_xy(nx_ny_key, cell):
    if nx_ny_key in _ref_cache:
        return _ref_cache[nx_ny_key]
    from ase.build import fcc111
    nx, ny = nx_ny_key
    ref = fcc111('Au', size=(nx, ny, 4), a=4.158, vacuum=None, orthogonal=False, periodic=True)
    p = ref.get_positions()
    top_z = p[:, 2].max()
    xy = p[np.isclose(p[:, 2], top_z, atol=0.05)][:, :2]
    _ref_cache[nx_ny_key] = xy
    return xy


for tag, (n_defect, kind) in TAGS.items():
    atoms = read(f"{SRC}/{tag}.poscar")
    pos = atoms.get_positions()
    cell = atoms.get_cell()
    layer_z_all = np.sort(np.unique(np.round(pos[:, 2], 2)))
    # detect if an adatom layer exists (a layer above the "normal" 4-layer stack, sparsely populated)
    is_add = kind == "add"
    normal_layers = layer_z_all[:4] if not is_add else layer_z_all[:4]
    top_normal_z = normal_layers[-1]

    # supercell key from atom density in a layer (approx): count atoms at bottom layer -> nx*ny
    bottom_n = np.sum(np.isclose(pos[:, 2], layer_z_all[0], atol=0.05))
    nkey = (6, 6) if bottom_n > 20 else (4, 4)

    colors, sizes = [], []
    adatom_positions = []
    for p in pos:
        if is_add and p[2] > top_normal_z + 0.3:
            colors.append(DEFECT_ADD); sizes.append(320)
            adatom_positions.append(p[:2])
            continue
        li = int(np.argmin(np.abs(normal_layers - p[2])))
        colors.append(LAYER_COLORS[min(li, 3)])
        sizes.append(190 + li * 12)

    vac_xy = []
    if kind == "vac":
        ref_xy = ref_top_xy(nkey, cell)
        cur_top_xy = pos[np.isclose(pos[:, 2], top_normal_z, atol=0.05)][:, :2]
        for xy in ref_xy:
            if not np.any(np.all(np.isclose(cur_top_xy, xy, atol=0.15), axis=1)):
                vac_xy.append(xy)

    fig, axes = plt.subplots(1, 2, figsize=(11, 5.6), dpi=230,
                              gridspec_kw={"width_ratios": [1.35, 1]})
    fig.patch.set_alpha(0)
    pad = 2.4 if nkey == (4, 4) else 3.2

    ax = axes[0]
    ax.set_facecolor("none")
    for i in range(-1, 2):
        for j in range(-1, 2):
            shift = i * cell[0][:2] + j * cell[1][:2]
            if i == 0 and j == 0:
                continue
            ax.scatter(pos[:, 0] + shift[0], pos[:, 1] + shift[1], c="none", s=[s * 0.6 for s in sizes],
                       edgecolors="#00000012", linewidths=0.6, zorder=1)
    order = np.argsort(pos[:, 2])
    ax.scatter(pos[order, 0], pos[order, 1], c=[colors[k] for k in order], s=[sizes[k] * 0.62 for k in order],
               edgecolors="#00000045", linewidths=0.6, zorder=3)
    for k, xy in enumerate(vac_xy):
        ax.scatter([xy[0]], [xy[1]], s=210, facecolors="none", edgecolors=VACANCY_RING, linewidths=2.6, zorder=5)
    if vac_xy:
        cx, cy = np.mean(vac_xy, axis=0)
        ax.annotate(f"{len(vac_xy)} vacanc{'y' if len(vac_xy)==1 else 'ies'}", (cx, cy),
                    xytext=(cx + 4.6, cy + 4.0), fontsize=13, color=VACANCY_RING, fontweight="bold",
                    arrowprops=dict(arrowstyle="-", color=VACANCY_RING, lw=1.2))
    if adatom_positions:
        cx, cy = np.mean(adatom_positions, axis=0)
        label = f"{len(adatom_positions)} adatom{'s' if len(adatom_positions)>1 else ''}"
        ax.annotate(label, (cx, cy), xytext=(cx + 4.6, cy + 4.0), fontsize=13, color=DEFECT_ADD,
                    fontweight="bold", arrowprops=dict(arrowstyle="-", color=DEFECT_ADD, lw=1.2))
    c0 = np.array([0, 0]); c1 = cell[0][:2]; c2 = cell[1][:2]; c3 = c1 + c2
    outline = np.array([c0, c1, c3, c2, c0])
    ax.plot(outline[:, 0], outline[:, 1], color="#00000050", lw=1.2, ls=(0, (4, 3)), zorder=6)
    ax.set_xlim(pos[:, 0].min() - pad, pos[:, 0].max() + pad + 6)
    ax.set_ylim(pos[:, 1].min() - pad, pos[:, 1].max() + pad)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("top view", fontsize=13, color="#6b6b6b", pad=10)

    ax2 = axes[1]
    ax2.set_facecolor("none")
    order2 = np.argsort(pos[:, 0])
    ax2.scatter(pos[order2, 0], pos[order2, 2], c=[colors[k] for k in order2], s=[sizes[k] * 0.75 for k in order2],
               edgecolors="#00000045", linewidths=0.6, zorder=3)
    ax2.set_xlim(pos[:, 0].min() - pad, pos[:, 0].max() + pad)
    ax2.set_ylim(pos[:, 2].min() - 2.5, pos[:, 2].max() + 5)
    ax2.set_aspect("equal"); ax2.axis("off")
    ax2.set_title("side view", fontsize=13, color="#6b6b6b", pad=10)

    plt.tight_layout()
    plt.savefig(f"{OUT}/{tag}.png", transparent=True, dpi=230)
    plt.close()
    print(tag, "rendered,", len(atoms), "atoms, cell", nkey)
