#!/usr/bin/env python3
"""Render bigger, clearer top-view and side-view PNGs of T/V1/A1_fcc for the report."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ase.io import read

OUT = "/anvil/scratch/x-rywang/Au_Cl/03_pilot/report_assets"
os.makedirs(OUT, exist_ok=True)

# depth gradient: layer 0 (bottom, fixed) -> darkest, top layer -> lightest gold
LAYER_COLORS = ["#7a5a1e", "#9c7526", "#c2952f", "#e3b545"]
DEFECT_ADD = "#2f8f8a"
VACANCY_RING = "#c2483d"

STRUCTS = {
    "T": ("03_pilot/structures/T/POSCAR", None),
    "V1": ("03_pilot/structures/V1/POSCAR", "vacancy"),
    "A1_fcc": ("03_pilot/structures/A1_fcc/POSCAR", "adatom"),
}

atoms_T = read("/anvil/scratch/x-rywang/Au_Cl/03_pilot/structures/T/POSCAR")
pos_T = atoms_T.get_positions()
top_z_T = pos_T[:, 2].max()


def layer_colors_sizes(atoms, kind, layer_z):
    pos = atoms.get_positions()
    colors, sizes = [], []
    for p in pos:
        li = int(np.argmin(np.abs(layer_z - p[2])))
        if kind == "adatom" and np.isclose(p[2], layer_z[-1] + (layer_z[-1] - layer_z[-2]), atol=0.3):
            colors.append(DEFECT_ADD); sizes.append(300)
        else:
            colors.append(LAYER_COLORS[min(li, len(LAYER_COLORS) - 1)])
            sizes.append(190 + li * 12)
    return colors, sizes


def find_vacancy_xy(atoms_defect):
    pos_d = atoms_defect.get_positions()
    top_z = pos_d[:, 2].max()
    top_T_xy = pos_T[np.isclose(pos_T[:, 2], top_z_T, atol=0.05)][:, :2]
    top_d_xy = pos_d[np.isclose(pos_d[:, 2], top_z, atol=0.05)][:, :2]
    for xy in top_T_xy:
        if not np.any(np.all(np.isclose(top_d_xy, xy, atol=0.1), axis=1)):
            return xy
    return None


def find_adatom_xy(atoms):
    pos = atoms.get_positions()
    return pos[np.argmax(pos[:, 2])][:2]


for tag, (path, kind) in STRUCTS.items():
    atoms = read(path)
    pos = atoms.get_positions()
    cell = atoms.get_cell()
    layer_z = np.sort(np.unique(np.round(pos[:, 2], 2)))
    if kind == "adatom":
        layer_z = layer_z[:-1]  # exclude the adatom's own z from the depth scale
    colors, sizes = layer_colors_sizes(atoms, kind, layer_z)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 5.2), dpi=230,
                              gridspec_kw={"width_ratios": [1.35, 1]})
    fig.patch.set_alpha(0)
    pad = 2.2

    ax = axes[0]
    ax.set_facecolor("none")
    for i in range(-1, 2):
        for j in range(-1, 2):
            shift = i * cell[0][:2] + j * cell[1][:2]
            if i == 0 and j == 0:
                continue
            ax.scatter(pos[:, 0] + shift[0], pos[:, 1] + shift[1], c="none", s=[s * 0.62 for s in sizes],
                       edgecolors="#00000014", linewidths=0.7, zorder=1)
    order = np.argsort(pos[:, 2])
    ax.scatter(pos[order, 0], pos[order, 1], c=[colors[k] for k in order], s=[sizes[k] * 0.68 for k in order],
               edgecolors="#00000045", linewidths=0.7, zorder=3)
    vac_xy = None
    if kind == "vacancy":
        vac_xy = find_vacancy_xy(atoms)
        if vac_xy is not None:
            ax.scatter([vac_xy[0]], [vac_xy[1]], s=230, facecolors="none",
                       edgecolors=VACANCY_RING, linewidths=2.8, zorder=5)
            ax.annotate("vacancy", (vac_xy[0], vac_xy[1]), xytext=(vac_xy[0] + 4.2, vac_xy[1] + 3.6),
                        fontsize=12, color=VACANCY_RING, fontweight="bold",
                        arrowprops=dict(arrowstyle="-", color=VACANCY_RING, lw=1.2))
    if kind == "adatom":
        axy = find_adatom_xy(atoms)
        ax.annotate("adatom", (axy[0], axy[1]), xytext=(axy[0] + 4.2, axy[1] + 3.6),
                    fontsize=12, color=DEFECT_ADD, fontweight="bold",
                    arrowprops=dict(arrowstyle="-", color=DEFECT_ADD, lw=1.2))
    c0 = np.array([0, 0]); c1 = cell[0][:2]; c2 = cell[1][:2]; c3 = c1 + c2
    outline = np.array([c0, c1, c3, c2, c0])
    ax.plot(outline[:, 0], outline[:, 1], color="#00000055", lw=1.3, ls=(0, (4, 3)), zorder=6)
    ax.set_xlim(pos[:, 0].min() - pad, pos[:, 0].max() + pad + 6)
    ax.set_ylim(pos[:, 1].min() - pad, pos[:, 1].max() + pad)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("top view", fontsize=13, color="#6b6b6b", pad=10)

    ax2 = axes[1]
    ax2.set_facecolor("none")
    order2 = np.argsort(pos[:, 0])
    ax2.scatter(pos[order2, 0], pos[order2, 2], c=[colors[k] for k in order2], s=[sizes[k] * 0.78 for k in order2],
               edgecolors="#00000045", linewidths=0.7, zorder=3)
    ax2.set_xlim(pos[:, 0].min() - pad, pos[:, 0].max() + pad)
    ax2.set_ylim(pos[:, 2].min() - 2.5, pos[:, 2].max() + 5)
    ax2.set_aspect("equal"); ax2.axis("off")
    ax2.set_title("side view", fontsize=13, color="#6b6b6b", pad=10)

    plt.tight_layout()
    plt.savefig(f"{OUT}/{tag}.png", transparent=True, dpi=230)
    plt.close()
    print(tag, "n_atoms=", len(atoms))
