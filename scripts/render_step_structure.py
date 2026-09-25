#!/usr/bin/env python3
"""Render the phase-2 step-edge concept structure, same visual style as
render_all_defects.py (depth-graded colors, top+side view, dashed cell
outline, periodic ghost neighbors)."""
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
STRIP_COLOR = "#2f8f8a"
STEP_MARK = "#1f6a66"

atoms = read(f"{SRC}/Step.poscar")
pos = atoms.get_positions()
cell = atoms.get_cell()
layer_z = np.sort(np.unique(np.round(pos[:, 2], 2)))  # includes the 5th (strip) level
base_layers = layer_z[:4]
strip_z = layer_z[4] if len(layer_z) > 4 else None

colors, sizes = [], []
for p in pos:
    if strip_z is not None and np.isclose(p[2], strip_z, atol=0.05):
        colors.append(STRIP_COLOR); sizes.append(210)
    else:
        li = int(np.argmin(np.abs(base_layers - p[2])))
        colors.append(LAYER_COLORS[min(li, 3)])
        sizes.append(190 + li * 12)

fig, axes = plt.subplots(1, 2, figsize=(11, 5.6), dpi=230, gridspec_kw={"width_ratios": [1.35, 1]})
fig.patch.set_alpha(0)
pad = 2.4

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

# step-edge boundaries: strip occupies fractional a1 in [0, 0.5); mark both edges (u=0 and u=0.5)
a1 = cell[0][:2]; a2 = cell[1][:2]
for u, label_dx in [(0.0, -1), (0.5, 1)]:
    p0 = u * a1
    p1 = u * a1 + a2
    ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=STEP_MARK, lw=2.2, ls=(0, (1, 1)), zorder=4)
c0 = np.array([0, 0]); c1 = cell[0][:2]; c2 = cell[1][:2]; c3 = c1 + c2
outline = np.array([c0, c1, c3, c2, c0])
ax.plot(outline[:, 0], outline[:, 1], color="#00000050", lw=1.2, ls=(0, (4, 3)), zorder=6)
ax.set_xlim(pos[:, 0].min() - pad, pos[:, 0].max() + pad + 6)
ax.set_ylim(pos[:, 1].min() - pad - 3.2, pos[:, 1].max() + pad)
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
ax2.set_title("side view — plateau profile", fontsize=13, color="#6b6b6b", pad=10)

fig.text(0.5, 0.015,
          "teal = extra Au layer covering half the top face  →  raised strip bounded by two step edges (dotted)",
          ha="center", fontsize=12, color=STEP_MARK)

plt.tight_layout(rect=[0, 0.05, 1, 1])
plt.savefig(f"{OUT}/Step.png", transparent=True, dpi=230)
plt.close()
print("Step rendered,", len(atoms), "atoms")
