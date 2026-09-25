#!/usr/bin/env python3
"""Rigorous re-check of Step-16x4's two edges. Does NOT use
atoms.get_distances(mic=True) for coordination counting -- that call
returns only the SINGLE nearest periodic image per atom index, which
silently misses cases where a small unit cell puts two different periodic
copies of the same base atom within the cutoff. Instead this enumerates
every (atom_j, in-plane lattice shift) pair explicitly, so the neighbor
list is auditable atom-by-atom."""
import numpy as np
from ase.io import read

a = read("03_pilot/all_defect_structures/Step-16x4.poscar")
pos = a.get_positions()
cell = a.get_cell()
a1, a2 = cell[0][:2], cell[1][:2]
M = np.array([a1, a2]).T
Minv = np.linalg.inv(M)
uv = (Minv @ pos[:, :2].T).T
u = uv[:, 0] % 1.0

zvals = np.round(pos[:, 2], 2)
zlevels = np.sort(np.unique(zvals))
base_top_z = zlevels[3]
strip_z = zlevels[4]
is_base_top = np.isclose(zvals, base_top_z, atol=0.05)
is_strip = np.isclose(zvals, strip_z, atol=0.05)


def full_neighbors(idx, cutoff, shift_range=2):
    p0 = pos[idx]
    out = []
    for j in range(len(a)):
        for i in range(-shift_range, shift_range + 1):
            for k in range(-shift_range, shift_range + 1):
                if j == idx and i == 0 and k == 0:
                    continue
                shift = i * a1 + k * a2
                pj_xy = pos[j][:2] + shift
                d = np.sqrt((pj_xy[0] - p0[0]) ** 2 + (pj_xy[1] - p0[1]) ** 2 + (pos[j][2] - p0[2]) ** 2)
                if d < cutoff:
                    out.append((j, (i, k), round(float(d), 4)))
    return sorted(out, key=lambda x: x[2])


base_top_idx = np.where(is_base_top)[0]
strip_idx = np.where(is_strip)[0]
u_bt = u[base_top_idx]
u_st = u[strip_idx]

groups = {
    "flat_control (base-top, mid-terrace u~0.75)": base_top_idx[np.argmin(np.abs(u_bt - 0.75))],
    "edge1_foot (base-top, u~1.0, just before strip starts)": base_top_idx[np.argmin(np.abs(u_bt - 1.0))],
    "edge1_top (strip, u~0.0, just after strip starts)": strip_idx[np.argmin(u_st)],
    "edge2_foot (base-top, u~0.5+, just after strip ends)": base_top_idx[np.argmin(np.abs(u_bt - 0.52))],
    "edge2_top (strip, u~0.5-, just before strip ends)": strip_idx[np.argmax(u_st)],
}

NN = 2.94
cutoffs = [3.2, 3.4, 3.6]

for name, idx in groups.items():
    print(f"\n=== {name}  (atom {idx}, u={u[idx]:.3f}, z={pos[idx,2]:.2f}) ===")
    for c in cutoffs:
        nbrs = full_neighbors(idx, c)
        print(f"  cutoff {c} A: {len(nbrs)} neighbors")
    # full listing at the middle cutoff, for audit
    nbrs = full_neighbors(idx, 3.4)
    for j, shift, d in nbrs:
        dz = pos[j, 2] - pos[idx, 2]
        layer_tag = "strip" if is_strip[j] else ("base-top" if is_base_top[j] else "base-lower")
        print(f"    atom {j:3d} shift={shift}  d={d:.3f}  dz={dz:+.2f}  ({layer_tag})")
