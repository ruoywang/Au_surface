#!/usr/bin/env python3
"""Per-column anion-response profile across the step for the Step-8x1 and Step-16x1
pairs (reference potential -> perturbed potential, same construction, same statistics
as step_potential_mechanism.py). For each x-column (across-step direction) in the
near-surface region (metal top to +15 A):

    K(x) = sum_z,y n_anion / (n_bulk * sum_z,y S_ion)          (whole-column S-weighted ratio)
    dK(x) = K_B(x) - K_A(x)

plotted against the true perpendicular across-step position s = u*L_perp (u = fractional
coordinate along a1), with the two step edges and the raised terrace marked -- NOT folded
onto a single edge distance (rev 1 did that and it merged the upper/lower sides of the two
edges into one zigzag). Purpose: see whether the terrace-interior response at 16x1 reaches
a plateau (a genuine 'far-from-edge' value) and how far the edge influence extends -- the
8x1 'far' region is only a ~4 A wide strip and may not be a converged terrace interior.

Usage (from 03_pilot/):  ../scripts/pyrun.sh ../scripts/step_width_profile.py
"""
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ase.io import read

KBT = 8.6173857e-5 * 298.0
N_BULK = 1.0 * 6.022e-4
R_ION = 4.0
N_MAX = 1.0 / (2 ** (5 / 6) * R_ION) ** 3
Lz = 44.603
NEAR_SURFACE_DEPTH = 15.0

PAIRS = [
    ("Step-8x1", "Step-8x1_muref", "Step-8x1_dUp02", "#a83d2f"),
    ("Step-16x1", "Step-16x1_muref_fastcfg", "Step-16x1_dUp02_fastcfg", "#2f6fa8"),
]


def read_field(path):
    with open(path) as f:
        lines = f.readlines()
    for i, l in enumerate(lines):
        parts = l.split()
        if len(parts) == 3 and all(p.isdigit() for p in parts):
            ngx, ngy, ngz = map(int, parts); dim_idx = i; break
    flat = np.fromstring("".join(lines[dim_idx + 1:]), sep=" ")[:ngx * ngy * ngz]
    return flat.reshape((ngz, ngy, ngx)), ngx, ngy, ngz


def g_occupation(u):
    u = np.clip(u, -50, 50)
    D = 1 + (2 * N_BULK / N_MAX) * (np.cosh(u) - 1)
    return np.exp(-u) / D


def load(D):
    atoms = read(f"{D}/POSCAR")
    cell = atoms.get_cell()
    a1 = cell[0][0]
    a2v = np.array(cell[1][:2]); a1v = np.array([cell[0][0], cell[0][1]])
    t_hat = a2v / np.linalg.norm(a2v)
    L_perp = np.linalg.norm(a1v - np.dot(a1v, t_hat) * t_hat)
    phi, ngx, ngy, ngz = read_field(f"{D}/PHI")
    sion, *_ = read_field(f"{D}/SION")
    z = (np.arange(ngz) / ngz) * Lz
    u = np.arange(ngx) / ngx                      # fractional coordinate along a1 (grid axis)
    ztop = atoms.get_positions()[:, 2].max()
    frac = atoms.get_scaled_positions(wrap=True)
    top = frac[np.abs(atoms.get_positions()[:, 2] - ztop) < 0.3, 0]   # raised-terrace atoms' u
    return dict(phi=phi, sion=sion, z=z, u=u, a1=a1, L_perp=L_perp, ngx=ngx, ztop=ztop,
                upper_u=(float(top.min()), float(top.max())))


def edge_distance(u, edges_frac, L_perp):
    d = np.full_like(u, np.inf)
    for ue in edges_frac:
        d = np.minimum(d, np.abs(((u - ue + 0.5) % 1.0) - 0.5))
    return d * L_perp


def column_K(state):
    sion_z = state["sion"].mean(axis=(1, 2))
    zmask = (sion_z > 0.01) & (state["z"] > state["ztop"]) & (state["z"] < state["ztop"] + NEAR_SURFACE_DEPTH)
    n = state["sion"] * N_BULK * g_occupation(state["phi"] / KBT)
    num = n[zmask].sum(axis=(0, 1))
    den = state["sion"][zmask].sum(axis=(0, 1))
    return num / (N_BULK * den), zmask


out = {}
fig, axes = plt.subplots(2, 2, figsize=(12, 7.5), dpi=200, sharey="row",
                         gridspec_kw=dict(width_ratios=[1, 2]))
for col, (label, dA, dB, color) in enumerate(PAIRS):
    A, B = load(dA), load(dB)
    assert A["ngx"] == B["ngx"]
    KA, zm = column_K(A); KB, _ = column_K(B)
    dK = KB - KA
    s = A["u"] * A["L_perp"]
    d = edge_distance(A["u"], [0.0, 0.5], A["L_perp"])
    dpsi = (B["phi"] - A["phi"])[zm].mean(axis=(0, 1))
    interior = d > 6.0
    plateau = float(dK[interior].mean()) if interior.any() else float("nan")
    # peak location
    ipk = int(np.argmax(dK))
    out[label] = dict(pair=[dA, dB], terrace_half_width_A=float(A["L_perp"] / 2),
                      upper_terrace_u=A["upper_u"],
                      dK_at_edge_mean_lt1A=float(dK[d < 1.0].mean()),
                      dK_interior_gt6A=plateau, n_columns_interior=int(interior.sum()),
                      dK_min=float(dK.min()), dK_max=float(dK.max()),
                      dK_peak_position_A=float(s[ipk]), dK_peak_edge_distance_A=float(d[ipk]),
                      s_A=[float(v) for v in s], dK=[float(v) for v in dK], dpsi_V=[float(v) for v in dpsi])
    print(f"{label}: half-terrace={A['L_perp']/2:.2f} A  raised terrace u in [{A['upper_u'][0]:.3f},{A['upper_u'][1]:.3f}]  "
          f"dK(edge, d<1A)={out[label]['dK_at_edge_mean_lt1A']:.4f}  dK(interior, d>6A)={plateau:.4f} "
          f"over {interior.sum()} columns  dK range=[{dK.min():.4f},{dK.max():.4f}]  "
          f"peak at s={s[ipk]:.2f} A (d_edge={d[ipk]:.2f} A)")
    for row, y, ylab in [(0, dK, "dK(x) = K_B(x) - K_A(x)  (per column, near-surface)"),
                         (1, dpsi, "mean d(psi) per column (V, unweighted)")]:
        ax = axes[row, col]
        ax.plot(s, y, "-", lw=1.4, color=color)
        for ue in (0.0, 0.5):
            ax.axvline(ue * A["L_perp"], color="#333", lw=0.8, ls="--")
        lo, hi = A["upper_u"]
        ax.axvspan(lo * A["L_perp"], hi * A["L_perp"], color="#d9a441", alpha=0.18,
                   label="raised terrace (atoms)" if row == 0 else None)
        ax.grid(alpha=0.25)
        if col == 0:
            ax.set_ylabel(ylab)
        if row == 0:
            ax.set_title(f"{label}: terrace {A['L_perp']/2:.1f} A/side, edges at dashed lines")
            ax.legend(fontsize=8, loc="lower right")
        else:
            ax.set_xlabel("across-step position s (A, perpendicular to the edge)")
    axes[0, col].axhline(plateau if not np.isnan(plateau) else np.nan, color=color, lw=0.8, ls=":")
plt.suptitle("Anion enrichment response to the -0.2 eV mu_e step, resolved across the step", y=0.995)
plt.tight_layout()
plt.savefig("report_assets/batch1/step_width_dK_profile.png", dpi=200)
json.dump(out, open("step_width_profile.json", "w"), indent=1)
print("wrote report_assets/batch1/step_width_dK_profile.png and step_width_profile.json")
