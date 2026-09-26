#!/usr/bin/env python3
"""Cross-width comparison on a window-independent footing.

The enrichment ratio K_D = sum(n_anion)/(n_bulk sum(S_ion)) is a *volume-diluted*
quantity: its value depends on how much bulk-like accessible fluid the fixed 15 A
near-surface window contains above the local SION boundary. That is harmless for a
near/far comparison inside one cell (same window), but it is NOT a clean way to compare
the flat T cell (metal top 12.20 A, SION boundary ~15 A) against the step cells (metal
top 14.60 A, corrugated boundary 18-20.6 A): the deeper accessible column in T dilutes
K differently. So here every response is expressed as an anion SURFACE EXCESS per
projected area, which converges once the window contains the whole double layer:

    Gamma_-(x) = sum_z,y [n_anion - n_bulk*S_ion] * dz*dy / a2      [e / A^2, per x-column]
    dGamma_-(x) = Gamma_-^B(x) - Gamma_-^A(x)

and the cell total  dQ_ion = sum over cell of d(n_anion - n_cation) dV  is checked against
the CPM electron-count change -dN_e (charge closure) so the numbers are anchored to the
code's own bookkeeping rather than to the analysis window.

Usage (from 03_pilot/):  ../scripts/pyrun.sh ../scripts/step_width_excess.py
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
EDGE_BAND = 3.0

# (label, dir_A, dir_B, N_e^A, N_e^B) -- electron counts from the CPM-ion lines
PAIRS = [
    ("T (flat)", "T_dUp00", "T_dUm02", 704.000000, 703.835517, "#555555"),
    ("Step-8x1", "Step-8x1_muref", "Step-8x1_dUp02", 396.014037, 395.923404, "#a83d2f"),
    ("Step-16x1", "Step-16x1_muref_fastcfg", "Step-16x1_dUp02_fastcfg", 792.038710, 791.862697, "#2f6fa8"),
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


def species(phi, sion):
    u = np.clip(phi / KBT, -50, 50)
    D = 1 + (2 * N_BULK / N_MAX) * (np.cosh(u) - 1)
    return sion * N_BULK * np.exp(-u) / D, sion * N_BULK * np.exp(u) / D


def load(D):
    atoms = read(f"{D}/POSCAR")
    cell = atoms.get_cell()
    a1 = cell[0][0]
    a2v = np.array(cell[1][:2]); a1v = np.array([cell[0][0], cell[0][1]])
    t_hat = a2v / np.linalg.norm(a2v)
    L_perp = np.linalg.norm(a1v - np.dot(a1v, t_hat) * t_hat)
    phi, ngx, ngy, ngz = read_field(f"{D}/PHI")
    sion, *_ = read_field(f"{D}/SION")
    V = atoms.get_volume()
    dV = V / (ngx * ngy * ngz)
    area = V / Lz
    n_an, n_cat = species(phi, sion)
    z = (np.arange(ngz) / ngz) * Lz
    x = (np.arange(ngx) / ngx) * a1
    ztop = atoms.get_positions()[:, 2].max()
    sz = sion.mean(axis=(1, 2))
    zmask = (sz > 0.01) & (z > ztop) & (z < ztop + NEAR_SURFACE_DEPTH)
    return dict(n_an=n_an, n_cat=n_cat, sion=sion, phi=phi, z=z, x=x, a1=a1, L_perp=L_perp,
                dV=dV, area=area, ngx=ngx, zmask=zmask, ztop=ztop)


def edge_distance(x, a1, edges_frac, L_perp):
    u = x / a1
    d = np.full_like(u, np.inf)
    for ue in edges_frac:
        d = np.minimum(d, np.abs(((u - ue + 0.5) % 1.0) - 0.5))
    return d * L_perp


def gamma_columns(st):
    """anion surface excess per projected area, per x column (e/A^2), near-surface window."""
    excess = (st["n_an"] - N_BULK * st["sion"])[st["zmask"]]      # e/A^3
    col = excess.sum(axis=(0, 1)) * st["dV"]                        # e per column
    return col / (st["area"] / st["ngx"])                           # e/A^2


def upper_terrace_u(D):
    atoms = read(f"{D}/POSCAR")
    ztop = atoms.get_positions()[:, 2].max()
    frac = atoms.get_scaled_positions(wrap=True)
    top = frac[np.abs(atoms.get_positions()[:, 2] - ztop) < 0.3, 0]
    return float(top.min()), float(top.max())


out = {}
fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), dpi=200, sharey=True,
                         gridspec_kw=dict(width_ratios=[1, 2]))
T_ref = None
print("charge closure: sum d(n_anion - n_cation) dV over the FULL cell vs -dN_e from CPM")
for label, dA, dB, NeA, NeB, color in PAIRS:
    A, B = load(dA), load(dB)
    dQ_ion = ((B["n_an"] - B["n_cat"]) - (A["n_an"] - A["n_cat"])).sum() * A["dV"]
    dNe = NeB - NeA
    dsig = -dNe / A["area"]
    gA, gB = gamma_columns(A), gamma_columns(B)
    dG = gB - gA
    dG_all = dG.mean()                                   # cell-average excess response, e/A^2
    dG_an_full = ((B["n_an"] - A["n_an"]).sum() * A["dV"]) / A["area"]   # full-cell anion part
    dG_cat_full = ((B["n_cat"] - A["n_cat"]).sum() * A["dV"]) / A["area"]
    rec = dict(pair=[dA, dB], area_A2=float(A["area"]), dN_e=float(dNe),
               induced_sigma_e_per_A2=float(dsig),
               dQ_ion_recon=float(dQ_ion), closure_rel_err=float((dQ_ion + dNe) / abs(dNe)),
               dGamma_anion_window_e_per_A2=float(dG_all),
               dGamma_anion_fullcell_e_per_A2=float(dG_an_full),
               dGamma_cation_fullcell_e_per_A2=float(dG_cat_full),
               anion_share_of_countercharge=float(dG_an_full / (dG_an_full - dG_cat_full)))
    print(f"{label:10s} area={A['area']:7.2f} A^2  dN_e={dNe:+.6f}  dQ_ion(recon)={dQ_ion:+.6f}  "
          f"closure err={rec['closure_rel_err']:+.1e}  induced sigma={dsig*1e3:+.4f} e-3 e/A^2  "
          f"dGamma_-(window)={dG_all*1e3:+.4f}  dGamma_-(full)={dG_an_full*1e3:+.4f}  "
          f"dGamma_+(full)={dG_cat_full*1e3:+.4f}  anion share={rec['anion_share_of_countercharge']:.3f}")
    if "Step" in label:
        d = edge_distance(A["x"], A["a1"], [0.0, 0.5], A["L_perp"])
        near = d < EDGE_BAND; far = ~near
        interior = d > 6.0
        rec.update(dGamma_near=float(dG[near].mean()), dGamma_far=float(dG[far].mean()),
                   ratio_near_far=float(dG[near].mean() / dG[far].mean()),
                   dGamma_interior_gt6A=float(dG[interior].mean()) if interior.any() else None,
                   dGamma_edge_lt1A=float(dG[d < 1.0].mean()))
        s = (A["x"] / A["a1"]) * A["L_perp"]
        ax = axes[0 if label == "Step-8x1" else 1]
        ax.plot(s, dG * 1e3, "-", lw=1.4, color=color, label=f"{label}")
        for ue in (0.0, 0.5):
            ax.axvline(ue * A["L_perp"], color="#333", lw=0.8, ls="--")
        lo, hi = upper_terrace_u(dA)
        ax.axvspan(lo * A["L_perp"], hi * A["L_perp"], color="#d9a441", alpha=0.18, label="raised terrace (atoms)")
        if T_ref is not None:
            ax.axhline(T_ref * 1e3, color="#555555", ls="--", lw=1.2, label="flat T: cell average")
        ax.set_title(f"{label}: terrace {A['L_perp']/2:.1f} A/side")
        ax.set_xlabel("across-step position s (A, perpendicular to the edge)")
        ax.grid(alpha=0.25); ax.legend(fontsize=8, loc="lower right")
        # binned print
        bins = np.arange(0, A["L_perp"] / 2 + 1e-6, 2.0)
        idx = np.digitize(d, bins) - 1
        prof = [f"{bins[i]:.0f}-{bins[i]+2:.0f}A:{dG[idx==i].mean()*1e3:.3f}" for i in range(len(bins)) if (idx == i).any()]
        print(f"    dGamma_-(x) x1e3 e/A^2 binned by edge distance: " + "  ".join(prof))
        print(f"    near(<{EDGE_BAND:.0f}A)={rec['dGamma_near']*1e3:.4f}  far={rec['dGamma_far']*1e3:.4f}  "
              f"ratio={rec['ratio_near_far']:.3f}  interior(>6A)={ (rec['dGamma_interior_gt6A'] or float('nan'))*1e3:.4f}")
    else:
        T_ref = dG_all
    out[label] = rec

axes[0].set_ylabel("dGamma_- per column  (1e-3 e / A^2, projected area)")
plt.suptitle("Anion surface-excess response to the -0.2 eV mu_e step (window-independent), resolved across the step", y=0.99)
plt.tight_layout()
plt.savefig("report_assets/batch1/step_width_dGamma_profile.png", dpi=200)
json.dump(out, open("step_width_excess.json", "w"), indent=1)
print("wrote report_assets/batch1/step_width_dGamma_profile.png and step_width_excess.json")
