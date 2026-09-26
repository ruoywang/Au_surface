#!/usr/bin/env python3
"""Where does the induced charge sit across the step, on the metal side and on the
electrolyte side?  For each pair (reference potential A -> perturbed potential B):

  metal side   : d sigma_e(x) = -sum_{z,y} [rho_e^B - rho_e^A] dV / (projected column area)
                 from the CHGCAR total-density blocks (CHGCAR stores rho*V on the fine grid,
                 so sum(rho*V)/Ngrid = N_e; the electron density lives on the metal, the
                 implicit ions are NOT in CHGCAR).  Sign: positive = electrons removed.
  fluid side   : d q_ion(x) = sum_{z,y} d(n_anion - n_cation) dV / (projected column area)
                 from PHI/SION (same reconstruction as step_width_excess.py), all z.
  potential    : S_ion-weighted mean potential shift per column in the near-surface window,
                 d psi_S(x) = sum S*(phi_B - phi_A) / sum S  -- the shift the ions actually
                 feel, unlike the unweighted column mean plotted in step_width_profile.py.

Checks: sum_x d sigma_e * column area = -dN_e (CPM) and sum_x d q_ion * column area = -dN_e.
Both profiles are plotted on the same axis, per projected area.  The metal and fluid grids
differ (fine FFT grid vs solvent grid), so each is resampled to fractional x before plotting.

Usage (from 03_pilot/):  ../scripts/pyrun.sh ../scripts/step_induced_charge_profile.py
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
    ("Step-8x1", "Step-8x1_muref", "Step-8x1_dUp02", 396.014037, 395.923404, "#a83d2f"),
    ("Step-16x1", "Step-16x1_muref_fastcfg", "Step-16x1_dUp02_fastcfg", 792.038710, 791.862697, "#2f6fa8"),
]


def read_grid_block(path):
    """First 3D block of a VASP grid file (CHGCAR / PHI / SION ...): returns (ngz,ngy,ngx) array.
    Reads exactly ngx*ngy*ngz values, so CHGCAR augmentation blocks are never touched."""
    import itertools
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) == 3 and all(p.isdigit() for p in parts):
                ngx, ngy, ngz = map(int, parts); n = ngx * ngy * ngz
                break
        first = next(f)
        per_line = len(first.split())
        n_lines = -(-n // per_line)                       # ceil
        chunk = first + "".join(itertools.islice(f, n_lines - 1))
    arr = np.fromstring(chunk, sep=" ")[:n]
    assert arr.size == n, f"{path}: read {arr.size} values, expected {n}"
    return arr.reshape((ngz, ngy, ngx))


def species(phi, sion):
    u = np.clip(phi / KBT, -50, 50)
    D = 1 + (2 * N_BULK / N_MAX) * (np.cosh(u) - 1)
    return sion * N_BULK * np.exp(-u) / D, sion * N_BULK * np.exp(u) / D


def geometry(D):
    atoms = read(f"{D}/POSCAR")
    cell = atoms.get_cell()
    a2v = np.array(cell[1][:2]); a1v = np.array([cell[0][0], cell[0][1]])
    t_hat = a2v / np.linalg.norm(a2v)
    L_perp = np.linalg.norm(a1v - np.dot(a1v, t_hat) * t_hat)
    V = atoms.get_volume(); area = V / Lz
    ztop = atoms.get_positions()[:, 2].max()
    frac = atoms.get_scaled_positions(wrap=True)
    top = frac[np.abs(atoms.get_positions()[:, 2] - ztop) < 0.3, 0]
    return dict(V=V, area=area, L_perp=L_perp, ztop=ztop, upper_u=(float(top.min()), float(top.max())))


out = {}
fig, axes = plt.subplots(2, 2, figsize=(12, 7.5), dpi=200, sharey="row",
                         gridspec_kw=dict(width_ratios=[1, 2]))
for col, (label, dA, dB, NeA, NeB, color) in enumerate(PAIRS):
    g = geometry(dA)
    dNe = NeB - NeA
    # --- metal side: CHGCAR (rho*V on the fine grid) ---
    rA = read_grid_block(f"{dA}/CHGCAR"); rB = read_grid_block(f"{dB}/CHGCAR")
    assert rA.shape == rB.shape, "CHGCAR grids differ between the two states"
    ngz_c, ngy_c, ngx_c = rA.shape
    NeA_chg = rA.sum() / rA.size; NeB_chg = rB.sum() / rB.size
    d_e_col = (rB - rA).sum(axis=(0, 1)) / rA.size          # electrons gained per column (e)
    dsig_e = -d_e_col / (g["area"] / ngx_c)                   # induced +charge per projected area
    # --- fluid side: PHI/SION ---
    phiA = read_grid_block(f"{dA}/PHI"); phiB = read_grid_block(f"{dB}/PHI")
    sA = read_grid_block(f"{dA}/SION"); sB = read_grid_block(f"{dB}/SION")
    ngz_s, ngy_s, ngx_s = phiA.shape
    dV_s = g["V"] / phiA.size
    anA, catA = species(phiA, sA); anB, catB = species(phiB, sB)
    dq_ion_col = ((anB - catB) - (anA - catA)).sum(axis=(0, 1)) * dV_s
    dq_ion = dq_ion_col / (g["area"] / ngx_s)
    dG_an = (anB - anA).sum(axis=(0, 1)) * dV_s / (g["area"] / ngx_s)
    dG_cat = (catB - catA).sum(axis=(0, 1)) * dV_s / (g["area"] / ngx_s)
    # --- S-weighted potential shift per column, near-surface window ---
    z = np.arange(ngz_s) / ngz_s * Lz
    sz = sA.mean(axis=(1, 2))
    zm = (sz > 0.01) & (z > g["ztop"]) & (z < g["ztop"] + NEAR_SURFACE_DEPTH)
    w = sA[zm]
    dpsi_S = (w * (phiB - phiA)[zm]).sum(axis=(0, 1)) / w.sum(axis=(0, 1))
    dpsi_unw = (phiB - phiA)[zm].mean(axis=(0, 1))

    s_c = np.arange(ngx_c) / ngx_c * g["L_perp"]
    s_s = np.arange(ngx_s) / ngx_s * g["L_perp"]
    lo, hi = g["upper_u"]
    upper_c = (np.arange(ngx_c) / ngx_c >= lo - 0.5 / ngx_c) & (np.arange(ngx_c) / ngx_c <= hi + 0.5 / ngx_c)
    upper_s = (np.arange(ngx_s) / ngx_s >= lo) & (np.arange(ngx_s) / ngx_s <= hi)
    lower_s = (np.arange(ngx_s) / ngx_s > 0.5) & (np.arange(ngx_s) / ngx_s < 1.0)
    rec = dict(pair=[dA, dB], dN_e_CPM=float(dNe),
               dN_e_CHGCAR=float(NeB_chg - NeA_chg),
               N_e_CHGCAR=[float(NeA_chg), float(NeB_chg)],
               sum_dsig_e_x_area=float(dsig_e.sum() * g["area"] / ngx_c),
               sum_dq_ion_x_area=float(dq_ion.sum() * g["area"] / ngx_s),
               metal_charge_share_upper_terrace_columns=float(dsig_e[upper_c].sum() / dsig_e.sum()),
               fraction_of_columns_upper=float(upper_c.mean()),
               ion_countercharge_share_upper_terrace_columns=float(dq_ion[upper_s].sum() / dq_ion.sum()),
               dsig_e_peak_s_A=float(s_c[np.argmax(dsig_e)]), dsig_e_peak=float(dsig_e.max()),
               dq_ion_peak_s_A=float(s_s[np.argmax(dq_ion)]), dq_ion_peak=float(dq_ion.max()),
               dpsi_S_mean_upper_V=float(dpsi_S[upper_s].mean()), dpsi_S_mean_lower_V=float(dpsi_S[lower_s].mean()),
               dpsi_unweighted_mean_upper_V=float(dpsi_unw[upper_s].mean()),
               dpsi_unweighted_mean_lower_V=float(dpsi_unw[lower_s].mean()),
               dG_anion_mean_upper=float(dG_an[upper_s].mean()), dG_anion_mean_lower=float(dG_an[lower_s].mean()),
               s_metal_A=[float(v) for v in s_c], dsig_e=[float(v) for v in dsig_e],
               s_fluid_A=[float(v) for v in s_s], dq_ion=[float(v) for v in dq_ion],
               dpsi_S_V=[float(v) for v in dpsi_S])
    out[label] = rec
    print(f"== {label}: CHGCAR grid {ngx_c}x{ngy_c}x{ngz_c}, solvent grid {ngx_s}x{ngy_s}x{ngz_s}")
    print(f"   N_e from CHGCAR: A={NeA_chg:.6f} B={NeB_chg:.6f}  dN_e(CHGCAR)={NeB_chg-NeA_chg:+.6f} vs CPM {dNe:+.6f}")
    print(f"   closure: sum dsig_e*area={rec['sum_dsig_e_x_area']:+.6f}  sum dq_ion*area={rec['sum_dq_ion_x_area']:+.6f}  (expect {-dNe:+.6f})")
    print(f"   raised-terrace columns ({upper_c.mean()*100:.0f}% of width) carry {rec['metal_charge_share_upper_terrace_columns']*100:.1f}% "
          f"of the induced METAL charge and {rec['ion_countercharge_share_upper_terrace_columns']*100:.1f}% of the ionic countercharge")
    print(f"   induced metal charge peak at s={rec['dsig_e_peak_s_A']:.2f} A ({rec['dsig_e_peak']*1e3:.3f} e-3 e/A^2); "
          f"ionic countercharge peak at s={rec['dq_ion_peak_s_A']:.2f} A ({rec['dq_ion_peak']*1e3:.3f} e-3 e/A^2)")
    print(f"   S-weighted dpsi: raised {rec['dpsi_S_mean_upper_V']*1e3:+.3f} mV, lower {rec['dpsi_S_mean_lower_V']*1e3:+.3f} mV   "
          f"| unweighted: raised {rec['dpsi_unweighted_mean_upper_V']*1e3:+.3f} mV, lower {rec['dpsi_unweighted_mean_lower_V']*1e3:+.3f} mV")
    print(f"   dGamma_- mean: raised {rec['dG_anion_mean_upper']*1e3:.4f}, lower {rec['dG_anion_mean_lower']*1e3:.4f} (e-3 e/A^2)")

    ax = axes[0, col]
    ax.plot(s_c, dsig_e * 1e3, "-", lw=1.4, color="#444444", label="induced metal charge  d(sigma_e)")
    ax.plot(s_s, dq_ion * 1e3, "-", lw=1.6, color=color, label="ionic countercharge  d(q_ion) = dGamma_- - dGamma_+")
    ax.plot(s_s, dG_an * 1e3, "--", lw=1.0, color=color, alpha=0.7, label="anion part  dGamma_-")
    ax.plot(s_s, -dG_cat * 1e3, ":", lw=1.0, color=color, alpha=0.7, label="cation part  -dGamma_+")
    ax.axhline(0, color="#999", lw=0.6)
    for ue in (0.0, 0.5):
        ax.axvline(ue * g["L_perp"], color="#333", lw=0.8, ls="--")
    ax.axvspan(lo * g["L_perp"], hi * g["L_perp"], color="#d9a441", alpha=0.18, label="raised terrace (atoms)")
    ax.set_title(f"{label}: terrace {g['L_perp']/2:.1f} A/side"); ax.grid(alpha=0.25)
    if col == 0:
        ax.set_ylabel("induced charge per projected area (1e-3 e/A^2)")
    ax.legend(fontsize=7, loc="upper left")
    ax = axes[1, col]
    ax.plot(s_s, dpsi_S * 1e3, "-", lw=1.4, color=color, label="S_ion-weighted  d(psi)")
    ax.plot(s_s, dpsi_unw * 1e3, "--", lw=1.0, color="#777777", label="unweighted column mean (for reference)")
    for ue in (0.0, 0.5):
        ax.axvline(ue * g["L_perp"], color="#333", lw=0.8, ls="--")
    ax.axvspan(lo * g["L_perp"], hi * g["L_perp"], color="#d9a441", alpha=0.18)
    ax.set_xlabel("across-step position s (A, perpendicular to the edge)"); ax.grid(alpha=0.25)
    if col == 0:
        ax.set_ylabel("potential shift in the near-surface fluid (mV)")
    ax.legend(fontsize=8, loc="lower right")
plt.suptitle("Induced charge across the step: metal side (CHGCAR difference) vs electrolyte side (ion reconstruction)", y=0.995)
plt.tight_layout()
plt.savefig("report_assets/batch1/step_induced_charge_profile.png", dpi=200)
json.dump(out, open("step_induced_charge_profile.json", "w"), indent=1)
print("wrote report_assets/batch1/step_induced_charge_profile.png and step_induced_charge_profile.json")
