#!/usr/bin/env python3
"""Cross-structure comparison of the potential-driven ion response on a footing that does
not depend on the analysis window's accessible volume.  Rev 3 (2026-09-27).

Quantity (ion NUMBER excess per projected area, 1e-3 A^-2; multiply by the ion charge for a
charge density):
        dGamma_-  = (1/A_proj) int_Omega [ dn_-  -  n_b dS_ion ] dV
S_ion changes slightly between the two potentials (the cavity follows the electron density),
so the -n_b dS term is kept; the plain anion-count change (1/A) int dn_- dV is reported
separately under its own name.  Electron counts come from the CPM-ion lines of each run.

Two integration conventions, both reported, because the response is NOT window-converged
inside the available electrolyte (it decays with ~3.2 A but sits on a residual far-field
potential offset of -0.2..-0.3 mV, and the far artificial ion window begins at 31.3 A):

  (A) ABSOLUTE:  Omega = { z < Z_UP }, Z_UP = 31.0 A = last plane with S_ion = 1 everywhere,
      the same for every cell.  No lower bound is needed (S_ion = 0 inside the metal; the
      script prints the below-metal contribution as a check).  Column values partition the
      cell total exactly, so this is the convention for totals, closure and the cell mean.
      Its bias: a column whose accessible boundary sits higher (raised terrace, +2.4 A)
      keeps less of the decaying tail than a flat surface (~3-4 % at lambda = 3.2 A).
  (B) BOUNDARY-RELATIVE: per column, Omega = { z_b(x) < z < z_b(x) + D } with z_b(x) the
      local S_ion = 0.5 crossing (y-averaged) and D = 10.4 A (the largest depth the raised
      terrace can accommodate below 31 A).  Every column keeps the same fraction of the
      tail, so regional ratios between cells/regions are comparable; totals are not exact
      partitions.  This is the convention for raised/lower/flat ratios.

A window scan of (A) is printed so the residual z_up dependence is visible.  Until a scan
plateaus, every number is "surface excess within the stated integration range".

Usage (from 03_pilot/):  ../scripts/pyrun.sh ../scripts/step_width_excess.py
"""
import json
import os
import re

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
ION_Z1 = 32.603
Z_UP = 31.0
DEPTH = 10.4
Z_SCAN = [22.0, 24.0, 26.0, 28.0, 29.6, 31.0, 32.0]
EDGE_BAND = 3.0

PAIRS = [
    ("T-old (4x4 cell, Accurate, 3x3x1)", "T_dUp00", "T_dUm02", "flat", "#777777"),
    ("Flat16x1 (Step-16x1 cell, Normal, 1x12x1)", "Flat16x1_muref", "Flat16x1_dUp02", "flat", "#2a7f62"),
    ("Step-8x1", "Step-8x1_muref", "Step-8x1_dUp02", "step", "#a83d2f"),
    ("Step-16x1", "Step-16x1_muref_fastcfg", "Step-16x1_dUp02_fastcfg", "step", "#2f6fa8"),
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


def cpm_electrons(D):
    m = re.findall(r"CPM-ion:.*?N_ele=\s*([-\d.]+)\s+mu_e=\s*([-\d.]+)", open(f"{D}/log.out").read())
    if not m:
        raise RuntimeError(f"{D}: no CPM-ion line -- run not finished?")
    return float(m[-1][0]), float(m[-1][1])


def species(phi, sion):
    u = np.clip(phi / KBT, -50, 50)
    D = 1 + (2 * N_BULK / N_MAX) * (np.cosh(u) - 1)
    return sion * N_BULK * np.exp(-u) / D, sion * N_BULK * np.exp(u) / D


def load(D):
    atoms = read(f"{D}/POSCAR")
    cell = atoms.get_cell()
    a2v = np.array(cell[1][:2]); a1v = np.array([cell[0][0], cell[0][1]])
    t_hat = a2v / np.linalg.norm(a2v)
    L_perp = np.linalg.norm(a1v - np.dot(a1v, t_hat) * t_hat)
    phi, ngx, ngy, ngz = read_field(f"{D}/PHI")
    sion, *_ = read_field(f"{D}/SION")
    V = atoms.get_volume()
    n_an, n_cat = species(phi, sion)
    z = (np.arange(ngz) / ngz) * Lz
    ztop = atoms.get_positions()[:, 2].max()
    frac = atoms.get_scaled_positions(wrap=True)
    top = frac[np.abs(atoms.get_positions()[:, 2] - ztop) < 0.3, 0]
    Ne, mu = cpm_electrons(D)
    # local accessible boundary per column: first z above the metal top where <S>_y >= 0.5
    s_xz = sion.mean(axis=1)                       # (ngz, ngx)
    zb = np.full(ngx, np.nan)
    for i in range(ngx):
        idx = np.where((s_xz[:, i] >= 0.5) & (z > ztop - 3.0) & (z < ION_Z1 - 5))[0]
        zb[i] = z[idx[0]] if len(idx) else np.nan
    return dict(n_an=n_an, n_cat=n_cat, sion=sion, z=z, u=np.arange(ngx) / ngx, L_perp=L_perp,
                dV=V / (ngx * ngy * ngz), area=V / Lz, ngx=ngx, ztop=ztop, zb=zb,
                upper_u=(float(top.min()), float(top.max())), Ne=Ne, mu=mu)


def edge_distance(u, edges_frac, L_perp):
    d = np.full_like(u, np.inf)
    for ue in edges_frac:
        d = np.minimum(d, np.abs(((u - ue + 0.5) % 1.0) - 0.5))
    return d * L_perp


def integrand(A, B, which):
    dn = (B["n_an"] - A["n_an"]) if which == "anion" else (B["n_cat"] - A["n_cat"])
    return dn - N_BULK * (B["sion"] - A["sion"])


def cols_absolute(A, B, z_up, which="anion"):
    zm = A["z"] < z_up
    return integrand(A, B, which)[zm].sum(axis=(0, 1)) * A["dV"] / (A["area"] / A["ngx"])


def cols_boundary(A, B, depth, which="anion"):
    f = integrand(A, B, which)
    out = np.zeros(A["ngx"])
    for i in range(A["ngx"]):
        zm = (A["z"] >= A["zb"][i]) & (A["z"] < A["zb"][i] + depth)
        out[i] = f[zm, :, i].sum() * A["dV"] / (A["area"] / A["ngx"])
    return out


def count_change_cols(A, B, z_up):
    zm = A["z"] < z_up
    return (B["n_an"] - A["n_an"])[zm].sum(axis=(0, 1)) * A["dV"] / (A["area"] / A["ngx"])


def regions(A, s, y):
    lo, hi = A["upper_u"]
    d = edge_distance(A["u"], [0.0, 0.5], A["L_perp"])
    upper = (A["u"] >= lo) & (A["u"] <= hi)
    lower = (A["u"] > 0.5) & (A["u"] < 1.0)
    near = d < EDGE_BAND
    rows = np.array([lo * A["L_perp"], hi * A["L_perp"]])
    low_c = (hi * A["L_perp"] + (lo + 1.0) * A["L_perp"]) / 2
    idx = [i for i in range(1, len(y) - 1) if y[i] >= y[i - 1] and y[i] >= y[i + 1]]
    idx = sorted(idx, key=lambda i: -y[i]); pk = []
    for i in idx:
        if all(abs(s[i] - s[j]) > 5.0 for j in pk): pk.append(i)
        if len(pk) == 2: break
    def dist_to_row(sv):
        c = np.concatenate([rows, rows + A["L_perp"], rows - A["L_perp"]]); return float(sv - c[np.argmin(np.abs(sv - c))])
    return dict(raised_terrace=float(y[upper].mean()), lower_terrace=float(y[lower].mean()),
                lower_terrace_centre_pm3A=float(y[np.abs(s - low_c) < 3].mean()),
                raised_terrace_centre_pm3A=float(y[np.abs(s - rows.mean()) < 3].mean()),
                peaks=[dict(s_A=float(s[i]), value=float(y[i]), offset_from_upper_edge_row_A=dist_to_row(s[i])) for i in pk],
                folded_near_lt3A=float(y[near].mean()), folded_far=float(y[~near].mean()),
                folded_ratio=float(y[near].mean() / y[~near].mean()),
                upper_edge_rows_s_A=rows.tolist())


out = {"conventions": {"A_absolute": f"z < {Z_UP} A for all cells", "B_boundary_relative": f"z_b(x) <= z < z_b(x)+{DEPTH} A, z_b = local S_ion=0.5 crossing"},
       "definition": "dGamma = (1/A_proj) int_Omega [dn_ion - n_b dS_ion] dV, ions per A^2", "z_scan_A": Z_SCAN}
flat_refs = {}
step_recs = []
print(f"(A) absolute z < {Z_UP} A (all cells);  (B) boundary-relative depth {DEPTH} A;  ION_Z1 = {ION_Z1} A\n")
for label, dA, dB, kind, color in PAIRS:
    if not (os.path.exists(f"{dA}/PHI") and os.path.exists(f"{dB}/PHI") and os.path.exists(f"{dB}/log.out")):
        print(f"[skip] {label}: fields not present yet ({dA}, {dB})\n"); continue
    try:
        A, B = load(dA), load(dB)
    except RuntimeError as e:
        print(f"[skip] {e}\n"); continue
    dNe = B["Ne"] - A["Ne"]; dsig = -dNe / A["area"]
    dmu = B["mu"] - A["mu"]                      # the ACTUAL applied mu_e step (FERMICONVERGE=0.01 lets it differ from -0.2 by up to 10 meV)
    dq_ion = ((B["n_an"] - B["n_cat"]) - (A["n_an"] - A["n_cat"])).sum() * A["dV"]
    dS_term = N_BULK * (B["sion"] - A["sion"]).sum() * A["dV"] / A["area"]
    below = A["z"] < A["ztop"]
    below_contrib = integrand(A, B, "anion")[below].sum() * A["dV"] / A["area"]
    gA = cols_absolute(A, B, Z_UP); gA_cat = cols_absolute(A, B, Z_UP, "cation"); nA = count_change_cols(A, B, Z_UP)
    gB = cols_boundary(A, B, DEPTH)
    scan = {f"{zu:g}": float(cols_absolute(A, B, zu).mean()) for zu in Z_SCAN}
    s = A["u"] * A["L_perp"]
    rec = dict(label=label, kind=kind, pair=[dA, dB], area_A2=float(A["area"]), N_e=[A["Ne"], B["Ne"]], mu_e=[A["mu"], B["mu"]],
               dN_e=float(dNe), dmu_e_actual_eV=float(dmu), induced_sigma_e_per_A2=float(dsig),
               induced_sigma_per_eV=float(dsig / abs(dmu)), dQ_ion_fullcell=float(dq_ion),
               closure_rel_err=float((dq_ion + dNe) / abs(dNe)), n_b_dS_term_fullcell_per_A2=float(dS_term),
               below_metal_top_contribution_per_A2=float(below_contrib), metal_top_z=float(A["ztop"]), L_perp_A=float(A["L_perp"]),
               boundary_z_range_A=[float(np.nanmin(A["zb"])), float(np.nanmax(A["zb"]))],
               A_dGamma_anion=float(gA.mean()), A_dGamma_cation=float(gA_cat.mean()), A_anion_count_change=float(nA.mean()),
               A_anion_share_of_countercharge=float(gA.mean() / (gA.mean() - gA_cat.mean())), A_zscan=scan,
               A_window_change_28_to_31_rel=float((scan["31"] - scan["28"]) / scan["28"]),
               B_dGamma_anion=float(gB.mean()))
    print(f"== {label}")
    print(f"   area={A['area']:.2f} A^2  N_e {A['Ne']:.6f} -> {B['Ne']:.6f}  dN_e={dNe:+.6f}  mu_e {A['mu']:.6f} -> {B['mu']:.6f} "
          f"(actual step {dmu:+.4f} eV)  induced sigma={dsig*1e3:+.4f} e-3 e/A^2 = {dsig/abs(dmu)*1e3:.4f} e-3 e/A^2/eV  "
          f"metal top {A['ztop']:.2f} A, boundary z_b {rec['boundary_z_range_A'][0]:.2f}-{rec['boundary_z_range_A'][1]:.2f} A")
    print(f"   closure: sum d(n_- - n_+) dV = {dq_ion:+.6f} vs -dN_e {-dNe:+.6f} (rel {rec['closure_rel_err']:+.1e});  "
          f"n_b int dS dV/A = {dS_term*1e3:+.5f};  below metal top {below_contrib*1e3:+.1e}   (1e-3 A^-2)")
    print(f"   (A) dGamma_- = {gA.mean()*1e3:.4f}  dGamma_+ = {gA_cat.mean()*1e3:+.4f}  anion count change/A = {nA.mean()*1e3:.4f}  "
          f"anion share {rec['A_anion_share_of_countercharge']:.3f}   |  (B) dGamma_- = {gB.mean()*1e3:.4f}")
    print("   (A) z_up scan: " + "  ".join(f"{k}:{v*1e3:.4f}" for k, v in scan.items()) + f"   | 28->31 A: {rec['A_window_change_28_to_31_rel']*100:+.2f}%")
    if kind == "flat":
        rec.update(A_profile_std_rel=float(gA.std() / gA.mean()), B_profile_std_rel=float(gB.std() / gB.mean()))
        print(f"   flat profile std/mean: (A) {rec['A_profile_std_rel']:.4f}  (B) {rec['B_profile_std_rel']:.4f}")
        flat_refs[label] = dict(A=float(gA.mean()), B=float(gB.mean()), color=color, sigma=float(dsig), dmu=float(dmu))
    else:
        rec["A_regions"] = regions(A, s, gA); rec["B_regions"] = regions(A, s, gB)
        rec.update(s_A=[float(v) for v in s], A_columns=[float(v) for v in gA], B_columns=[float(v) for v in gB],
                   z_boundary_columns=[float(v) for v in A["zb"]])
        for conv in ("A", "B"):
            r = rec[f"{conv}_regions"]
            print(f"   ({conv}) raised {r['raised_terrace']*1e3:.4f}  lower {r['lower_terrace']*1e3:.4f}  lower centre {r['lower_terrace_centre_pm3A']*1e3:.4f}  "
                  f"raised centre {r['raised_terrace_centre_pm3A']*1e3:.4f}  folded near/far {r['folded_ratio']:.3f}  peaks "
                  + "; ".join(f"{p['value']*1e3:.4f} @ {p['offset_from_upper_edge_row_A']:+.2f} A from edge row" for p in r["peaks"]))
        step_recs.append((label, A, gA, gB, color))
    out[label] = rec
    print()

print("ratios step/flat.  'raw' = ratio of the responses as computed;  'per eV' = each response first divided by its own actual |dmu_e|\n"
      "(linear-response normalisation; the two pairs' steps differ when FERMICONVERGE=0.01 stops the CP loop short of the target).")
for label, A, gA, gB, color in step_recs:
    r = out[label]
    for fl, ref in flat_refs.items():
        rat = {}
        f_norm = abs(ref["dmu"]) / abs(r["dmu_e_actual_eV"])       # multiply a raw ratio by this to get the per-eV ratio
        for conv, g, refval in (("A", gA, ref["A"]), ("B", gB, ref["B"])):
            reg = r[f"{conv}_regions"]
            raw = dict(cell_mean=float(g.mean() / refval), raised=reg["raised_terrace"] / refval, lower=reg["lower_terrace"] / refval,
                       lower_centre=reg["lower_terrace_centre_pm3A"] / refval, raised_centre=reg["raised_terrace_centre_pm3A"] / refval,
                       peaks=[p["value"] / refval for p in reg["peaks"]])
            per_eV = {k: (v * f_norm if not isinstance(v, list) else [x * f_norm for x in v]) for k, v in raw.items()}
            rat[conv] = dict(raw=raw, per_eV=per_eV)
            for nm, rr in (("raw", raw), ("per eV", per_eV)):
                print(f"{label} / [{fl}] ({conv}, {nm:6s}): cell mean {rr['cell_mean']:.3f}  raised {rr['raised']:.3f}  lower {rr['lower']:.3f}  "
                      f"lower centre {rr['lower_centre']:.3f}  peaks " + ", ".join("%.3f" % v for v in rr["peaks"]))
        rat["induced_sigma_raw"] = float(r["induced_sigma_e_per_A2"] / ref["sigma"])
        rat["induced_sigma_per_eV"] = float(rat["induced_sigma_raw"] * f_norm)
        rat["dmu_ratio_step_over_flat"] = float(1.0 / f_norm)
        print(f"{label} / [{fl}]: induced sigma ratio raw {rat['induced_sigma_raw']:.3f}, per eV {rat['induced_sigma_per_eV']:.3f}   "
              f"(dmu step {r['dmu_e_actual_eV']:+.4f} vs flat {ref['dmu']:+.4f} eV)")
        out[label].setdefault("ratios_to_flat", {})[fl] = rat

if step_recs:
    fig, axes = plt.subplots(2, len(step_recs), figsize=(6 * len(step_recs), 8.2), dpi=200, sharey="row",
                             gridspec_kw=dict(width_ratios=[r[1]["L_perp"] for r in step_recs]))
    axes = np.atleast_2d(axes) if len(step_recs) > 1 else np.array([[axes[0]], [axes[1]]])
    for j, (label, A, gA, gB, color) in enumerate(step_recs):
        s = A["u"] * A["L_perp"]; lo, hi = A["upper_u"]
        for row, g, conv, title in ((0, gA, "A", f"(A) absolute, z < {Z_UP:g} A"), (1, gB, "B", f"(B) boundary-relative, {DEPTH:g} A above local S=0.5")):
            ax = axes[row, j]
            ax.plot(s, g * 1e3, "-", lw=1.4, color=color, label=f"{label}")
            for ue in (0.0, 0.5):
                ax.axvline(ue * A["L_perp"], color="#333", lw=0.8, ls="--")
            ax.axvspan(lo * A["L_perp"], hi * A["L_perp"], color="#d9a441", alpha=0.18, label="raised terrace (atoms)")
            for fl, ref in flat_refs.items():
                ax.axhline(ref[conv] * 1e3, color=ref["color"], ls="--", lw=1.2, label=f"{fl.split(' (')[0]} cell mean")
            ax.grid(alpha=0.25); ax.legend(fontsize=7.5, loc="lower right")
            ax.set_title(f"{label} ({A['L_perp']/2:.1f} A/side) -- {title}", fontsize=10)
            if row == 1:
                ax.set_xlabel("across-step position s (A, perpendicular to the edge)")
    axes[0, 0].set_ylabel("dGamma_- per column (1e-3 ions/A^2)"); axes[1, 0].set_ylabel("dGamma_- per column (1e-3 ions/A^2)")
    plt.suptitle("Anion surface-excess response to the -0.2 eV mu_e step, resolved across the step (ion number per projected area)", y=0.995)
    plt.tight_layout()
    plt.savefig("report_assets/batch1/step_width_dGamma_profile.png", dpi=200)
json.dump(out, open("step_width_excess.json", "w"), indent=1)
print("wrote report_assets/batch1/step_width_dGamma_profile.png and step_width_excess.json")
