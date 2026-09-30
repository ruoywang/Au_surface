#!/usr/bin/env python3
"""Morphology -> charging response, from the converged CP records alone (no 3D field post-processing).

For every distinct geometry, using the ACTUAL converged electron count and electron chemical potential:

    q_e = N_e - N_e^0          (N_e^0 = 11 Z_Au per atom; q_e > 0 means excess electrons)
    Q   = -e q_e               (the metal's net charge)
    sigma(U) = Q / A_proj      A_proj = |a1 x a2| of the cell
    U   = mu_0 - mu_e(actual)  the project's internal relative potential, mu_0 = -4.9071 eV
                               (NOT vs RHE, NOT a per-structure PZC; U > 0 = the "positive" direction,
                                matching the campaign's dU labels)

Three quantities are reported per geometry:
  * sigma(U) at each sampled point                    -> how much net charge this morphology carries at a given U
  * the SECANT slope dsigma/dU between adjacent points -> how readily it keeps charging over that interval
  * the U where sigma = 0                              -> its potential of zero charge relative to the common reference,
                                                          reported ONLY when the sampled points bracket sigma = 0;
                                                          otherwise the distance to the nearest sampled point is given
                                                          and the row is marked "not bracketed" (no extrapolated PZC).

SCOPE. Every cross-structure number here uses only the three potentials every geometry has (-5.1071 / -4.9071 /
-4.7071). The +-0.5 V extension is still running and covers some geometries and not others; folding it in as it
arrives would compare a geometry sampled over +-0.5 V against one sampled over +-0.2 V and would move the published
numbers every time a job finished. The extension points are kept on each geometry as `all_points` for plotting.
No curvature fit is attempted on three points, and no PZC is extrapolated outside the sampled range.

It also evaluates the grand-potential comparison for pairs of geometries with the SAME composition and cell. With
dOmega/dmu = -N_e and U = mu_0 - mu_e, the potential-induced change across a symmetric window is

    D = [Omega_A - Omega_B]_{U=+w} - [Omega_A - Omega_B]_{U=-w} = + int_{mu_0-w}^{mu_0+w} [N_A(mu) - N_B(mu)] dmu

(the sign follows because U = +w is the LOWER mu). D < 0 means A is relatively stabilised as U moves positive.
This is the potential-INDUCED change only. It does NOT rank stability at the reference potential -- that needs the
relative Omega there, which this dataset has not fixed -- so it cannot say whether the potential re-orders two
morphologies. It is never applied across different Au counts (that needs a reservoir term).

Usage (from Au_Cl/):  scripts/pyrun.sh scripts/analysis_charging.py
Outputs: analysis/charging/{charging.json, charging_table.md, sigma_U.png, capacitance.png, pzc.png, domega.png}
"""
import collections
import itertools
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
OUT = f"{ROOT}/analysis/charging"
MU0 = -4.9071
E_PER_A2_TO_UC_PER_CM2 = 1.602176634e-19 / 1e-16 * 1e6      # 1 e/A^2 -> uC/cm^2
FAMILY_ORDER = ["flat Au(111)", "point defect", "reconstruction-related", "strip step", "vicinal step face",
                "kink / edge rearrangement", "single-layer island", "single-layer pit", "composite"]
FAMILY_COLOR = dict(zip(FAMILY_ORDER, ["#4c78a8", "#f58518", "#54a24b", "#b279a2", "#e45756",
                                       "#72b7b2", "#eeca3b", "#9d755d", "#bab0ac"]))


BASE_MU = {-5.1071, -4.9071, -4.7071}      # the three potentials EVERY geometry has; the +-0.5 V extension is partial
BASE_LABEL = "+-0.2 V subset (complete for every geometry)"


def load_geometries():
    """Every geometry, with its points split into the complete base set and the still-incomplete extension.

    Cross-structure numbers (PZC, secant capacitance, the scale comparison, the grand-potential pairs) are computed
    on the BASE set only. Mixing in the +-0.5 V states as they arrive would silently compare a geometry sampled over
    +-0.5 V with one sampled over +-0.2 V, and would change the published numbers every time a job finished."""
    S = json.load(open(f"{ROOT}/dataset_v1/states.json"))["states"]
    g = collections.defaultdict(lambda: dict(points=[]))
    for s in S.values():
        es = s["electronic_state"]
        cell = np.array(s["structure"]["cell_A"])
        d = g[s["geometry_id"]]
        d.update(geometry_id=s["geometry_id"], structure_id=s["structure_id"], family=s["family"], tier=s["tier"],
                 config=("relaxed" if s["config"] in ("relax", "relaxed") else s["config"]), n_atoms=s["structure"]["n_atoms"],
                 A_proj=float(np.linalg.norm(np.cross(cell[0], cell[1]))), cell=cell.tolist(),
                 N_neutral=es["N_neutral"])
        d["points"].append(dict(state_id=s["state_id"], TARGETMU=es["TARGETMU_eV"], mu_e=es["mu_e_actual_eV"],
                                N_e=es["N_e_final"], q_e=es["delta_N_e"], campaign=s["campaign"]))
    for d in g.values():
        A = d["A_proj"]
        for p in d["points"]:
            p["U"] = MU0 - p["mu_e"]                                   # internal relative potential, volts
            p["sigma_e_per_A2"] = -p["q_e"] / A
            p["sigma_uC_per_cm2"] = -p["q_e"] / A * E_PER_A2_TO_UC_PER_CM2
            p["base"] = round(p["TARGETMU"], 4) in BASE_MU
        d["points"].sort(key=lambda p: p["U"])
        d["all_points"] = d["points"]
        d["points"] = [p for p in d["points"] if p["base"]]            # everything below uses the base set only
        d["n_extension_points"] = len(d["all_points"]) - len(d["points"])
    return dict(g)


def secants(d):
    out = []
    for a, b in zip(d["points"], d["points"][1:]):
        dU = b["U"] - a["U"]
        if abs(dU) < 1e-6: continue
        out.append(dict(U_lo=a["U"], U_hi=b["U"], U_mid=0.5 * (a["U"] + b["U"]),
                        C_uF_per_cm2=(b["sigma_uC_per_cm2"] - a["sigma_uC_per_cm2"]) / dU,
                        C_e_per_V_per_A2=(b["sigma_e_per_A2"] - a["sigma_e_per_A2"]) / dU,
                        dN_e_per_V=(b["N_e"] - a["N_e"]) / (a["mu_e"] - b["mu_e"]) if abs(a["mu_e"] - b["mu_e"]) > 1e-9 else None,
                        from_states=[a["state_id"], b["state_id"]]))
    return out


def pzc(d):
    """U where sigma = 0, by linear interpolation, ONLY if the sampled points bracket it."""
    P = d["points"]; s = [p["sigma_e_per_A2"] for p in P]; u = [p["U"] for p in P]
    for i in range(len(P) - 1):
        if s[i] == 0.0: return dict(bracketed=True, U_pzc=u[i], note="a sampled point sits exactly at sigma = 0")
        if s[i] * s[i + 1] < 0:
            f = -s[i] / (s[i + 1] - s[i])
            return dict(bracketed=True, U_pzc=u[i] + f * (u[i + 1] - u[i]),
                        note=f"linear interpolation between {P[i]['state_id']} and {P[i+1]['state_id']}")
    j = int(np.argmin(np.abs(s)))
    side = "all sampled points carry the same sign of sigma"
    return dict(bracketed=False, U_pzc=None, nearest_U=u[j], nearest_sigma_uC_per_cm2=P[j]["sigma_uC_per_cm2"],
                note=f"{side}; NOT extrapolated (nearest point {P[j]['state_id']})")


U_WINDOW = 0.2      # the window every geometry has; the +-0.5 V extension will allow U_WINDOW = 0.5 later


def delta_omega_pairs(G, U_win=U_WINDOW):
    """Potential-induced relative stabilisation for same-composition, same-cell geometry pairs, on a FIXED U window
    so every pair is comparable (mixing different integration ranges would not be)."""
    mu_hi, mu_lo = MU0 + U_win, MU0 - U_win                      # mu = mu0 - U, so U=+win is the LOWER mu
    by = collections.defaultdict(list)
    for gid, d in G.items():
        if min(p["mu_e"] for p in d["points"]) > mu_lo + 1e-3 or max(p["mu_e"] for p in d["points"]) < mu_hi - 1e-3:
            continue                                             # does not cover the window: excluded, never extrapolated
        by[(d["n_atoms"], tuple(np.round(np.array(d["cell"]).ravel(), 3)))].append(gid)
    grid = np.linspace(mu_lo, mu_hi, 201)
    pairs = []
    for key, gids in by.items():
        if len(gids) < 2: continue
        for a, b in itertools.combinations(sorted(gids), 2):
            A, B = G[a], G[b]
            # np.interp needs an INCREASING x: points are stored sorted by U, i.e. by DECREASING mu_e, so reverse.
            pa = sorted(A["points"], key=lambda p: p["mu_e"]); pb = sorted(B["points"], key=lambda p: p["mu_e"])
            NA = np.interp(grid, [p["mu_e"] for p in pa], [p["N_e"] for p in pa])
            NB = np.interp(grid, [p["mu_e"] for p in pb], [p["N_e"] for p in pb])
            # SIGN. dOmega/dmu = -N, and U = mu0 - mu, so U = +win is the LOWER mu (mu_lo) and U = -win the higher.
            #   D = [Om_A-Om_B](U=+win) - [Om_A-Om_B](U=-win)
            #     = [Om_A-Om_B](mu_lo) - [Om_A-Om_B](mu_hi)
            #     = -int_{mu_hi}^{mu_lo}(N_A-N_B) dmu  =  +int_{mu_lo}^{mu_hi}(N_A-N_B) dmu
            # grid runs mu_lo -> mu_hi, so the integral is taken as-is. An earlier version negated it as well and
            # therefore reported the opposite sign: with N_A-N_B = 1 over a +-0.2 V window it gave -0.4 eV
            # instead of +0.4 eV, i.e. it named the wrong geometry as the one the potential stabilises.
            dOmega = np.trapz(NA - NB, grid)                     # eV, change in (Omega_A - Omega_B) across the window
            pairs.append(dict(a=a, b=b, n_atoms=A["n_atoms"], family=A["family"], U_window_V=U_win,
                              same_structure=A["structure_id"] == B["structure_id"],
                              d_relative_Omega_eV=float(dOmega), mean_dN_e=float(np.mean(NA - NB)),
                              dU_pzc_mV=(1000 * (A["pzc"]["U_pzc"] - B["pzc"]["U_pzc"])
                                         if A["pzc"]["bracketed"] and B["pzc"]["bracketed"] else None),
                              note="(Omega_A - Omega_B) at U=+win minus the same at U=-win; negative = A is relatively "
                                   "stabilised as U becomes more positive. Potential-INDUCED change only."))
    pairs.sort(key=lambda p: -abs(p["d_relative_Omega_eV"]))
    return pairs


CELL_REFERENCE = {                 # the flat / undefected member of each cell, when the cell has one
    "T-4x4__ideal": "flat (111) terrace, 4x4 cell",
    "Flat-16x1__ideal": "flat (111) terrace, 16x1 cell",
    "Flat-8x2__ideal": "flat (111) terrace, 8x2 cell",
}


def cell_groups(G):
    """Group geometries by their EXACT cell. Within a group the cell shape, k-mesh and grid are identical, so
    differences in U_pzc and C are morphology, not the numerical offset the audit records between cells
    (K.12: two 'flat' cells differ by 72 meV in neutral mu_e purely from cell shape / k-points / PREC)."""
    by = collections.defaultdict(list)
    for gid, d in G.items(): by[tuple(np.round(np.array(d["cell"]).ravel(), 2))].append(gid)
    groups = []
    for key, gids in by.items():
        ds = [G[g] for g in gids]
        ref = next((g for g in gids if f"{G[g]['structure_id']}__{G[g]['config']}" in CELL_REFERENCE), None)
        C = [s["C_uF_per_cm2"] for d in ds for s in d["secants"]] or [float("nan")]
        Z = [d["pzc"]["U_pzc"] for d in ds if d["pzc"]["bracketed"]]
        groups.append(dict(cell_key=str(key[:3]), A_proj=ds[0]["A_proj"], n=len(ds),
                           structures=sorted({d["structure_id"] for d in ds}), geometry_ids=sorted(gids),
                           reference=ref, reference_note=CELL_REFERENCE.get(
                               f"{G[ref]['structure_id']}__{G[ref]['config']}" if ref else "", None),
                           C_median=float(np.median(C)), C_min=float(min(C)), C_max=float(max(C)),
                           U_pzc_median=float(np.median(Z)) if Z else None,
                           U_pzc_min=float(min(Z)) if Z else None, U_pzc_max=float(max(Z)) if Z else None,
                           U_pzc_spread_mV=float(1000 * (max(Z) - min(Z))) if len(Z) > 1 else None))
        for d in ds:
            d["cell_group"] = groups[-1]["cell_key"]
            if ref and d["pzc"]["bracketed"] and G[ref]["pzc"]["bracketed"]:
                d["dU_pzc_vs_cell_ref_mV"] = 1000 * (d["pzc"]["U_pzc"] - G[ref]["pzc"]["U_pzc"])
                d["dC_vs_cell_ref_pct"] = 100 * (np.median([s["C_uF_per_cm2"] for s in d["secants"]]) /
                                                 np.median([s["C_uF_per_cm2"] for s in G[ref]["secants"]]) - 1)
    groups.sort(key=lambda g: -g["n"])
    return groups


def decompose(G, groups):
    """Compare the SCALE of the two effects. Not a variance decomposition, and not a share of the variation.

    With sigma_i = C_i (U - Z_i), expanding about a reference (C_0, Z_0) gives
        d sigma_i = -C_0 dZ_i + (U - Z_0) dC_i - dC_i dZ_i,
    i.e. a PZC term, a capacitance term AND a cross term, and across geometries dC and dZ are correlated. What is
    computed here is only the size of the first two terms evaluated on the max-min SPANS:
        A = <C> (max Z - min Z)        B = (max C - min C) |U - <Z>|
    A/B says which effect operates on the larger scale in this sample. It must NOT be read as "A explains
    A/(A+B) of the charge difference" -- that would need the cross term and the correlations, which this does not
    touch. Reported as a scale ratio for exactly that reason."""
    out = {}
    for tag, sel in [("all geometries", list(G.values()))] + \
                    [(f"within cell A={g['A_proj']:.0f} A^2 ({', '.join(g['structures'][:3])}{'...' if len(g['structures'])>3 else ''})",
                      [G[x] for x in g["geometry_ids"]]) for g in groups if g["n"] >= 4 and g["U_pzc_spread_mV"]]:
        sel = [d for d in sel if d["secants"] and d["pzc"]["bracketed"]]      # need >= 2 sampled points for either
        Z = [d["pzc"]["U_pzc"] for d in sel]
        C = [float(np.median([s["C_uF_per_cm2"] for s in d["secants"]])) for d in sel]
        if len(Z) < 2: continue
        Cm, Zm = float(np.median(C)), float(np.median(Z))
        for U in (0.2, -0.2):
            out.setdefault(tag, {})[f"U={U:+.1f}V"] = dict(
                scale_from_pzc_span_uC_per_cm2=Cm * (max(Z) - min(Z)),
                scale_from_C_span_uC_per_cm2=(max(C) - min(C)) * abs(U - Zm),
                scale_ratio=(Cm * (max(Z) - min(Z))) / max(1e-9, (max(C) - min(C)) * abs(U - Zm)))
        out[tag].update(n=len(sel), U_pzc_spread_mV=1000 * (max(Z) - min(Z)), C_spread_pct=100 * (max(C) - min(C)) / Cm,
                        C_median=Cm, U_pzc_median=Zm)
    return out


G = load_geometries()
for gid, d in G.items():
    d["secants"] = secants(d); d["pzc"] = pzc(d); d["n_points"] = len(d["points"])
groups = cell_groups(G)
decomp = decompose(G, groups)
pairs = delta_omega_pairs(G)
os.makedirs(OUT, exist_ok=True)
json.dump(dict(created=__import__("time").strftime("%Y-%m-%d %H:%M"), mu_reference_eV=MU0,
               convention="U = mu0 - mu_e(actual) in volts; sigma = -e (N_e - N_e^0)/A_proj",
               n_geometries=len(G), geometries=G, cell_groups=groups, decomposition=decomp,
               same_composition_pairs=pairs), open(f"{OUT}/charging.json", "w"), indent=1)

# ---------------------------------------------------------------- table
rows = sorted(G.values(), key=lambda d: (FAMILY_ORDER.index(d["family"]) if d["family"] in FAMILY_ORDER else 99,
                                         d["structure_id"], d["config"]))
L = ["# Morphology -> charging response (from the converged CP records)", "",
     f"Generated {__import__('time').strftime('%Y-%m-%d %H:%M')} from `dataset_v1/states.json`. "
     f"{len(G)} distinct geometries. U = mu0 - mu_e(actual), mu0 = {MU0} eV (internal reference, not vs RHE, not a "
     "per-structure PZC). sigma = -e (N_e - N_e^0) / A_proj with N_e^0 = 11 per Au atom. Every value uses the ACTUAL "
     "converged mu_e and N_e, never TARGETMU.", "",
     "`C_sec` is the secant slope between adjacent sampled points, not a differential capacitance fit. A PZC is given "
     "only where the sampled points bracket sigma = 0; otherwise the row says how far the nearest point is.", "",
     "| family | structure | config | N | A_proj (A^2) | pts | sigma at U (uC/cm^2) | C_sec (uF/cm^2) | U(sigma=0) |",
     "|---|---|---|---|---|---|---|---|---|"]
for d in rows:
    sig = " ; ".join(f"{p['U']:+.3f}V: {p['sigma_uC_per_cm2']:+.2f}" for p in d["points"])
    cap = " ; ".join(f"{s['C_uF_per_cm2']:.1f}@{s['U_mid']:+.2f}" for s in d["secants"])
    z = d["pzc"]
    zs = f"{z['U_pzc']:+.4f} V" if z["bracketed"] else f"not bracketed (nearest {z['nearest_U']:+.3f} V, {z['nearest_sigma_uC_per_cm2']:+.2f})"
    L.append(f"| {d['family']} | {d['structure_id']} | {d['config']} | {d['n_atoms']} | {d['A_proj']:.1f} | "
             f"{d['n_points']} | {sig} | {cap} | {zs} |")

nb = [d for d in rows if not d["pzc"]["bracketed"]]
L += ["", f"- PZC bracketed by the sampled range: {len(rows)-len(nb)} of {len(rows)} geometries; not bracketed: {len(nb)}.", "",
      "## Within one cell: morphology only", "",
      "Cross-cell PZC comparisons carry a numerical offset: audit K.12 records that two *flat* cells differ by 72 meV in "
      "neutral mu_e purely from cell shape, k-mesh and PREC. Inside one cell all of that is common, so the columns below "
      "are morphology. `ref` is the flat member of that cell where one exists.", "",
      "| A_proj (A^2) | n | structures | reference | secant C med (range) | U_pzc med (range, mV) | spread (mV) |",
      "|---|---|---|---|---|---|---|"]
for g in groups:
    st = ", ".join(g["structures"])[:60]
    z = (f"{1000*g['U_pzc_median']:+.0f} ({1000*g['U_pzc_min']:+.0f}..{1000*g['U_pzc_max']:+.0f})"
         if g["U_pzc_median"] is not None else "-")
    sp = f"{g['U_pzc_spread_mV']:.0f}" if g["U_pzc_spread_mV"] else "-"
    L.append(f"| {g['A_proj']:.1f} | {g['n']} | {st} | {g['reference'] or '-'} | "
             f"{g['C_median']:.2f} ({g['C_min']:.2f}-{g['C_max']:.2f}) | {z} | {sp} |")
L += ["", "## Which matters more at a given potential: the PZC shift or the capacitance?", "",
      "Writing sigma(U) = C (U - U_pzc), the spread of sigma across geometries at a fixed U splits into "
      "`<C> x spread(U_pzc)` and `spread(C) x |U - <U_pzc>|`.", "",
      "| set | n | spread U_pzc (mV) | spread C (%) | sigma spread from PZC | from C | ratio |", "|---|---|---|---|---|---|---|"]
for tag, v in decomp.items():
    w = v["U=+0.2V"]
    L.append(f"| {tag} | {v['n']} | {v['U_pzc_spread_mV']:.0f} | {v['C_spread_pct']:.1f} | "
             f"{w['scale_from_pzc_span_uC_per_cm2']:.2f} uC/cm2 | {w['scale_from_C_span_uC_per_cm2']:.2f} uC/cm2 | "
             f"{w['scale_ratio']:.1f}x |")
L += ["", "(evaluated at U = +0.2 V; the ratio grows as U approaches the median PZC and shrinks far from it.)", ""]
cross = [p for p in pairs if not p["same_structure"]]
same = [p for p in pairs if p["same_structure"]]
L += [f"## Potential-induced relative stabilisation, same composition and cell (fixed U window $\\pm${U_WINDOW} V)", "",
      "D = +int_{mu0-w}^{mu0+w} [N_A(mu) - N_B(mu)] dmu over the SAME window for every pair, so the values are "
      "comparable. The sign follows from dOmega/dmu = -N and U = mu0 - mu, which puts U = +w at the LOWER mu; an "
      "earlier version negated the integral as well and so named the opposite geometry. D < 0 means A is relatively "
      "stabilised as U becomes more positive. This is the potential-INDUCED change ONLY: without the relative Omega "
      "at the reference potential it cannot say whether the potential re-orders the two morphologies, and it is never "
      "taken across different Au counts. N_A - N_B is nearly constant over this window, so D tracks the PZC offset, "
      "shown alongside.", "",
      f"{len(pairs)} pairs cover the window; {len(cross)} are between DIFFERENT structures. Those, largest first:", "",
      "| A | B | N | d(Omega_A-Omega_B) (eV) | mean N_A-N_B (e) | dU_pzc (mV) |", "|---|---|---|---|---|---|"]
for p in cross[:40]:
    z = f"{p['dU_pzc_mV']:+.1f}" if p["dU_pzc_mV"] is not None else "-"
    L.append(f"| {p['a']} | {p['b']} | {p['n_atoms']} | {p['d_relative_Omega_eV']:+.4f} | {p['mean_dN_e']:+.4f} | {z} |")
if same:
    v = [abs(p["d_relative_Omega_eV"]) for p in same]
    L += ["", f"- Same-structure pairs (ideal / relaxed / perturbed / collective / path images of one parent, "
          f"{len(same)} pairs): |dOmega| median {np.median(v):.4f} eV, max {max(v):.4f} eV. That is the scale on which "
          "sampling configurations of ONE morphology already differ, i.e. the bar a cross-morphology difference must "
          "clear to be meaningful.", ""]
open(f"{OUT}/charging_table.md", "w").write("\n".join(L) + "\n")

# ---------------------------------------------------------------- figures
def famfig(fname, ykey, ylabel, title):
    fig, ax = plt.subplots(figsize=(7.2, 4.6), dpi=160)
    for fam in FAMILY_ORDER:
        ds = [d for d in G.values() if d["family"] == fam]
        if not ds: continue
        first = True
        for d in ds:
            x = [p["U"] for p in d["points"]]; y = [p[ykey] for p in d["points"]]
            ax.plot(x, y, "-o", ms=2.6, lw=1.0, color=FAMILY_COLOR[fam], alpha=0.75,
                    label=fam if first else None); first = False
    ax.axhline(0, color="#999", lw=0.8, zorder=0); ax.axvline(0, color="#999", lw=0.8, zorder=0)
    ax.set_xlabel("U = $\\mu_0-\\mu_e$  (V, internal reference)"); ax.set_ylabel(ylabel); ax.set_title(title, fontsize=10)
    ax.legend(fontsize=6.5, ncol=2, frameon=False); fig.tight_layout(); fig.savefig(f"{OUT}/{fname}"); plt.close(fig)

famfig("sigma_U.png", "sigma_uC_per_cm2", "$\\sigma$  ($\\mu$C/cm$^2$)", "Surface charge density vs internal potential, all geometries")

fig, ax = plt.subplots(figsize=(7.2, 4.6), dpi=160)
for fam in FAMILY_ORDER:
    pts = [(s["U_mid"], s["C_uF_per_cm2"]) for d in G.values() if d["family"] == fam for s in d["secants"]]
    if pts: ax.scatter(*zip(*pts), s=14, color=FAMILY_COLOR[fam], alpha=0.8, label=fam, edgecolors="none")
ax.set_xlabel("interval midpoint U (V)"); ax.set_ylabel("secant $\\Delta\\sigma/\\Delta U$  ($\\mu$F/cm$^2$)")
ax.set_title("Secant capacitance per sampled interval", fontsize=10); ax.legend(fontsize=6.5, ncol=2, frameon=False)
fig.tight_layout(); fig.savefig(f"{OUT}/capacitance.png"); plt.close(fig)

fig, ax = plt.subplots(figsize=(7.2, 4.2), dpi=160)
br = [(d["pzc"]["U_pzc"], d["family"]) for d in G.values() if d["pzc"]["bracketed"]]
for fam in FAMILY_ORDER:
    v = [b[0] for b in br if b[1] == fam]
    if v: ax.scatter(v, [FAMILY_ORDER.index(fam)] * len(v), s=18, color=FAMILY_COLOR[fam], alpha=0.85, edgecolors="none")
ax.set_yticks(range(len(FAMILY_ORDER))); ax.set_yticklabels(FAMILY_ORDER, fontsize=7)
ax.set_xlabel("U where $\\sigma=0$  (V, internal reference)"); ax.axvline(0, color="#999", lw=0.8)
ax.set_title(f"Potential of zero charge, {len(br)} geometries where it is bracketed by the sampled points", fontsize=10)
fig.tight_layout(); fig.savefig(f"{OUT}/pzc.png"); plt.close(fig)

def short(gid):
    p = gid.split("__"); return f"{p[0]} {p[1]}"


if cross:
    sel = cross[:22]
    fig, ax = plt.subplots(figsize=(8.4, max(2.4, 0.3 * len(sel))), dpi=160)
    ax.barh(range(len(sel)), [p["d_relative_Omega_eV"] for p in sel], color="#4c78a8")
    ax.set_yticks(range(len(sel))); ax.set_yticklabels([f"{short(p['a'])}  vs  {short(p['b'])}" for p in sel], fontsize=6)
    ax.invert_yaxis(); ax.axvline(0, color="#333", lw=0.8)
    ax.set_xlabel(f"$\\Delta(\\Omega_A-\\Omega_B)$ across U = $\\pm${U_WINDOW} V  (eV)")
    ax.set_title("Potential-induced relative stabilisation, same composition and cell, different morphology", fontsize=9)
    fig.tight_layout(); fig.savefig(f"{OUT}/domega.png"); plt.close(fig)

print(f"{len(G)} geometries -> {OUT}/")
print(f"  PZC bracketed: {sum(1 for d in G.values() if d['pzc']['bracketed'])}/{len(G)}")
print(f"  same-composition pairs: {len(pairs)}")
caps = [s["C_uF_per_cm2"] for d in G.values() for s in d["secants"]]
print(f"  secant capacitance over all intervals: median {np.median(caps):.1f}, range {min(caps):.1f}..{max(caps):.1f} uF/cm^2")
