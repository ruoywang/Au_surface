#!/usr/bin/env python3
"""Turn the per-column reductions into the two spatial answers, plus the figures and the web data pack.

(1) WHERE the anions go, resolved by morphology. For every reduced state, region-resolved
        K_Omega   = int n_- / (n_b int S_ion)            concentration relative to bulk inside that pocket
        Gamma_-   = (1/A) int [n_- - n_b S_ion]          extra anions per unit projected area
    in four coordination classes (CN <= 6 kink/adatom, 7-8 edge/rim, 9 terrace, >= 10 sub-surface / foot), under BOTH
    integration conventions, so no conclusion rests on the window choice.

(2) HOW MUCH of the metal's lateral structure survives into the ion distribution. For each geometry, between the two
    end potentials, the per-column changes
        d n_e(x,y)  (metal)  ->  d rho_b(x,y)  (dielectric bound charge)  ->  d Gamma_-(x,y)  (anions)
    are compared for: peak position (the lateral OFFSET between the metal response and the ion response, from the 2D
    cross-correlation, minimum-image), lateral CONTRAST (peak-to-peak over the mean), and the contrast RATIO between
    the two sides, which is what "how much is screened" means quantitatively.

Usage (from Au_Cl/):  scripts/pyrun.sh scripts/analysis_spatial_report.py [--figures]
Outputs: analysis/spatial/{regions.json, transmission.json}, analysis/maps/<geometry>.png, analysis/web/data.json
"""
import argparse
import collections
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "/anvil/scratch/x-rywang/Au_Cl/scripts")
from analysis_spatial import N_BULK, Z_UP, DEPTH  # noqa: E402

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
PC = f"{ROOT}/analysis/spatial/percolumn"
MAPS = f"{ROOT}/analysis/maps"
WEB = f"{ROOT}/analysis/web"
MU0 = -4.9071
CLASSES = [("kink/adatom", lambda l: l <= 6), ("edge/rim", lambda l: (l >= 7) & (l <= 8)),
           ("terrace", lambda l: l == 9), ("sub-surface/foot", lambda l: l >= 10)]


def load(sid):
    z = np.load(f"{PC}/{sid}.npz")
    r = {k: z[k] for k in z.files}
    r["A_proj"] = float(np.linalg.norm(np.cross(r["cell"][0], r["cell"][1])))   # scalars live in the npz's cell array
    return r


def regions_of(r):
    out = {}
    lab = r["cn_label"]
    for name, f in CLASSES:
        m = f(lab)
        if not m.any(): continue
        out[name] = dict(
            area_fraction=float(m.mean()), zb_mean=float(r["zb"][m].mean()),
            K_rel=float(r["nminus_rel"][m].sum() / (N_BULK * r["sion_rel"][m].sum())),
            K_abs=float(r["nminus_abs"][m].sum() / (N_BULK * r["sion_abs"][m].sum())),
            Gamma_rel=float(r["gamma_rel"][m].mean()), Gamma_abs=float(r["gamma_abs"][m].mean()),
            Gamma_share_rel=float(r["gamma_rel"][m].sum() / r["gamma_rel"].sum()) if r["gamma_rel"].sum() else None)
    out["all"] = dict(area_fraction=1.0, zb_mean=float(r["zb"].mean()),
                      K_rel=float(r["nminus_rel"].sum() / (N_BULK * r["sion_rel"].sum())),
                      K_abs=float(r["nminus_abs"].sum() / (N_BULK * r["sion_abs"].sum())),
                      Gamma_rel=float(r["gamma_rel"].mean()), Gamma_abs=float(r["gamma_abs"].mean()),
                      Gamma_share_rel=1.0)
    return out


def contrast(a):
    """peak-to-peak of the lateral pattern, normalised by its mean magnitude (dimensionless)."""
    m = float(np.abs(a).mean())
    return float(a.max() - a.min()), (float((a.max() - a.min()) / m) if m > 0 else float("nan"))


MIN_CONTRAST = 0.4      # below this the metal map has no lateral structure worth locating (a flat terrace is ~0.75
                        # only because of the atomic corrugation, which the ion side does not resolve at all)


def dominant_mode(B):
    """(h, k) of the strongest non-uniform Fourier component of the metal pattern, as minimum-image integers."""
    F = np.abs(np.fft.fft2(B - B.mean())) ** 2
    F[0, 0] = 0
    ngy, ngx = B.shape
    iy, ix = np.unravel_index(np.argmax(F), F.shape)
    h = ix - ngx * (ix > ngx // 2); k = iy - ngy * (iy > ngy // 2)
    return int(h), int(k), float(F.max() / F.sum())


def lateral_offset(ion, metal, cell):
    """Offset between the ION response and the METAL response, along the direction the metal actually modulates.

    A step or a stripe modulates in ONE direction; a 2D cross-correlation is then degenerate along the step and its
    peak position along that axis is noise. So: take the strongest Fourier mode (h,k) of the metal map, project both
    maps onto the scalar phase s = h f1 + k f2 (the only coordinate either pattern depends on), and read the shift
    from the phase difference of that mode. The real-space displacement is ds / |G|, G = h b1 + k b2 with b_i the
    reciprocal vectors (b_i . a_j = delta_ij), directed along G."""
    mpp, mc = contrast(metal)
    if mc < MIN_CONTRAST or not np.isfinite(mc):
        return dict(offset_A=None, reason=f"metal lateral contrast {mc:.3f} below {MIN_CONTRAST}: no pattern to locate")
    h, k, frac = dominant_mode(metal)
    ngy, ngx = metal.shape
    Fi = np.fft.fft2(ion - ion.mean())[k % ngy, h % ngx]
    Fm = np.fft.fft2(metal - metal.mean())[k % ngy, h % ngx]
    if abs(Fm) == 0 or abs(Fi) == 0:
        return dict(offset_A=None, reason="the ion map carries no amplitude in the metal's dominant mode")
    dphase = np.angle(Fi / Fm) / (2 * np.pi)                      # in units of the mode's period, in (-0.5, 0.5]
    a1, a2 = cell[0][:2], cell[1][:2]
    M = np.array([a1, a2]); B = np.linalg.inv(M).T                # rows b1, b2 with b_i . a_j = delta_ij
    G = h * B[0] + k * B[1]
    period = 1.0 / np.linalg.norm(G)
    d = dphase * period
    amp_ratio = float(abs(Fi) / abs(Fm))
    return dict(offset_A=float(abs(d)), signed_offset_A=float(d), period_A=float(period), mode=[h, k],
                mode_power_fraction=frac, direction=[float(x) for x in G / np.linalg.norm(G)],
                amplitude_ratio_ion_over_metal=amp_ratio,
                reason=None)


LAMBDA_D = 3.04     # A, Debye length of a 1:1 1 M aqueous electrolyte at 298 K (the model's bulk concentration)


def transfer_function(ion, metal, cell, nbins=14):
    """Empirical wavelength-resolved transmission |F_ion(k)| / |F_metal(k)|.

    The point of item 3: a metal charge pattern of lateral wavevector k reaches the ion-accessible region attenuated.
    Linearised Poisson-Boltzmann over a flat, uniform dielectric gives d psi_k(z) ~ exp(-sqrt(k^2 + kappa_D^2) z), so
    SHORT-wavelength (atomic) structure should be filtered out far more strongly than long-wavelength (defect-scale)
    structure. This measures that filter directly, without assuming it: every Fourier mode of the two per-column maps
    is compared, and the ratio is binned by |k|. It is a diagnostic of the full nonlinear, space-varying-cavity model,
    so departures from the flat-PB reference are the interesting part, not an error."""
    ngy, ngx = metal.shape
    Fm = np.fft.fft2(metal - metal.mean()); Fi = np.fft.fft2(ion - ion.mean())
    a1, a2 = cell[0][:2], cell[1][:2]
    B = np.linalg.inv(np.array([a1, a2])).T                       # b_i . a_j = delta_ij
    h = np.fft.fftfreq(ngx) * ngx; k = np.fft.fftfreq(ngy) * ngy
    H, K = np.meshgrid(h, k)
    Gx = H * B[0][0] + K * B[1][0]; Gy = H * B[0][1] + K * B[1][1]
    kmag = 2 * np.pi * np.hypot(Gx, Gy)                            # rad/A
    am, ai = np.abs(Fm), np.abs(Fi)
    keep = (kmag > 1e-6) & (am > 0.05 * am[kmag > 1e-6].max())     # only modes the metal actually carries
    if keep.sum() < 8: return None
    kk, ratio = kmag[keep], (ai[keep] / am[keep])
    edges = np.geomspace(max(kk.min(), 0.05), kk.max(), nbins + 1)
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (kk >= lo) & (kk < hi)
        if m.sum() < 2: continue
        out.append(dict(k=float(np.sqrt(lo * hi)), lam=float(2 * np.pi / np.sqrt(lo * hi)), n=int(m.sum()),
                        T=float(np.median(ratio[m])), T_lo=float(np.percentile(ratio[m], 25)),
                        T_hi=float(np.percentile(ratio[m], 75))))
    if len(out) < 4: return None
    # T compares two DIFFERENT quantities (electrons/A^2 against anions/A^2), so its absolute scale carries no
    # meaning -- only the wavelength dependence does. Normalise each curve by its longest-wavelength bin, which is
    # the closest thing to the local "flat capacitor" limit this cell can resolve.
    T_ref = max(o["T"] for o in out[:2]) if out else 1.0
    for o in out:
        o["Tn"] = o["T"] / T_ref if T_ref > 0 else float("nan")
    half = None
    for p, q in zip(out[:-1], out[1:]):                            # first crossing of Tn = 0.5, walking to shorter lambda
        if p["Tn"] >= 0.5 > q["Tn"]:
            f = (p["Tn"] - 0.5) / (p["Tn"] - q["Tn"])
            half = float(np.exp(np.log(p["lam"]) + f * (np.log(q["lam"]) - np.log(p["lam"])))); break
    kap = 1.0 / LAMBDA_D
    x = np.array([np.sqrt(o["k"] ** 2 + kap ** 2) - kap for o in out])
    y = np.log(np.array([max(o["Tn"], 1e-12) for o in out]))
    A = np.polyfit(x, y, 1)
    return dict(bins=out, lambda_half_A=half, d_eff_A=float(-A[0]), kappa_D=kap,
                r2=float(1 - np.var(y - np.polyval(A, x)) / max(np.var(y), 1e-30)))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--figures", action="store_true"); a = ap.parse_args()
    os.makedirs(MAPS, exist_ok=True); os.makedirs(WEB, exist_ok=True)
    S = json.load(open(f"{ROOT}/dataset_v1/states.json"))["states"]
    ch = json.load(open(f"{ROOT}/analysis/charging/charging.json"))
    have = {f[:-4] for f in os.listdir(PC)}
    reg, byg = {}, collections.defaultdict(dict)
    for sid in sorted(have):
        s = S.get(sid)
        if s is None: continue
        r = load(sid)
        U = MU0 - s["electronic_state"]["mu_e_actual_eV"]
        reg[sid] = dict(geometry_id=s["geometry_id"], structure_id=s["structure_id"], family=s["family"],
                        config=("relaxed" if s["config"] in ("relax", "relaxed") else s["config"]),
                        U=U, n_atoms=s["structure"]["n_atoms"], A_proj=float(r["A_proj"]),
                        cols_fitting_depth=float((r["zb"] + DEPTH <= Z_UP).mean()), regions=regions_of(r))
        byg[s["geometry_id"]][round(U, 3)] = sid
    json.dump(dict(convention=dict(Z_UP=Z_UP, DEPTH=DEPTH, n_bulk_A3=N_BULK), n_states=len(reg), states=reg),
              open(f"{ROOT}/analysis/spatial/regions.json", "w"), indent=1)
    print(f"region analysis: {len(reg)} states")

    # ---------------- transmission: metal -> bound charge -> anions, between the two end potentials
    trans = {}
    for gid, pts in byg.items():
        us = sorted(pts)
        lo, hi = us[0], us[-1]
        if hi - lo < 0.3: continue
        A, B = load(pts[hi]), load(pts[lo])
        if "ne_col" not in A or "ne_col" not in B: continue
        cell = A["cell"]
        d_ne = A["ne_col"] - B["ne_col"]                              # metal electrons per A^2 (more positive U -> fewer)
        d_gam = A["gamma_rel"] - B["gamma_rel"]                       # extra anions per A^2
        d_rhob = A["rhob_col"] - B["rhob_col"] if "rhob_col" in A else None
        d_phi = A["phi_plane"] - B["phi_plane"] if "phi_plane" in A else None
        e = {}
        for name, arr in [("metal_dne", d_ne), ("bound_charge_drhob", d_rhob), ("potential_dphi", d_phi),
                          ("anion_dGamma", d_gam)]:
            if arr is None: continue
            pp, c = contrast(arr)
            e[name] = dict(peak_to_peak=pp, contrast=c, mean=float(arr.mean()))
        off = lateral_offset(d_gam, d_ne, cell)
        tf = transfer_function(d_gam, d_ne, cell)
        s0 = S[pts[hi]]
        trans[gid] = dict(structure_id=s0["structure_id"], family=s0["family"],
                          config=("relaxed" if s0["config"] in ("relax", "relaxed") else s0["config"]),
                          U_lo=lo, U_hi=hi, n_atoms=s0["structure"]["n_atoms"], fields=e, offset=off, transfer=tf,
                          contrast_ratio_ion_over_metal=(e["anion_dGamma"]["contrast"] / e["metal_dne"]["contrast"]
                                                         if "metal_dne" in e and e["metal_dne"]["contrast"] else None))
    json.dump(dict(n_geometries=len(trans), geometries=trans), open(f"{ROOT}/analysis/spatial/transmission.json", "w"), indent=1)
    print(f"transmission analysis: {len(trans)} geometries")

    # ---------------- per-structure maps
    if a.figures:
        gal = json.load(open(f"{ROOT}/analysis/gallery/gallery.json"))
        n = 0
        for gid, pts in sorted(byg.items()):
            sid_struct = gid.split("__")[0]; cfg = gid.split("__")[1]
            if cfg not in ("ideal", "relaxed") or sid_struct not in gal: continue
            if gal[sid_struct]["geometry_id"] != gid: continue
            us = sorted(pts)
            if len(us) < 2: continue
            A = load(pts[us[-1]]); B = load(pts[us[0]])
            make_map(sid_struct, gid, A, B, us[-1], us[0]); n += 1
        print(f"maps: {n} figures -> {MAPS}/")

    pack(ch, reg, trans)


def make_map(struct, gid, A, B, U_hi, U_lo):
    cell = A["cell"]; ngy, ngx = A["gamma_rel"].shape
    fx, fy = np.meshgrid((np.arange(ngx) + 0.5) / ngx, (np.arange(ngy) + 0.5) / ngy)
    X = fx * cell[0][0] + fy * cell[1][0]
    Y = fx * cell[0][1] + fy * cell[1][1]
    panels = [("anion excess $\\Gamma_-$ at U = %+.2f V" % U_hi, A["gamma_rel"], "RdBu_r", None),
              ("change $\\Delta\\Gamma_-$, U %+.2f $\\to$ %+.2f V" % (U_lo, U_hi), A["gamma_rel"] - B["gamma_rel"], "RdBu_r", None),
              ("metal $\\Delta n_e$ over the same step", (A["ne_col"] - B["ne_col"]) if "ne_col" in A else None, "PuOr", None),
              ("surface coordination number", A["cn_label"].astype(float), "viridis", (5.5, 12.5))]
    panels = [p for p in panels if p[1] is not None]
    fig, axes = plt.subplots(1, len(panels), figsize=(3.5 * len(panels), 3.8), dpi=170)
    for ax, (t, D, cm, lim) in zip(np.atleast_1d(axes), panels):
        if lim: v0, v1 = lim
        else:
            m = np.abs(D - (0 if D.min() * D.max() < 0 else D.mean())).max()
            v0, v1 = (-m, m) if D.min() * D.max() < 0 else (D.min(), D.max())
        im = ax.pcolormesh(X, Y, D, cmap=cm, vmin=v0, vmax=v1, shading="nearest", rasterized=True)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(t, fontsize=7.6, color="#41474d", pad=3)
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02); cb.ax.tick_params(labelsize=6)
    fig.suptitle(f"{struct} · {gid.split('__')[1]}", fontsize=11, y=0.99, color="#14181c")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(f"{MAPS}/{struct}.png", facecolor="white"); plt.close(fig)


def pack(ch, reg, trans):
    """Compact JSON for the web page: only what the charts need."""
    G = ch["geometries"]
    geoms = {gid: dict(s=d["structure_id"], c=d["config"], f=d["family"], n=d["n_atoms"], A=round(d["A_proj"], 1),
                       pts=[[round(p["U"], 4), round(p["sigma_uC_per_cm2"], 4)] for p in d["points"]],
                       C=round(float(np.median([x["C_uF_per_cm2"] for x in d["secants"]])), 3) if d["secants"] else None,
                       z=(round(d["pzc"]["U_pzc"], 5) if d["pzc"]["bracketed"] else None),
                       dz=(round(d["dU_pzc_vs_cell_ref_mV"], 1) if "dU_pzc_vs_cell_ref_mV" in d else None),
                       dC=(round(d["dC_vs_cell_ref_pct"], 2) if "dC_vs_cell_ref_pct" in d else None))
             for gid, d in G.items()}
    gal = json.load(open(f"{ROOT}/analysis/gallery/gallery.json"))
    structs = {k: dict(n=v["n_atoms"], f=v["family"], A=round(v.get("A_proj", 0), 1), cn=v["cn_counts_dict"],
                       z=v.get("U_pzc_V"), C=v.get("C_median_uF_cm2"), dz=v.get("dU_pzc_vs_cell_ref_mV"),
                       dC=v.get("dC_vs_cell_ref_pct"), gid=v["geometry_id"], cfg=v["cfg"]) for k, v in gal.items()}
    # region summary per structure at the most positive U available
    rsum = {}
    for sid, d in reg.items():
        if d["config"] not in ("ideal", "relaxed"): continue
        k = d["structure_id"]
        if k not in structs or structs[k]["gid"] != d["geometry_id"]: continue
        cur = rsum.get(k)
        if cur is None or d["U"] > cur["U"]:
            rsum[k] = dict(U=round(d["U"], 4), regions={n: dict(a=round(v["area_fraction"], 4),
                                                                K=round(v["K_rel"], 5), Ka=round(v["K_abs"], 5),
                                                                G=v["Gamma_rel"], sh=v["Gamma_share_rel"])
                                                        for n, v in d["regions"].items()})
    tr = {g: dict(s=v["structure_id"], f=v["family"], cfg=v["config"], off=v["offset"].get("offset_A"),
                  period=v["offset"].get("period_A"), amp=v["offset"].get("amplitude_ratio_ion_over_metal"),
                  why=v["offset"].get("reason"), ratio=v["contrast_ratio_ion_over_metal"],
                  tf=(dict(d=round(v["transfer"]["d_eff_A"], 3), r2=round(v["transfer"]["r2"], 3),
                           half=(round(v["transfer"]["lambda_half_A"], 2) if v["transfer"]["lambda_half_A"] else None),
                           b=[[round(o["lam"], 3), round(o["Tn"], 6)] for o in v["transfer"]["bins"]])
                      if v.get("transfer") else None),
                  c={k: round(x["contrast"], 4) for k, x in v["fields"].items()}) for g, v in trans.items()}
    out = dict(mu0=MU0, decomposition=ch["decomposition"], cell_groups=ch["cell_groups"],
               pairs=[p for p in ch["same_composition_pairs"] if not p["same_structure"]][:60],
               same_structure_pair_scale=dict(
                   median=float(np.median([abs(p["d_relative_Omega_eV"]) for p in ch["same_composition_pairs"] if p["same_structure"]])),
                   max=float(max([abs(p["d_relative_Omega_eV"]) for p in ch["same_composition_pairs"] if p["same_structure"]]))),
               geometries=geoms, structures=structs, regions=rsum, transmission=tr,
               convention=dict(Z_UP=Z_UP, DEPTH=DEPTH))
    json.dump(out, open(f"{WEB}/data.json", "w"), separators=(",", ":"))
    print(f"web pack: {os.path.getsize(WEB+'/data.json')/1024:.0f} kB -> {WEB}/data.json")


if __name__ == "__main__":
    main()
