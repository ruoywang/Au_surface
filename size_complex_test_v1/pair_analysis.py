#!/usr/bin/env python3
"""Compare the finished members of each size pair: label completeness of each cell on its own, local force differences
on the truly matched atoms, potential and ion profiles in the same real-space region above the core, and measured cost.

A. each member alone: CP closed on target, 15 fields parsed value by value, charge closure, CONTCAR = POSCAR (fieldio);
   thermal forces are not a failure, nor is a different total energy.
B. forces: only atoms present in both cells whose 6 A neighbourhood (as parent atom ids) is identical in both -- the
   matched set C -- enter RMS dF = sqrt( sum_i |F_i(6x6) - F_i(8x8)|^2 / (3 |C|) ); also max |dF|, where it sits and
   its distance to the nearest 6x6 seam; a protected atom whose own neighbours lie outside the protected set is NOT
   assumed matched -- it is checked like every other atom.
C. fields: both cells aligned on the core; PHI, S_ion, n- and n+ (analysis_spatial.species, the validated convention
   and zero) averaged over a 4 A disc above the centre atom as functions of height; Gamma- = integral of n- over the
   window; the accessible boundary (first height with S_ion > 0.5) of each; sigma = -e (N_e - 11 N) / A as a
   supplement only (different defect densities; its inequality proves nothing by itself).
D. cost: atoms, k-points, FFT grid, NBANDS, electronic steps, CP rounds, SCF time (sum of LOOP real times), elapsed,
   non-SCF time (elapsed - SCF, includes field writing), peak memory (sacct MaxRSS), output size.

Usage (from Au_Cl/):  rough_sampling_v1/pyrun_rs.sh size_complex_test_v1/pair_analysis.py
Output: size_complex_test_v1/paired_dft_results.csv, paired_dft_results.md, profiles_<pair>.png
"""
import csv
import json
import os
import re
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from ase.io import read  # noqa: E402
from scipy.ndimage import map_coordinates  # noqa: E402

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
R = f"{ROOT}/rough_sampling_v1"; T = f"{ROOT}/size_complex_test_v1"
sys.path.insert(0, R); sys.path.insert(0, f"{ROOT}/scripts")
import fieldio  # noqa: E402
import analysis_spatial as AS  # noqa: E402  read_field (ngz, ngy, ngx), species(phi, sion) -> (n_minus, n_plus)
from extract_cells import neighbours, R_CORE  # noqa: E402
from rough_dft import FIELD_FILES, MU_TOL, ZVAL_AU  # noqa: E402

DISC_R = 4.0; Z_WIN = 12.0; DZ = 0.25
E_PER_UC_CM2 = 1.0 / 1602.18


def parse_outcar(d, natoms):
    txt = open(f"{d}/OUTCAR").read()
    blk = txt.split("TOTAL-FORCE (eV/Angst)")[-1].split("total drift")[0].splitlines()[2:2 + natoms]
    F = np.array([[float(x) for x in l.split()[3:6]] for l in blk]) if "TOTAL-FORCE" in txt else None
    grid = re.search(r"NGX =\s*(\d+) NGY =\s*(\d+) NGZ =\s*(\d+)", txt); nb = re.search(r"NBANDS=\s*(\d+)", txt)
    loops = [float(x) for x in re.findall(r"LOOP:\s+cpu time\s+[\d.]+: real time\s+([\d.]+)", txt)]
    el = re.search(r"Elapsed time \(sec\):\s+([\d.]+)", txt)
    return dict(F=F, grid="x".join(grid.groups()) if grid else "", nbands=int(nb.group(1)) if nb else None, n_loops=len(loops), scf_time_s=sum(loops),
                elapsed_s=float(el.group(1)) if el else None, finished="General timing" in txt)


def cp_state(d):
    ion = re.findall(r"CPM-ion:.*?N_ele=\s*([-\d.]+)\s+mu_e=\s*([-\d.]+)\s+TARGETMU=\s*([-\d.]+)", open(f"{d}/log.out").read()) if os.path.exists(f"{d}/log.out") else []
    return dict(n_cp=len(ion), N_e=float(ion[-1][0]) if ion else None, mu_e=float(ion[-1][1]) if ion else None, TARGETMU=float(ion[-1][2]) if ion else None)


def sacct(job_id):
    if not job_id: return {}
    out = subprocess.run(["sacct", "-j", str(job_id), "-n", "-o", "MaxRSS,Elapsed,State"], capture_output=True, text=True).stdout.split("\n")
    rss = [l.split()[0] for l in out if l.strip() and l.split()[0][-1] in "KMG"]
    return dict(maxrss=max(rss, key=lambda s: float(s[:-1]) * {"K": 1, "M": 1e3, "G": 1e6}[s[-1]]) if rss else "", elapsed=out[0].split()[1] if out and len(out[0].split()) > 1 else "")


def accept(d, cell_at):
    """label acceptance of one member (as rough_dft.status would judge it)."""
    st = cp_state(d); oc = parse_outcar(d, len(cell_at))
    notes = []; ok = True
    if st["mu_e"] is None: return False, "no CP closure", st, oc
    # closure is judged at the tolerance the run was actually given: the six 2026-10-02 pair runs carried FERMICONVERGE 0.001
    # (judged at 0.0015 eV); everything else, including every run from 2026-10-03 on (user rule: always 0.01), at 0.011
    fermi = re.search(r"FERMICONVERGE\s*=\s*([\d.]+)", open(f"{d}/INCAR").read()); tol = 0.0015 if fermi and float(fermi.group(1)) <= 0.001 else MU_TOL
    if abs(st["mu_e"] - st["TARGETMU"]) > tol: ok = False; notes.append(f"mu_e off target by {st['mu_e'] - st['TARGETMU']:+.4f} (tolerance {tol})")
    if not oc["finished"]: ok = False; notes.append("VASP not finished")
    bad = [f for f in FIELD_FILES if not fieldio.check(f"{d}/{f}", f)[0]]
    if bad: ok = False; notes.append("fields failing: " + ",".join(bad))
    okc, _, stc = fieldio.check(f"{d}/CHGCAR", "CHGCAR")
    if okc and abs(stc["sum_over_grid"] - st["N_e"]) > 0.02: ok = False; notes.append(f"charge closure {stc['sum_over_grid'] - st['N_e']:+.3f} e")
    g, gn = fieldio.contcar_matches_poscar(f"{d}/POSCAR", f"{d}/CONTCAR")
    if not g: ok = False; notes.append("CONTCAR != POSCAR")
    return ok, "; ".join(notes) or "complete", st, oc


def dft_positions(d):
    return read(f"{d}/POSCAR").get_positions()


def seam_distance(at, P):
    c = at.get_cell().array; C = c[:2, :2]; g = P[:, :2] @ np.linalg.inv(C); A = abs(np.linalg.det(C))
    h1 = A / np.linalg.norm(c[1][:2]); h2 = A / np.linalg.norm(c[0][:2])
    return np.minimum(np.minimum(g[:, 0], 1 - g[:, 0]) * h1, np.minimum(g[:, 1], 1 - g[:, 1]) * h2)


def profiles(d, at, centre_xy, z0):
    """disc-averaged PHI, SION, n-, n+ against height above z0 (DFT frame), plus Gamma- and the accessible boundary."""
    phi = AS.read_field(f"{d}/PHI"); sion = AS.read_field(f"{d}/SION"); nm, npl = AS.species(phi, sion)
    cell = read(f"{d}/POSCAR").get_cell().array; inv = np.linalg.inv(cell)
    zs = np.arange(0.0, Z_WIN + 1e-9, DZ); out = {k: [] for k in ("PHI", "SION", "n_minus", "n_plus")}
    r = np.arange(-DISC_R, DISC_R + 1e-9, 0.5); gx, gy = np.meshgrid(r, r); m = gx ** 2 + gy ** 2 <= DISC_R ** 2; pts = np.c_[gx[m], gy[m]]
    for z in zs:
        xyz = np.c_[pts + centre_xy, np.full(len(pts), z0 + z)]; f = (xyz @ inv) % 1.0
        coords = np.array([f[:, 2] * phi.shape[0], f[:, 1] * phi.shape[1], f[:, 0] * phi.shape[2]])
        for k, arr in (("PHI", phi), ("SION", sion), ("n_minus", nm), ("n_plus", npl)):
            out[k].append(float(map_coordinates(arr, coords, order=1, mode="grid-wrap").mean()))
    out = {k: np.array(v) for k, v in out.items()}
    gamma_minus = float(np.trapezoid(out["n_minus"], zs))              # A^-2 (number per area over the window)
    acc = zs[np.argmax(out["SION"] > 0.5)] if (out["SION"] > 0.5).any() else np.nan
    return zs, out, gamma_minus, float(acc)


def main():
    plan = json.load(open(f"{T}/pair_plan.json")); pq = json.load(open(f"{T}/pair_queue.json"))["tasks"] if os.path.exists(f"{T}/pair_queue.json") else []
    rows = []; md = ["# Size-pair results", ""]
    for pr in plan["pairs"]:
        mem = {}
        for key, size in (("six", "6x6"), ("eight", "8x8")):
            if pr[key].get("ref_task_dir"):
                d = pr[key]["ref_task_dir"]; cell_dir = f"{R}/cells/" + os.path.basename(os.path.dirname(d))
                job = next((t["job_id"] for t in json.load(open(f"{R}/dft/queue.json")) if t["dir"] == d), None)
            else:
                t = next((t for t in pq if t["pair"] == pr["id"] and t["member"] == size and "dir" in t), None)
                if t is None: continue
                d = t["dir"]; cell_dir = t["cell_dir"]; job = t.get("job_id")
            if not os.path.exists(f"{d}/OUTCAR"): mem[size] = dict(status="not run"); continue
            at = read(f"{cell_dir}/cell.extxyz"); ok, note, st, oc = accept(d, at)
            mem[size] = dict(status="complete" if ok else note, dir=d, at=at, st=st, oc=oc, job=job, P=dft_positions(d))
        row = dict(pair=pr["id"], target=pr["target"], U_V=pr["U_V"], status_6x6=mem.get("6x6", {}).get("status", "not run"), status_8x8=mem.get("8x8", {}).get("status", "not run"))
        md += [f"## {pr['id']} — {pr['target']} (U = {pr['U_V']:+.2f} V)", "", f"6x6: {row['status_6x6']} | 8x8: {row['status_8x8']}", ""]
        both = all(mem.get(s, {}).get("status") == "complete" for s in ("6x6", "8x8"))
        for s in ("6x6", "8x8"):
            m = mem.get(s)
            if not m or "at" not in m: continue
            sa = sacct(m["job"]); N = len(m["at"]); A = abs(np.linalg.det(m["at"].get_cell().array[:2, :2]))
            sigma = -(m["st"]["N_e"] - ZVAL_AU * N) / A / E_PER_UC_CM2 if m["st"]["N_e"] else None
            row.update({f"atoms_{s}": N, f"kpoints_{s}": open(f"{m['dir']}/KPOINTS").read().splitlines()[3].strip(), f"fft_{s}": m["oc"]["grid"], f"nbands_{s}": m["oc"]["nbands"],
                        f"scf_steps_{s}": m["oc"]["n_loops"], f"cp_rounds_{s}": m["st"]["n_cp"], f"mu_e_{s}": m["st"]["mu_e"], f"N_e_{s}": m["st"]["N_e"], f"sigma_uC_cm2_{s}": round(sigma, 3) if sigma else "",
                        f"scf_time_h_{s}": round(m["oc"]["scf_time_s"] / 3600, 2), f"elapsed_h_{s}": round((m["oc"]["elapsed_s"] or 0) / 3600, 2),
                        f"nonscf_time_h_{s}": round(((m["oc"]["elapsed_s"] or 0) - m["oc"]["scf_time_s"]) / 3600, 2), f"node_h_{s}": round((m["oc"]["elapsed_s"] or 0) / 3600, 2),
                        f"maxrss_{s}": sa.get("maxrss", ""), f"output_GB_{s}": round(sum(os.path.getsize(os.path.join(m["dir"], f)) for f in os.listdir(m["dir"])) / 1e9, 1)})
        if both:
            a6, a8 = mem["6x6"]["at"], mem["8x8"]["at"]; p6, p8 = a6.get_array("parent_id"), a8.get_array("parent_id")
            i6 = {int(p): k for k, p in enumerate(p6)}; i8 = {int(p): k for k, p in enumerate(p8)}
            common = sorted(set(i6) & set(i8)); C = []
            for pid in common:
                n6 = set(int(p6[j]) for j in neighbours(a6, i6[pid], R_CORE)[0]); n8 = set(int(p8[j]) for j in neighbours(a8, i8[pid], R_CORE)[0])
                if n6 == n8: C.append(pid)
            F6, F8 = mem["6x6"]["oc"]["F"], mem["8x8"]["oc"]["F"]
            dF = np.array([F6[i6[p]] - F8[i8[p]] for p in C]); mag = np.linalg.norm(dF, axis=1)
            rms = float(np.sqrt((dF ** 2).sum() / (3 * len(C)))); imax = int(np.argmax(mag))
            seam6 = seam_distance(a6, mem["6x6"]["P"]); dseam = np.array([seam6[i6[p]] for p in C])
            near = dseam < 2 * 2.55
            F8m = np.array([F8[i8[p]] for p in C]); rmsF8 = float(np.sqrt((F8m ** 2).sum() / (3 * len(C))))
            row.update(rms_F_8x8_matched_eV_A=round(rmsF8, 4), rms_dF_over_rms_F=round(rms / rmsF8, 3) if rmsF8 else "")
            # the matched criterion is identity of neighbour IDS; the 6x6 seam band was FLARE-minimised while the same atoms are
            # thermal in the 8x8, so a matched atom's neighbours can sit at different relative positions -> measure that shift
            shift = []
            for pid in C:
                nb = [int(j) for j in neighbours(a6, i6[pid], R_CORE)[0]]; nb8 = [i8[int(p6[j])] for j in nb]
                v6 = a6.get_distances(i6[pid], nb, mic=True, vector=True); v8 = a8.get_distances(i8[pid], nb8, mic=True, vector=True)
                shift.append(float(np.linalg.norm(v6 - v8, axis=1).max()) if nb else 0.0)
            shift = np.array(shift); rigid = shift < 0.02
            row.update(n_matched_rigid_env=int(rigid.sum()),
                       rms_dF_rigid_env=round(float(np.sqrt((dF[rigid] ** 2).sum() / (3 * max(rigid.sum(), 1)))), 4) if rigid.any() else "",
                       rms_dF_shifted_env=round(float(np.sqrt((dF[~rigid] ** 2).sum() / (3 * max((~rigid).sum(), 1)))), 4) if (~rigid).any() else "")
            row.update(n_matched=len(C), rms_dF_eV_A=round(rms, 4), max_dF_eV_A=round(float(mag[imax]), 4), max_dF_parent_atom=C[imax], max_dF_seam_dist_A=round(float(dseam[imax]), 2),
                       rms_dF_near_seam=round(float(np.sqrt((dF[near] ** 2).sum() / (3 * max(near.sum(), 1)))), 4), rms_dF_interior=round(float(np.sqrt((dF[~near] ** 2).sum() / (3 * max((~near).sum(), 1)))), 4),
                       n_near_seam=int(near.sum()))
            # fields in the same region above the core (aligned on the centre atom; the core is identical by construction)
            c6 = json.load(open(f"{pr['six']['cell_dir']}/cut_info.json"))["centre_in_cell"] if os.path.exists(f"{pr['six']['cell_dir']}/cut_info.json") else None
            if c6 is not None:
                pid_c = int(p6[c6]); k8 = i8[pid_c]
                prof = {}
                for s, m, k in (("6x6", mem["6x6"], c6), ("8x8", mem["8x8"], k8)):
                    P = m["P"]; z0 = P[k, 2]; zs, out, gm, acc = profiles(m["dir"], m["at"], P[k, :2], z0); prof[s] = (zs, out); row[f"gamma_minus_A-2_{s}"] = round(gm, 5); row[f"access_boundary_A_{s}"] = acc
                # electrolyte-side differences (6x6 - 8x8) on the common height grid, from 4 A above the centre atom
                zs6, o6 = prof["6x6"]; zs8, o8 = prof["8x8"]; sel = zs6 >= 4.0
                diff = {k: o6[k] - o8[k] for k in o6}
                row["dPHI_max_meV_above4A"] = round(float(np.abs(diff["PHI"][sel]).max()) * 1000, 1)
                row["dSION_max_above4A"] = round(float(np.abs(diff["SION"][sel]).max()), 3)
                # the n- difference is quoted against the 8x8 peak in the window (a point-wise ratio at the onset, where both are ~0, is meaningless)
                row["dn_minus_max_A-3_above4A"] = float(np.abs(diff["n_minus"][sel]).max())
                row["dn_minus_max_over_peak8x8"] = round(row["dn_minus_max_A-3_above4A"] / max(float(np.abs(o8["n_minus"][sel]).max()), 1e-12), 3)
                row["gamma_minus_rel_diff"] = round((row["gamma_minus_A-2_6x6"] - row["gamma_minus_A-2_8x8"]) / row["gamma_minus_A-2_8x8"], 3) if row.get("gamma_minus_A-2_8x8") else ""
                with plt.rc_context({"font.size": 15, "axes.titlesize": 17, "axes.labelsize": 15, "legend.fontsize": 14}):
                    fig, axes = plt.subplots(2, 4, figsize=(24, 11), sharex=True)
                    labs = (("PHI", "PHI (eV)", "dPHI (meV)", 1000.0), ("SION", "S_ion", "dS_ion", 1.0), ("n_minus", "n- (A^-3)", "dn- (A^-3)", 1.0), ("n_plus", "n+ (A^-3)", "dn+ (A^-3)", 1.0))
                    for col, (k, lab, dlab, scale) in enumerate(labs):
                        ax = axes[0, col]
                        for s, (zs, out) in prof.items(): ax.plot(zs, out[k], label=s, lw=2.2)
                        ax.set_title(lab); ax.legend(); ax.grid(alpha=0.3)
                        if k == "PHI": ax.set_ylim(-1.0, 0.5); ax.set_title("PHI (eV), electrolyte side; metal dip clipped")
                        ax = axes[1, col]; ax.plot(zs6, diff[k] * scale, color="k", lw=2.2); ax.axhline(0, color="0.6", lw=1); ax.axvline(4.0, color="0.6", lw=1, ls="--")
                        ax.set_title(f"{dlab}: 6x6 - 8x8"); ax.set_xlabel("height above the centre atom (A)"); ax.grid(alpha=0.3)
                        if k == "PHI": ax.set_ylim(-60, 60)
                    fig.suptitle(f"{pr['id']}: disc r = {DISC_R} A above the centre atom, U = {pr['U_V']:+.2f} V  (bottom row: differences; dashed line = 4 A, where the comparison starts)", fontsize=17)
                    fig.tight_layout(); fig.savefig(f"{T}/profiles_{pr['id']}.png", dpi=100); plt.close(fig)
                md.append(f"![profiles](profiles_{pr['id']}.png)")
                md += ["", f"electrolyte-side differences (z >= 4 A above the centre atom): max |dPHI| {row['dPHI_max_meV_above4A']} meV, max |dS_ion| {row['dSION_max_above4A']}, "
                       f"max |dn-| {row['dn_minus_max_A-3_above4A']:.1e} A^-3 = {row['dn_minus_max_over_peak8x8']:.0%} of the 8x8 peak in the window; Gamma- differs by {row['gamma_minus_rel_diff']:+.1%}."]
            order = np.argsort(-mag)
            md += ["", f"matched atoms |C| = {len(C)} (of {len(common)} shared; identical neighbour sets within R_CORE = {R_CORE} A in both cells); RMS dF = {rms:.4f} eV/A "
                   f"against RMS |F| = {rmsF8:.4f} eV/A on the same atoms in the 8x8 (ratio {rms / rmsF8:.2f}), max {mag[imax]:.4f} at parent atom {C[imax]} "
                   f"({dseam[imax]:.1f} A from the nearest 6x6 seam); RMS within 2 rows of a seam {row['rms_dF_near_seam']:.4f} ({near.sum()} atoms) vs interior {row['rms_dF_interior']:.4f}.",
                   f"Neighbour positions: {int(rigid.sum())} matched atoms have every neighbour within {R_CORE} A at the same relative position in both cells (< 0.02 A), "
                   f"RMS dF {row['rms_dF_rigid_env']} eV/A; the other {int((~rigid).sum())} have neighbours that the 6x6 seam repair moved (max shift up to {shift.max():.2f} A), RMS dF {row['rms_dF_shifted_env']} eV/A.", "",
                   "| parent atom | layer (6x6) | dist. to 6x6 seam (A) | max neighbour shift (A) | |F| 8x8 (eV/A) | |dF| (eV/A) |", "|---|---|---|---|---|---|"]
            lay6 = a6.get_array("layer") if "layer" in a6.arrays else None
            md += [f"| {C[j]} | {int(lay6[i6[C[j]]]) if lay6 is not None else ''} | {dseam[j]:.1f} | {shift[j]:.3f} | {np.linalg.norm(F8m[j]):.3f} | {mag[j]:.3f} |" for j in order]
            md += ["", f"sigma 6x6 {row.get('sigma_uC_cm2_6x6')} vs 8x8 {row.get('sigma_uC_cm2_8x8')} uC/cm2 (supplement); Gamma- {row.get('gamma_minus_A-2_6x6')} vs {row.get('gamma_minus_A-2_8x8')} A^-2; "
                   f"accessible boundary {row.get('access_boundary_A_6x6')} vs {row.get('access_boundary_A_8x8')} A above the centre.",
                   f"cost: 6x6 {row.get('scf_steps_6x6')} SCF steps / {row.get('cp_rounds_6x6')} CP rounds, SCF {row.get('scf_time_h_6x6')} h, elapsed {row.get('elapsed_h_6x6')} h, non-SCF {row.get('nonscf_time_h_6x6')} h, "
                   f"MaxRSS {row.get('maxrss_6x6')}, output {row.get('output_GB_6x6')} GB | 8x8 {row.get('scf_steps_8x8')} / {row.get('cp_rounds_8x8')}, SCF {row.get('scf_time_h_8x8')} h, elapsed {row.get('elapsed_h_8x8')} h, "
                   f"non-SCF {row.get('nonscf_time_h_8x8')} h, MaxRSS {row.get('maxrss_8x8')}, output {row.get('output_GB_8x8')} GB -> 6x6 saves {(row.get('elapsed_h_8x8') or 0) - (row.get('elapsed_h_6x6') or 0):.1f} node-h on this pair."]
        rows.append(row); md.append("")
    keys = sorted({k for r in rows for k in r}, key=lambda k: (k not in ("pair", "target", "U_V"), k))
    with open(f"{T}/paired_dft_results.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); [w.writerow(r) for r in rows]
    open(f"{T}/paired_dft_results.md", "w").write("\n".join(md) + "\n"); print("\n".join(md))


if __name__ == "__main__":
    main()
