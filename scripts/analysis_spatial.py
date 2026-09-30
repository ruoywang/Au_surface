#!/usr/bin/env python3
"""Where do the anions actually pile up? Morphology-resolved spatial maps from the self-consistent electrolyte fields.

For every requested state this reduces the 3D fields to compact per-column (x,y) arrays and a few scalars, so the
cross-structure comparison never has to hold a full grid. Nothing here adds explicit ions: the enrichment comes from
the converged electrostatic potential PHI and the ion-accessibility mask SION through the model's own constitutive
relation (the same reconstruction validated in scripts/validate_ion_reconstruction.py).

    u          = PHI / kT
    D          = 1 + (2 n_b / n_max) (cosh u - 1)          finite-size (lattice-gas) saturation
    n_-(r)     = S_ion n_b exp(-u) / D                     anion number density
    n_+(r)     = S_ion n_b exp(+u) / D

Two DIFFERENT questions are reported, because they are not interchangeable:

    K_Omega    = int_Omega n_- dV / ( n_b int_Omega S_ion dV )
                 "inside this accessible pocket, how many times the bulk concentration?"
    Gamma_-    = (1/A_proj) int_Omega [ n_- - n_b S_ion ] dV
                 "per unit projected electrode area, how many EXTRA anions does this region hold?"
A tight pocket can have a large K and hold very few extra ions; a broad region can have K ~ 1 and dominate Gamma.

Integration region. S_ion changes slightly between potentials (the cavity follows the electron density), so the
-n_b S_ion term is kept. Two conventions, both recorded, following the convention already established in
scripts/step_width_excess.py:
  (A) ABSOLUTE   Omega = {z < Z_UP}: column values partition the cell exactly -> use for totals and closure.
  (B) BOUNDARY-RELATIVE  per column Omega = {z_b(x,y) < z < z_b(x,y) + DEPTH}, z_b = the local S_ion = 0.5 crossing
      -> every column keeps the same fraction of the decaying tail -> use for comparing regions and structures whose
      surfaces sit at different heights (an island top is ~2.4 A above its terrace).

Morphology resolution without bespoke per-structure masks: each column is assigned to the nearest surface atom in xy
(minimum image) and labelled by that atom's bulk-coordination number, computed at 3.4 A. CN then separates terrace
(CN 9), step edge and island/pit rim (CN 7-8), and kink / corner / adatom (CN <= 6) uniformly across every structure.
The per-column maps are also saved so bespoke regions can be cut later without re-reading the 3D fields.

Usage (from Au_Cl/):
    scripts/pyrun.sh scripts/analysis_spatial.py --states <n>|all [--workers N] [--transmission]
Outputs: analysis/spatial/percolumn/<state_id>.npz  and  analysis/spatial/spatial.json
"""
import argparse
import collections
import json
import os
import sys
import time

import numpy as np
from ase.io import read

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
OUT = f"{ROOT}/analysis/spatial"
KBT = 8.6173857e-5 * 298.0
N_BULK = 1.0 * 6.022e-4                    # 1 M in A^-3
R_ION = 4.0
N_MAX = 1.0 / (2 ** (5 / 6) * R_ION) ** 3
Z_UP = 31.0                                # last plane with S_ion = 1 everywhere in this window (audit K / step_width_excess)
DEPTH = 10.4                               # boundary-relative depth the raised terraces can still accommodate below Z_UP
CN_CUT = 3.4                               # A, first-shell cutoff for the coordination number (a = 4.158 -> d_nn = 2.94)


def read_field(path, dtype=np.float32):
    with open(path) as f:
        lines = f.readlines()
    for i, l in enumerate(lines):
        p = l.split()
        if len(p) == 3 and all(x.isdigit() for x in p):
            ngx, ngy, ngz = map(int, p); di = i; break
    else:
        raise ValueError(f"no grid dimensions in {path}")
    flat = np.fromstring("".join(lines[di + 1:]), sep=" ", dtype=np.float64)[:ngx * ngy * ngz]
    return flat.astype(dtype).reshape((ngz, ngy, ngx))


def species(phi, sion):
    u = np.clip(phi / KBT, -50, 50)
    D = 1 + (2 * N_BULK / N_MAX) * (np.cosh(u) - 1)
    return sion * N_BULK * np.exp(-u) / D, sion * N_BULK * np.exp(u) / D


def coordination(atoms):
    pos = atoms.get_positions(); cell = atoms.get_cell().array
    n = len(pos); cn = np.zeros(n, int)
    shifts = [i * cell[0] + j * cell[1] for i in (-1, 0, 1) for j in (-1, 0, 1)]
    for s in shifts:
        d = np.linalg.norm(pos[:, None, :] - (pos[None, :, :] + s), axis=-1)
        cn += ((d < CN_CUT) & (d > 0.1)).sum(1)
    return cn


def column_labels(atoms, ngx, ngy):
    """For each grid column, the coordination number of the nearest SURFACE atom in xy (minimum image)."""
    cell = atoms.get_cell().array; pos = atoms.get_positions(); cn = coordination(atoms)
    top = pos[:, 2].max()
    sel = np.where(pos[:, 2] > top - 3.0)[0]                       # the outermost ~one layer, incl. adatoms and island tops
    fx, fy = np.meshgrid((np.arange(ngx) + 0.5) / ngx, (np.arange(ngy) + 0.5) / ngy, indexing="ij")
    gxy = fx[..., None] * cell[0][:2] + fy[..., None] * cell[1][:2]        # (ngx, ngy, 2) cartesian
    best = np.full(gxy.shape[:2], np.inf); lab = np.zeros(gxy.shape[:2], int); hgt = np.zeros(gxy.shape[:2])
    for i in sel:
        for si in (-1, 0, 1):
            for sj in (-1, 0, 1):
                p = pos[i, :2] + si * cell[0][:2] + sj * cell[1][:2]
                d = np.linalg.norm(gxy - p, axis=-1)
                m = d < best
                best[m] = d[m]; lab[m] = cn[i]; hgt[m] = pos[i, 2]
    return lab.T, hgt.T, best.T, cn, sel                                   # transposed to (ngy, ngx) = field order


def reduce_state(d, want_metal=False):
    """Reduce one run directory to per-column arrays + scalars. Returns a dict of numpy arrays and scalars."""
    atoms = read(f"{d}/CONTCAR" if os.path.basename(d).startswith("relax__") and os.path.exists(f"{d}/CONTCAR")
                 else f"{d}/POSCAR")
    cell = atoms.get_cell().array
    A_proj = float(np.linalg.norm(np.cross(cell[0], cell[1])))
    V = float(atoms.get_volume())
    phi = read_field(f"{d}/PHI"); sion = read_field(f"{d}/SION")
    ngz, ngy, ngx = phi.shape
    dz = cell[2][2] / ngz
    dV = V / (ngx * ngy * ngz)
    z = (np.arange(ngz) + 0.5) * cell[2][2] / ngz
    nm, npl = species(phi, sion)
    excess = nm - N_BULK * sion                                   # e/A^3 -> anion NUMBER excess density

    kup = int(np.searchsorted(z, Z_UP))
    col_area = A_proj / (ngx * ngy)
    # (A) absolute: integrate each column from 0 to Z_UP
    gam_abs = excess[:kup].sum(0) * dz                            # (ngy, ngx), anions per A^2 of column
    nm_abs = nm[:kup].sum(0) * dz
    si_abs = sion[:kup].sum(0) * dz

    # (B) boundary-relative: per column, from the S_ion = 0.5 crossing up by DEPTH
    zb_idx = np.argmax(sion > 0.5, axis=0)                        # first index where the column becomes accessible
    zb = z[zb_idx]
    ndep = int(round(DEPTH / dz))
    kk = np.arange(ndep)[:, None, None] + zb_idx[None, :, :]
    kk = np.clip(kk, 0, ngz - 1)
    gam_rel = np.take_along_axis(excess, kk, 0).sum(0) * dz
    nm_rel = np.take_along_axis(nm, kk, 0).sum(0) * dz
    si_rel = np.take_along_axis(sion, kk, 0).sum(0) * dz
    fits = float(np.mean(zb + DEPTH <= Z_UP))

    lab, hgt, dist, cn, sel = column_labels(atoms, ngx, ngy)
    out = dict(A_proj=A_proj, n_atoms=len(atoms), ngx=ngx, ngy=ngy, ngz=ngz, dz=dz, Z_UP=Z_UP, DEPTH=DEPTH,
               cols_fitting_depth=fits, cell=cell,
               gamma_abs=gam_abs, gamma_rel=gam_rel, nminus_abs=nm_abs, nminus_rel=nm_rel,
               sion_abs=si_abs, sion_rel=si_rel, zb=zb, cn_label=lab, top_height=hgt, nn_distance=dist,
               metal_top_z=float(atoms.get_positions()[:, 2].max()),
               cn_values=cn, surface_atom_idx=sel,
               # cell totals (absolute convention partitions exactly)
               Gamma_total=float(gam_abs.mean() * 1.0),                       # per A^2 of projected area
               K_cell=float(nm_abs.sum() / (N_BULK * si_abs.sum())) if si_abs.sum() > 0 else float("nan"),
               phi_max=float(phi[:kup].max()), phi_min=float(phi[:kup].min()))
    if want_metal:
        chg = read_field(f"{d}/CHGCAR") / V                        # CHGCAR holds rho*V
        rhob = read_field(f"{d}/RHOB") / V
        kmet = int(np.searchsorted(z, out["metal_top_z"] + 2.0))
        out["ne_col"] = chg[:kmet].sum(0) * dz                     # electrons per A^2 in the metal + 2 A
        out["rhob_col"] = rhob.sum(0) * dz
        zref = out["metal_top_z"] + 4.0
        out["phi_plane"] = phi[int(np.searchsorted(z, zref))]
        out["phi_plane_z"] = zref
    return out


def aggregate(r):
    """Region-resolved K and Gamma by coordination class, boundary-relative convention."""
    lab = r["cn_label"]; ncol = lab.size
    classes = {"kink/adatom (CN<=6)": lab <= 6, "edge/rim (CN 7-8)": (lab >= 7) & (lab <= 8),
               "terrace (CN 9)": lab == 9, "over-coordinated (CN>=10)": lab >= 10}
    out = {}
    for name, m in classes.items():
        if not m.any(): continue
        out[name] = dict(area_fraction=float(m.mean()),
                         K=float(r["nminus_rel"][m].sum() / (N_BULK * r["sion_rel"][m].sum())),
                         Gamma_per_A2=float(r["gamma_rel"][m].mean()),
                         Gamma_share=float(r["gamma_rel"][m].sum() / r["gamma_rel"].sum()) if r["gamma_rel"].sum() else None)
    out["__all__"] = dict(area_fraction=1.0,
                          K=float(r["nminus_rel"].sum() / (N_BULK * r["sion_rel"].sum())),
                          Gamma_per_A2=float(r["gamma_rel"].mean()), Gamma_share=1.0)
    return out


def work(job):
    sid, d, want_metal = job
    try:
        t0 = time.time()
        r = reduce_state(d, want_metal)
        np.savez_compressed(f"{OUT}/percolumn/{sid}.npz", **{k: v for k, v in r.items() if isinstance(v, np.ndarray)})
        scal = {k: v for k, v in r.items() if not isinstance(v, np.ndarray)}
        scal["regions"] = aggregate(r); scal["seconds"] = round(time.time() - t0, 1)
        return sid, scal, None
    except Exception as e:
        return sid, None, f"{type(e).__name__}: {e}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--states", default="all")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--transmission", action="store_true", help="also read CHGCAR and RHOB (metal and bound charge)")
    ap.add_argument("--campaign", default="dataset_plan_v1_rev2")
    a = ap.parse_args()
    os.makedirs(f"{OUT}/percolumn", exist_ok=True)
    S = json.load(open(f"{ROOT}/dataset_v1/states.json"))["states"]
    sel = [s for s in S.values() if a.campaign in ("any", s["campaign"])]
    sel.sort(key=lambda s: (s["structure_id"], s["config"], s["electronic_state"]["TARGETMU_eV"]))
    if a.states != "all": sel = sel[:int(a.states)]
    done = {f[:-4] for f in os.listdir(f"{OUT}/percolumn")} if os.path.exists(f"{OUT}/percolumn") else set()
    jobs = [(s["state_id"], s["source_dir"], a.transmission) for s in sel if s["state_id"] not in done]
    print(f"{len(sel)} states selected, {len(sel)-len(jobs)} already reduced, {len(jobs)} to do, {a.workers} workers")
    res, errs = {}, {}
    if jobs:
        import multiprocessing as mp
        with mp.Pool(a.workers) as pool:
            for i, (sid, scal, err) in enumerate(pool.imap_unordered(work, jobs), 1):
                if err: errs[sid] = err
                else: res[sid] = scal
                if i % 20 == 0 or i == len(jobs): print(f"  {i}/{len(jobs)}  last={sid} {'ERR '+err if err else ''}", flush=True)
    path = f"{OUT}/spatial.json"
    old = json.load(open(path)) if os.path.exists(path) else dict(states={}, errors={})
    old["states"].update(res); old["errors"].update(errs); old["updated"] = time.strftime("%Y-%m-%d %H:%M")
    old["convention"] = dict(Z_UP=Z_UP, DEPTH=DEPTH, CN_CUT=CN_CUT, n_bulk_A3=N_BULK, R_ION=R_ION)
    json.dump(old, open(path, "w"), indent=1)
    print(f"reduced {len(res)} states, {len(errs)} errors -> {path}")
    for k, v in list(errs.items())[:5]: print("  ERR", k, v)


if __name__ == "__main__":
    main()
