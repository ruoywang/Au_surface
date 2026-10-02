#!/usr/bin/env python3
"""Cut a small periodic electrode out of a large parent surface around a chosen centre, keeping the core exactly
and repairing only the periphery. This is the step the whole plan turns on (section 6).

Core: the centre atom and everything within 6 A of it in 3D (the union over several key atoms if given). The
core is PROTECTED: its atoms are carried over with their IDs and must survive to 1e-3 A after one global
translation, and the 6 A neighbourhood of the centre must be the same set of atoms afterwards. 6 A is a
protection radius, not a physical interaction cutoff, and nothing claims that keeping the core keeps the
parent's potential.

Cell: 8 x 8 (256 base Au) first, then 10 x 8 (320), both commensurate with the parent's 32 x 32 lattice. The
cut origin is chosen among lattice-compatible origins that keep the core away from the seam, by the smallest
mismatch between the two sides of each seam (site occupancy and height); the whole cell is kept continuous
metal with periodic in-plane boundaries, never a metal block with vacuum on four sides.

Repair, in order: try other origins -> try the other cell size -> drop the candidate. The buffer (everything
in the top two levels farther than the protection radius from the core and within 2 rows of a seam) is the
only place occupancy may change, every added / removed atom is recorded, and a minimisation with the core and
the bottom two layers fixed relaxes the buffer with the candidate-generation potential. A smaller force is not
accepted as proof of a correct occupancy: the geometric checks below decide.

Checks that must pass: core preserved; neighbourhood 1:1; no duplicate atoms or contacts < 2.5 A; no vacuum
slit or floating cluster (every atom has >= 3 neighbours within 3.4 A, bottom layer complete); pit depth <= 1
layer; mapping invariant under r -> r + a, r + b; the repaired periphery is no more unusual than the core
(no atom in the buffer with fewer neighbours than the least-coordinated core atom minus one).

The CP-DFT label later attached to this cell belongs to THIS reconstructed periodic electrode, not to the
parent, and its total energy is not a "core energy".

Usage (from Au_Cl/):  env PYTHONPATH=rough_sampling_v1/env/pylib:.pyshim <python> rough_sampling_v1/extract_cells.py
                      [--centres rough_sampling_v1/centres/centres.jsonl] [--limit N] [--no-repair]
Output: rough_sampling_v1/cells/<cell_id>/POSCAR, cells_manifest.jsonl, cells_summary.md
"""
import argparse
import json
import os
import subprocess
import sys

import numpy as np
from ase import Atoms
from ase.io import read, write

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
R = f"{ROOT}/rough_sampling_v1"
sys.path.insert(0, R)
from select_centres import surface_atoms  # noqa: E402  (same exposure rule as the centre selection)
from lmpio import write_data, group_lines, read_data, mic_displacement  # noqa: E402
POT = f"{R}/potential/Au_training/lmp_t0.0001_no_bulk_vac_fix3.header2025.flare"
LMP = f"{R}/env/src/lammps-22Jul2025/build_mpi/lmp"
A0 = 4.158; NN = A0 / np.sqrt(2); D111 = A0 / np.sqrt(3)
R_CORE = 6.0
CELLS = [(8, 8), (10, 8)]
SEAM_ROWS = 2
MIN_CONTACT = 2.5
MAX_CORE_SHIFT = 1e-3
CELLDIR = f"{R}/cells"


def lattice(at_parent):
    """The parent's primitive in-plane vectors (32 x 32 of them) and its site origin."""
    c = at_parent.get_cell().array
    return c[0][:2] / 32.0, c[1][:2] / 32.0, c


def boundary_offset(at, a1, a2):
    """Fractional offset (in the primitive basis) that puts a cut boundary midway between occupied lattice
    rows. fcc(111) layers sit at three positions a third of a primitive vector apart, so a boundary through an
    atom row would split that row between the two sides under thermal noise; the offset is measured from the
    structure, not assumed."""
    f = (at.get_positions()[:, :2] @ np.linalg.inv(np.array([a1, a2]))) % 1.0
    off = []
    for ax in (0, 1):
        h, edges = np.histogram(f[:, ax], bins=60, range=(0, 1))
        empty = h == 0
        if not empty.any(): off.append(0.0); continue
        # longest circular run of empty bins -> its middle
        idx = np.flatnonzero(empty); runs = np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1)
        if len(runs) > 1 and runs[0][0] == 0 and runs[-1][-1] == 59: runs[0] = np.r_[runs[-1] - 60, runs[0]]; runs = runs[:-1]
        r = max(runs, key=len); off.append(((r[0] + r[-1] + 1) / 2 / 60) % 1.0)
    return np.array(off)


def cut(at, centre_xy, n1, n2, origin_ij, a1, a2, cell, off):
    """Atoms whose in-plane position (minimum image in the parent) falls inside the n1 x n2 parallelogram with
    the given lattice origin. Returns the sub-Atoms with a 'parent_id' array and the new cell."""
    P = at.get_positions()
    o = (origin_ij[0] + off[0]) * a1 + (origin_ij[1] + off[1]) * a2
    Cp = np.array([cell[0][:2], cell[1][:2]]); Cpi = np.linalg.inv(Cp)
    d = P[:, :2] - o
    f = d @ Cpi; f -= np.floor(f + 1e-9)                   # wrap into the parent cell measured from the origin
    d = f @ Cp
    S = np.array([n1 * a1, n2 * a2]); Si = np.linalg.inv(S)
    g = d @ Si
    inside = (g[:, 0] >= -1e-6) & (g[:, 0] < 1 - 1e-6) & (g[:, 1] >= -1e-6) & (g[:, 1] < 1 - 1e-6)
    idx = np.flatnonzero(inside)
    newcell = np.array([[*(n1 * a1), 0.0], [*(n2 * a2), 0.0], cell[2]])
    pos = np.c_[d[idx], P[idx, 2]]
    sub = Atoms("Au" * len(idx), positions=pos, cell=newcell, pbc=(True, True, False))
    sub.set_array("parent_id", idx)
    for key in ("layer", "fixed"):
        if key in at.arrays: sub.set_array(key, at.get_array(key)[idx])
    return sub


def neighbours(at, i, r):
    """Indices within r of atom i, minimum image in-plane."""
    P = at.get_positions(); cell = at.get_cell().array
    best = np.full(len(P), np.inf)
    for si in (-1, 0, 1):
        for sj in (-1, 0, 1):
            sh = si * cell[0] + sj * cell[1]
            best = np.minimum(best, np.linalg.norm(P + sh - P[i], axis=1))
    best[i] = np.inf
    return np.flatnonzero(best < r), best


def all_min_dist_and_cn(at, rcut=3.4):
    P = at.get_positions(); cell = at.get_cell().array; n = len(P)
    best = np.full(n, np.inf); cn = np.zeros(n, int)
    for si in (-1, 0, 1):
        for sj in (-1, 0, 1):
            sh = si * cell[0] + sj * cell[1]
            for a in range(0, n, 256):
                d = np.linalg.norm(P[a:a + 256, None, :] - (P[None, :, :] + sh), axis=-1)
                d[d < 0.1] = np.inf
                best[a:a + 256] = np.minimum(best[a:a + 256], d.min(1)); cn[a:a + 256] += (d < rcut).sum(1)
    return best, cn


def seam_mismatch(sub, n1, n2):
    """How different the two sides of each periodic seam are: compares the occupancy and height of the top
    two levels in a SEAM_ROWS-wide band at each edge with the band at the opposite edge, image-shifted."""
    P = sub.get_positions(); cell = sub.get_cell().array
    zt = P[:, 2].max(); top = P[:, 2] > zt - 2.5 * D111
    S = np.array([cell[0][:2], cell[1][:2]]); g = P[:, :2] @ np.linalg.inv(S)
    score = 0.0
    for ax, n in ((0, n1), (1, n2)):
        w = SEAM_ROWS / n
        lo = top & (g[:, ax] < w); hi = top & (g[:, ax] >= 1 - w)
        # height histograms per level on each side, compared after shifting the hi side by one cell
        hl = np.round((P[lo, 2] - zt) / D111); hh = np.round((P[hi, 2] - zt) / D111)
        for lev in (0, -1, -2):
            score += abs(int((hl == lev).sum()) - int((hh == lev).sum()))
    return score


def core_ok(parent, sub, core_parent_ids, centre_parent_id):
    """Core carried over exactly (after one translation) and the centre's 6 A neighbourhood is the same set."""
    pid = sub.get_array("parent_id"); where = {int(p): k for k, p in enumerate(pid)}
    if any(int(c) not in where for c in core_parent_ids): return False, "core atom missing from the cut"
    Pp = parent.get_positions(); Ps = sub.get_positions()
    cs = Ps[[where[int(c)] for c in core_parent_ids]]; cp = Pp[core_parent_ids]
    # one global translation: compare pairwise differences (translation-free)
    dp = cp - cp[0]; ds = cs - cs[0]
    # in-plane minimum image of the differences in the PARENT cell
    cellp = parent.get_cell().array; C = np.array([cellp[0][:2], cellp[1][:2]]); Ci = np.linalg.inv(C)
    def mic(v): f = v[:, :2] @ Ci; f -= np.round(f); v = v.copy(); v[:, :2] = f @ C; return v
    err = np.abs(mic(dp) - mic(ds)).max()
    if err > MAX_CORE_SHIFT: return False, f"core distorted by {err:.2e} A"
    nb_p, _ = neighbours(parent, int(centre_parent_id), R_CORE)
    nb_s, _ = neighbours(sub, where[int(centre_parent_id)], R_CORE)
    set_p = set(int(x) for x in nb_p); set_s = set(int(pid[k]) for k in nb_s)
    if set_p != set_s:
        return False, f"neighbourhood changed: {len(set_p - set_s)} lost, {len(set_s - set_p)} gained"
    return True, "core preserved to %.1e A, neighbourhood identical (%d atoms)" % (err, len(set_p))


def image_invariance(sub, centre_idx):
    """Shifting every atom by a lattice vector must not change the neighbourhood set."""
    base, _ = neighbours(sub, centre_idx, R_CORE)
    cell = sub.get_cell().array
    for sh in (cell[0], cell[1], cell[0] + cell[1]):
        t = sub.copy(); t.set_positions(t.get_positions() + sh)
        nb, _ = neighbours(t, centre_idx, R_CORE)
        if set(nb.tolist()) != set(base.tolist()): return False
    return True


def layers_from_z(at):
    """Layer index of every atom from its height above the (fixed, hence exactly known) bottom layer. After
    MD the parent's original layer labels are stale, so they are re-derived from z for every frame."""
    z = at.get_positions()[:, 2]
    return np.round((z - z.min()) / D111).astype(int)


def geometry_checks(sub, core_idx, n1, n2):
    md, cn = all_min_dist_and_cn(sub)
    lay = layers_from_z(sub)
    problems = []
    if md.min() < MIN_CONTACT: problems.append(f"contact {md.min():.2f} A < {MIN_CONTACT}")
    if (cn < 3).any(): problems.append(f"{int((cn < 3).sum())} atoms with < 3 neighbours (floating or slit)")
    n_bottom = int((lay == 0).sum())
    if n_bottom != n1 * n2: problems.append(f"bottom layer incomplete ({n_bottom}/{n1 * n2})")
    # pits may be one layer deep: the terrace is layer 3 (of 0..3 base layers, 4 = adatom level), so no
    # exposed atom may lie in layer 1 or below
    exposed = surface_atoms(sub)
    if len(exposed) and lay[exposed].min() < 2: problems.append(f"exposed atom in layer {int(lay[exposed].min())} (pit deeper than one layer)")
    cn_core_min = cn[core_idx].min() if len(core_idx) else 0
    buf = np.setdiff1d(np.arange(len(sub)), core_idx)
    if len(buf) and cn[buf].min() < cn_core_min - 1:
        problems.append(f"periphery has an atom less coordinated ({cn[buf].min()}) than the core minimum ({cn_core_min}) minus one")
    return problems, dict(min_dist_A=float(md.min()), cn_min=int(cn.min()), cn_core_min=int(cn_core_min))


def repair_minimise(sub, fixed_mask, workdir):
    """Minimise the buffer with the candidate potential, core and bottom two layers fixed. Returns the relaxed
    Atoms, or None if the LAMMPS binary is not available."""
    if not os.path.exists(LMP): return None
    os.makedirs(workdir, exist_ok=True)
    write_data(sub, f"{workdir}/data.lammps", comment="reconstructed cell; core and bottom two layers frozen by id")
    script = f"""units metal
atom_style atomic
boundary p p f
newton on
read_data data.lammps
pair_style flare
pair_coeff * * {POT}
{group_lines("frozen", fixed_mask)}
group mobile subtract all frozen
fix hold frozen setforce 0.0 0.0 0.0
min_style cg
minimize 1e-8 1e-4 500 5000
write_data min.data
"""
    open(f"{workdir}/in.lammps", "w").write(script)
    p = subprocess.run([LMP, "-in", "in.lammps", "-log", "log.lammps"], cwd=workdir, capture_output=True, text=True)
    if p.returncode != 0: return None
    a_in = read_data(f"{workdir}/data.lammps"); m = read_data(f"{workdir}/min.data")
    # the parents are lower-triangular so the LAMMPS frame is the ASE frame; apply the (unwrapped) displacement
    r = sub.copy(); r.set_positions(sub.get_positions() + mic_displacement(a_in.get_cell().array, a_in.get_positions(), m.get_positions()))
    return r


def extract_one(parent, centre_atom, cid, repair=True, sizes=CELLS, first_only=True):
    """Returns ([(size, Atoms, info), ...] for the sizes that passed, trial log). With first_only the smaller
    cell is taken as soon as it passes; otherwise every size is tried (test centres, for size pairs)."""
    a1, a2, cellp = lattice(parent)
    off = boundary_offset(parent, a1, a2)
    Pp = parent.get_positions()
    core_ids, _ = neighbours(parent, centre_atom, R_CORE); core_ids = np.append(core_ids, centre_atom)
    cxy = Pp[centre_atom, :2]
    log = []; passed = []
    for (n1, n2) in sizes:
        # lattice-compatible origins that put the centre near the middle of the cut: scan a few around it
        f_c = (cxy @ np.linalg.inv(np.array([a1, a2])))
        base_ij = np.round(f_c - off - np.array([n1 / 2, n2 / 2])).astype(int)
        trials = []
        for di in range(-2, 3):
            for dj in range(-2, 3):
                oij = base_ij + np.array([di, dj])
                sub = cut(parent, cxy, n1, n2, oij, a1, a2, cellp, off)
                ok, why = core_ok(parent, sub, core_ids, centre_atom)
                if not ok: continue
                trials.append((seam_mismatch(sub, n1, n2), tuple(int(x) for x in oij), sub))
        trials.sort(key=lambda t: t[0])
        for score, oij, sub in trials[:5]:
            pid = sub.get_array("parent_id"); where = {int(p): k for k, p in enumerate(pid)}
            core_idx = np.array([where[int(c)] for c in core_ids]); centre_idx = where[int(centre_atom)]
            lay = layers_from_z(sub)
            fixed = (lay < 2); fixed[core_idx] = True
            problems, stats = geometry_checks(sub, core_idx, n1, n2)
            repaired = None; rep_note = "no repair needed" if not problems else "unrepaired"
            if problems and repair:
                repaired = repair_minimise(sub, fixed, f"{CELLDIR}/{cid}/repair_{n1}x{n2}_{oij[0]}_{oij[1]}")
                if repaired is not None:
                    ok2, why2 = core_ok(parent, repaired, core_ids, centre_atom)
                    problems2, stats2 = geometry_checks(repaired, core_idx, n1, n2)
                    if ok2 and not problems2:
                        sub, problems, stats, rep_note = repaired, [], stats2, "buffer minimised with core and bottom layers fixed; occupancy unchanged"
                    else:
                        rep_note = f"minimisation did not clear: {problems2 or why2}"
            inv = image_invariance(sub, centre_idx)
            log.append(dict(cell=f"{n1}x{n2}", origin=oij, seam_score=score, problems=problems, invariant=inv, repair=rep_note))
            if not problems and inv:
                sub.set_array("layer", layers_from_z(sub)); sub.set_array("fixed", (layers_from_z(sub) < 2).astype(int))
                passed.append((f"{n1}x{n2}", sub, dict(cell=f"{n1}x{n2}", origin=oij, seam_score=score, n_atoms=len(sub), n_core=len(core_idx),
                                                      centre_in_cell=int(centre_idx), repair=rep_note, stats=stats, status="PASS")))
                break                       # first acceptable origin for this size
        if passed and first_only: break     # the smaller cell is enough unless both sizes were asked for
    return passed, log


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--centres", default=f"{R}/centres/centres.jsonl"); ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--no-repair", action="store_true")
    ap.add_argument("--out", default="cells", help="output folder under rough_sampling_v1")
    a = ap.parse_args()
    global CELLDIR; CELLDIR = f"{R}/{a.out}"
    os.makedirs(CELLDIR, exist_ok=True)
    C = [json.loads(l) for l in open(a.centres)]
    if a.limit: C = C[:a.limit]
    parents = {}
    rows = []; n_pass = 0; n_fail = 0
    for k, c in enumerate(C):
        pid = c["parent_id"]
        key = (pid, c["frame"], c["source"])
        if key not in parents:
            if c["source"] == "md":
                fr = read(f"{R}/md/{pid}/traj.lammpstrj", index=c["frame"], format="lammps-dump-text")
                fr.set_chemical_symbols(["Au"] * len(fr)); fr.set_pbc((True, True, False))
                p0 = read(f"{R}/parents/{pid}.extxyz")
                for key_ in ("layer", "fixed"): fr.set_array(key_, p0.get_array(key_))
                parents[key] = fr
            else:
                parents[key] = read(f"{R}/parents/{pid}.extxyz")
        parent = parents[key]
        base_id = f"{pid}_f{c['frame']:03d}_a{c['atom']:04d}"
        # test-split centres are cut in BOTH sizes when possible, so the final selection can hold size pairs
        passed, log = extract_one(parent, c["atom"], base_id, repair=not a.no_repair, first_only=(c["split"] != "test"))
        meta = {k_: c[k_] for k_ in ("parent_id", "cls", "split", "frame", "source", "atom", "cn", "novelty_vs_old")}
        if not passed:
            rows.append(dict(cell_id=base_id, centre_id=base_id, **meta, status="FAIL", trials=log)); n_fail += 1
        for size, sub, info in passed:
            cid = f"{base_id}_{size}"; d = f"{CELLDIR}/{cid}"; os.makedirs(d, exist_ok=True)
            write(f"{d}/POSCAR", sub, format="vasp", direct=False, sort=False)
            write(f"{d}/cell.extxyz", sub, format="extxyz")
            rows.append(dict(cell_id=cid, centre_id=base_id, **meta, **info, trials=log)); n_pass += 1
        if (k + 1) % 25 == 0: print(f"  {k+1}/{len(C)}  cells {n_pass}  failed centres {n_fail}", flush=True)
    with open(f"{CELLDIR}/cells_manifest.jsonl", "w") as f:
        for r in rows: f.write(json.dumps(r) + "\n")
    import collections
    fails = collections.Counter(p for r in rows if r["status"] == "FAIL" for t in r["trials"] for p in t["problems"])
    ok = [r for r in rows if r["status"] == "PASS"]
    L = [f"# Extracted cells", "", f"{len(C)} centres tried: {len(C) - n_fail} reconstructed into a periodic electrode that passes every check "
         f"({n_pass} cells, test centres in both sizes where both pass), {n_fail} dropped.",
         "", "Cell sizes: " + str(dict(collections.Counter(r["cell"] for r in ok))),
         "Repair outcomes: " + str(dict(collections.Counter(r["repair"] for r in ok))), "",
         "Most common reasons a trial failed:"] + [f"- {k}: {v}" for k, v in fails.most_common(8)]
    open(f"{CELLDIR}/cells_summary.md", "w").write("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    main()
