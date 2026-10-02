#!/usr/bin/env python3
"""CP-DFT inputs, queue and state machine for the rough states. SUBMITS NOTHING unless `submit --confirm` is given.

prepare   For every state in rough200/rough200_manifest.json write dft/<cell_id>/rough__mu<TARGETMU>/{POSCAR, INCAR,
          KPOINTS, POTCAR, job-run} in the production standard, imported from scripts/production.py so it cannot drift:
          INCAR_SP (K.8: PREC Normal, ENCUT 500, ALGO Fast, EDIFF 1e-7, ISMEAR 1 / 0.2, LSOL ISOL 2 C_MOLAR 1 R_ION 4,
          LCEP NESCHEME 3 FERMICONVERGE 0.01, IDIPOL 3, LVAC single-sided window SOL_Z0 9.801 / SOL_Z1 34.603), Lz 44.603,
          metal bottom at 5.0 A with the layers at 5.0 and 7.4 A fixed, k-mesh ceil(35 A / |a_i|) x 1, JOBRUN 16 MPI x 8 OpenMP.
          Added: a NELECT start guess from the dataset-wide capacitance (recorded in the INCAR comment; the CP loop decides
          the final N_e; the state is labelled by its converged mu_e, never by the guess or TARGETMU).
          Also writes dft/queue.json (all 'pending') and dft/budget.md (node-hours and disk from measured large-cell costs).
status    Refresh every task's state from its files and SLURM:
              pending -> submitted -> running -> scf_converged -> cp_converged -> writing_fields -> complete
                                                                              \\-> partial (ended while writing / a field truncated)
                                                                   \\-> failed (no OUTCAR, CP never closed on target, VASP error)
          Field completeness is checked by counting values against the grid in each file's header (as the dataset export
          did when it found C2's truncated RHOB), so a killed-mid-write run is 'partial', not 'complete'.
submit    --first N [--confirm]: choose N pending tasks (balanced over class, size and |U|; printed) and sbatch them ONLY
          with --confirm. Without --confirm it prints the list and exits. Default N = 8, the plan's first batch.

Usage (from Au_Cl/):  rough_sampling_v1/pyrun_rs.sh rough_sampling_v1/rough_dft.py prepare | status | submit --first 8 [--confirm]
"""
import argparse
import collections
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time

import numpy as np
from ase.io import read

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
R = f"{ROOT}/rough_sampling_v1"
DFT = f"{R}/dft"
QUEUE = f"{DFT}/queue.json"
sys.path.insert(0, f"{ROOT}/scripts")
import production as P  # noqa: E402   INCAR_SP, JOBRUN, write_poscar, kmesh, LZ, SOL_Z0, SOL_Z1, ZMIN, POTCAR_SRC, PARTITION
sys.path.insert(0, R)
import fieldio  # noqa: E402

ZVAL_AU = 11.0
MU0_REF = -4.9071              # U = MU0_REF - TARGETMU (sign of U from a task id)
# 2026-10-02: 'standard' and 'wholenode' are the same nodes (a[250-999]), QoS (part-standard) and billing (CPU=1.0), but
# wholenode had 1096 pending jobs (estimated start 5 days) and standard 16 (estimated start 30 min); a job may list both
ROUGH_PARTITION = "standard,wholenode"
MU_TOL = 0.011                 # same acceptance as production.evaluate: |mu_e - TARGETMU| within FERMICONVERGE
E_PER_UC_CM2 = 1.0 / 1602.18   # 1 uC/cm^2 = 1/1602.18 e/A^2
# measured single-point wall time of the eight > 200-atom reference states (dataset_v1, 231-300 atoms), minutes:
#   cold start  (13 states): median 1011, max 1346        NELECT-seeded (4 states): median 399, max 563
# scaled with (n/275)^1.5 -- an assumption, not a measurement, above 300 atoms (nothing larger has been run here)
COST = dict(seeded_median=399.0, seeded_max=563.0, cold_median=1011.0, cold_max=1346.0, n_ref=275.0, exponent=1.5)
DISK_GB_PER_STATE = 11.0       # du of Island-7-8x8 and Pit-19-8x8 single points with all fields


def scale(n): return (n / COST["n_ref"]) ** COST["exponent"]


def walltime_minutes(n):
    """Requested walltime: the COLD envelope with 20 % margin, because a run killed by the walltime loses everything.
    Recalibrate from the first batch before the remaining 192 are prepared (the plan's reason for running 8 first)."""
    return int(math.ceil((COST["cold_max"] * scale(n) * 1.2) / 60.0)) * 60


def hhmmss(minutes): return f"{minutes // 60:02d}:{minutes % 60:02d}:00"


def charging_constants():
    d = json.load(open(f"{ROOT}/analysis/charging/charging.json"))["decomposition"]["all geometries"]
    return float(d["C_median"]), float(d["U_pzc_median"])


def run_ready(cell_dir):
    at = read(f"{cell_dir}/cell.extxyz")
    pos = at.get_positions(); pos[:, 2] += P.ZMIN - pos[:, 2].min()
    cell = at.get_cell().array.copy(); cell[2] = [0, 0, P.LZ]
    at.set_cell(cell, scale_atoms=False); at.set_positions(pos); at.set_pbc((True, True, True))
    assert np.linalg.det(cell) > 0, "left-handed cell"
    zt = pos[:, 2].max(); assert zt < P.SOL_Z1 - 15.0, f"metal top {zt:.2f} A too high for the fixed solvent window"
    movable = pos[:, 2] > 8.0                               # layers at 5.0 and 7.4 A fixed, as for the (111) family
    assert 0 < movable.sum() < len(at)
    return at, movable


def load_queue(): return json.load(open(QUEUE)) if os.path.exists(QUEUE) else []
def save_queue(q): json.dump(q, open(QUEUE + ".tmp", "w"), indent=1); os.replace(QUEUE + ".tmp", QUEUE)


def prepare(args):
    man = json.load(open(f"{R}/{args.manifest}"))
    C_uF, U_pzc = charging_constants()
    os.makedirs(DFT, exist_ok=True); q = load_queue(); have = {t["task_id"] for t in q}; made = 0
    for s in man["states"]:
        mu = s["TARGETMU_eV"]; task = f"{s['cell_id']}__rough__mu{mu:.4f}"
        if task in have: continue
        at, movable = run_ready(s["cell_dir"])
        d = f"{DFT}/{s['cell_id']}/rough__mu{mu:.4f}"; os.makedirs(d, exist_ok=True)
        P.write_poscar(f"{d}/POSCAR", at, movable)
        kp, n = P.kmesh(at); open(f"{d}/KPOINTS", "w").write(kp)
        # NELECT start guess: sigma = C (U - U_pzc), N_e = N_neutral - sigma A / e
        a, b = at.cell[0][:2], at.cell[1][:2]; A = abs(a[0] * b[1] - a[1] * b[0]); n_neutral = ZVAL_AU * len(at)
        dN = -C_uF * (s["U_V"] - U_pzc) * E_PER_UC_CM2 * A
        incar = P.INCAR_SP.replace("dataset_plan_v1 rev 2 production single point", "rough_sampling_v1 single point (production standard K.8)")
        incar = incar.format(task=task, mu=f"{mu:.4f}", sol_z0=P.SOL_Z0, sol_z1=P.SOL_Z1)
        incar += (f"\nNELECT = {n_neutral + dN:.4f}   # start guess only: N_neutral {n_neutral:.0f} + {dN:+.4f} e from the dataset-wide "
                  f"C = {C_uF:.2f} uF/cm2, U_pzc = {U_pzc:+.4f} V, A = {A:.1f} A^2, U = {s['U_V']:+.2f} V; the CP loop decides the final N_e\n")
        open(f"{d}/INCAR", "w").write(incar)
        shutil.copy(P.POTCAR_SRC, f"{d}/POTCAR")
        wall = walltime_minutes(len(at))
        open(f"{d}/job-run", "w").write(P.JOBRUN.format(task=task, walltime=hhmmss(wall), partition=ROUGH_PARTITION))
        q.append(dict(task_id=task, state_id=s["state_id"], cell_id=s["cell_id"], cls=s["cls"], split=s["split"], cell=s["cell"], pair=s.get("pair"),
                      U_V=s["U_V"], TARGETMU=mu, dir=d, n_atoms=len(at), n_movable=int(movable.sum()), kpoints=f"{n[0]}x{n[1]}x1",
                      nelect_guess=round(n_neutral + dN, 4), walltime_min=wall, status="pending", job_id=None, history=[["pending", time.strftime("%Y-%m-%d %H:%M")]]))
        made += 1
    # sync: pending tasks whose state is no longer in the manifest were replaced before the freeze -> withdrawn (dirs kept)
    live = {f"{s['cell_id']}__rough__mu{s['TARGETMU_eV']:.4f}" for s in man["states"]}; n_wd = 0
    for t in q:
        if t["task_id"] not in live and t["status"] == "pending":
            t["status"] = "withdrawn"; t["note"] = f"replaced in the candidate list before the freeze ({time.strftime('%Y-%m-%d %H:%M')})"; n_wd += 1
    save_queue(q)
    budget([t for t in q if t["status"] != "withdrawn"], C_uF, U_pzc)
    print(f"{made} new tasks written, {n_wd} withdrawn; queue {len(q)} entries, {sum(t['status'] == 'pending' for t in q)} pending. NOTHING SUBMITTED.")


def budget(q, C_uF, U_pzc):
    L = ["# Rough CP-DFT budget (estimate, nothing submitted)", "",
         f"{len(q)} single points. Production standard K.8 from scripts/production.py; one node (16 MPI x 8 OpenMP) per task on `{P.PARTITION}`.",
         f"NELECT start guess from C = {C_uF:.2f} uF/cm2, U_pzc = {U_pzc:+.4f} V (dataset-wide medians; the per-geometry PZC spread is 254 mV, so the guess is off by up to ~0.5 e in an 8x8 cell).", "",
         "Measured cost basis (dataset_v1, 231-300 atoms, minutes per single point): cold start median 1011 / max 1346 (13 states); "
         "NELECT-seeded median 399 / max 563 (4 states, not a controlled comparison). Scaled by (n/275)^1.5 above that range: an ASSUMPTION.", "",
         "| cell | tasks | atoms | node-h if seeded median | if seeded max | if cold median | if cold max | walltime requested |", "|---|---|---|---|---|---|---|---|"]
    tot = collections.defaultdict(float)
    for size in sorted({t["cell"] for t in q}):
        T = [t for t in q if t["cell"] == size]; ns = np.array([t["n_atoms"] for t in T])
        row = {k: float(sum(COST[k] * scale(n) for n in ns)) / 60 for k in ("seeded_median", "seeded_max", "cold_median", "cold_max")}
        for k, v in row.items(): tot[k] += v
        L.append(f"| {size} | {len(T)} | {ns.min()}-{ns.max()} | {row['seeded_median']:.0f} | {row['seeded_max']:.0f} | {row['cold_median']:.0f} | {row['cold_max']:.0f} | "
                 f"{hhmmss(walltime_minutes(int(ns.min())))}-{hhmmss(walltime_minutes(int(ns.max())))} |")
    L.append(f"| **all** | {len(q)} | | **{tot['seeded_median']:.0f}** | {tot['seeded_max']:.0f} | {tot['cold_median']:.0f} | **{tot['cold_max']:.0f}** | |")
    L += ["", f"For scale: dataset_v1 (505 states) cost 757 node-h in total.",
          f"Disk: ~{DISK_GB_PER_STATE:.0f} GB per state with all fields (measured on the 8x8 references) -> ~{DISK_GB_PER_STATE * len(q) / 1000:.1f} TB for {len(q)} states; scratch usage 1.9 TB of 100 TB (2026-10-02).",
          "", "Plan: submit 8 first (`rough_dft.py submit --first 8 --confirm`), recalibrate COST and the walltime from them, then the remaining 192."]
    open(f"{DFT}/budget.md", "w").write("\n".join(L) + "\n"); print("\n".join(L))


# ---------------------------------------------------------------- state machine
# The full-label field set of dataset_v1 (scripts/export_dataset_v1.FIELD_FILES). 'complete' needs ALL of them parsed
# value by value (fieldio.check: grid header, >= nx*ny*nz finite floats, N_ions augmentation blocks for CHGCAR/POT),
# the CHGCAR charge closing on the CP electron count, and CONTCAR repeating POSCAR. There is no shorter list.
FIELD_FILES = ["CHGCAR", "LOCPOT", "PHI", "PHISOLV", "VSOLV", "RHOB", "RHOION", "ELOC", "P", "SVDW", "SION", "SSOLV", "SCAV", "SDIEL", "POT"]
CHARGE_TOL = 0.02      # e: |sum(CHGCAR)/N_grid - N_e(CP)|


def cpm_lines(log_out):
    txt = open(log_out).read() if os.path.exists(log_out) else ""
    ion = re.findall(r"CPM-ion:.*?N_ele=\s*([-\d.]+)\s+mu_e=\s*([-\d.]+)", txt)
    return txt, [(float(a), float(b)) for a, b in ion]


def squeue_mine():
    out = subprocess.run(["squeue", "-u", os.environ["USER"], "-h", "-o", "%i %j %T"], capture_output=True, text=True).stdout
    return {l.split()[1]: (l.split()[0], l.split()[2]) for l in out.splitlines() if l.strip()}


def evaluate(t, active):
    """Return (status, note). `active` = squeue name -> (job id, state)."""
    d = t["dir"]; log_out = f"{d}/log.out"; outcar = f"{d}/OUTCAR"
    if t["task_id"] in active:
        st = active[t["task_id"]][1]
        if st == "PENDING": return "submitted", "in queue"
        live = True
    else:
        live = False
        if t["status"] == "pending": return "pending", ""
    if not os.path.exists(outcar): return ("running", "started, no OUTCAR yet") if live else ("failed", "no OUTCAR")
    txt, ion = cpm_lines(log_out)
    if not ion: return ("running", "first SCF running") if live else ("failed", "no CP round closed")
    ne, mu = ion[-1]
    on_target = abs(mu - t["TARGETMU"]) <= MU_TOL
    oc = open(outcar).read()
    finished = "General timing and accounting" in oc
    if not on_target:
        return ("scf_converged", f"{len(ion)} CP rounds, mu_e {mu:.4f} vs target {t['TARGETMU']:.4f}") if live else ("failed", f"ended with mu_e {mu:.4f} off target after {len(ion)} rounds")
    if live and not finished: return "writing_fields", f"CP closed (N_e {ne:.4f}, mu_e {mu:.4f}); VASP still writing"
    if not finished: return "partial", f"CP closed but VASP ended before its timing summary (killed while writing fields?)"
    qc = {}; bad = []
    for f in FIELD_FILES:
        ok, note, st = fieldio.check(f"{d}/{f}", f); qc[f] = note
        if not ok: bad.append(f"{f} ({note})")
        elif f == "CHGCAR": qc["charge_closure_e"] = st["sum_over_grid"] - ne
    geom_ok, geom_note = fieldio.contcar_matches_poscar(f"{d}/POSCAR", f"{d}/CONTCAR"); qc["contcar_matches_poscar"] = geom_note
    t["qc"] = qc
    if bad: return "partial", f"CP closed, N_e {ne:.4f} mu_e {mu:.4f}; fields failing the value check: " + "; ".join(bad)
    if "charge_closure_e" in qc and abs(qc["charge_closure_e"]) > CHARGE_TOL: return "partial", f"CHGCAR integrates to N_e {qc['charge_closure_e']:+.4f} e (tolerance {CHARGE_TOL})"
    if not geom_ok: return "partial", f"CONTCAR differs from POSCAR ({geom_note})"
    return "complete", f"N_e {ne:.4f} mu_e {mu:.4f}, {len(ion)} CP rounds, 15 fields parsed and finite, charge closure {qc['charge_closure_e']:+.4f} e, geometry unchanged"


def status(args):
    q = load_queue(); active = squeue_mine(); changed = 0
    for t in q:
        if t["status"] == "withdrawn": continue                                   # replaced before the freeze; never evaluated
        if t["status"] in ("complete", "failed") and not args.recheck: continue
        st, note = evaluate(t, active)
        if st != t["status"] or note != t.get("note"):
            t["status"], t["note"] = st, note; t.setdefault("history", []).append([st, time.strftime("%Y-%m-%d %H:%M")]); changed += 1
    save_queue(q)
    c = collections.Counter(t["status"] for t in q)
    print("queue:", dict(c), f"({changed} updated)")
    for t in q:
        if t["status"] not in ("pending", "complete", "withdrawn"): print(f"  {t['status']:15s} {t['task_id']:50s} job={t['job_id']} {t.get('note', '')[:110]}")


def pick_first(q, n):
    """The calibration batch: two states per class (A, B, C, D), one at U > 0 and one at U < 0. Cell sizes alternate
    with class and sign so both 8x8 and 10x8 (both atom-count ranges) appear; within a (class, sign, size) the
    largest |U| is taken: the longest CP walk, i.e. the conservative calibration case. Deterministic and printed."""
    pend = [t for t in q if t["status"] == "pending"]
    sizes = sorted({t["cell"] for t in pend}); out = []
    for k, cls in enumerate("ABCD"):
        for s, sgn in enumerate((+1, -1)):
            want = sizes[(k + s) % len(sizes)] if sizes else None
            c = [t for t in pend if t["cls"] == cls and (t["U_V"] > 0) == (sgn > 0) and t not in out]
            cs = [t for t in c if t["cell"] == want] or c
            if cs: out.append(max(cs, key=lambda t: abs(t["U_V"])))
    for t in sorted(pend, key=lambda t: (t["cls"], -abs(t["U_V"]))):
        if len(out) >= n: break
        if t not in out: out.append(t)
    return out[:n]


def submit(args):
    man = json.load(open(f"{R}/{args.manifest}"))
    if not man.get("frozen"):
        print(f"{args.manifest} is a CANDIDATE list (frozen = false). Nothing can be submitted until it is frozen by approval "
              f"(finalize_200.py freeze). Listing the batch that WOULD be submitted:")
    q = load_queue()
    # the batch is pinned once (dft/first_batch.json) so 'the defined 8' stay the defined 8; a pinned task that is no
    # longer pending (replaced before the freeze) is substituted by the rule and the substitution is printed
    pin = f"{DFT}/first_batch.json"
    if os.path.exists(pin):
        pinned = json.load(open(pin)); byid = {t["task_id"]: t for t in q}
        chosen = [byid[i] for i in pinned if i in byid and byid[i]["status"] == "pending"]
        lost = [i for i in pinned if i not in byid or byid[i]["status"] == "withdrawn"]
        if lost:
            print(f"pinned tasks no longer available: {lost}")
            # a substitute must fill the SAME slot (class, sign of U, cell size) with the largest |U| available; the
            # earlier version took the head of pick_first() and put a class-A task into a class-C slot (job 21015731)
            for i in lost:
                m = re.match(r"P([ABCD])\d_s\d+_f\d+_a\d+_(\d+x\d+)__rough__mu(-?[\d.]+)", i); cls, size, mu = m.group(1), m.group(2), float(m.group(3))
                sgn = (MU0_REF - mu) > 0
                c = [t for t in q if t["status"] == "pending" and t["cls"] == cls and (t["U_V"] > 0) == sgn and t["cell"] == size and t not in chosen]
                c = c or [t for t in q if t["status"] == "pending" and t["cls"] == cls and (t["U_V"] > 0) == sgn and t not in chosen]
                if c and len(chosen) < args.first:
                    t = max(c, key=lambda t: abs(t["U_V"])); chosen.append(t); print(f"  slot {cls}/{'U>0' if sgn else 'U<0'}/{size} -> {t['task_id']} ({t['cell']})")
    else:
        chosen = pick_first(q, args.first); json.dump([t["task_id"] for t in chosen], open(pin, "w"), indent=1); print(f"first batch pinned to {pin}")
    print(f"{'WOULD submit' if not args.confirm else 'Submitting'} {len(chosen)} tasks:")
    for t in chosen: print(f"  {t['task_id']:50s} {t['cls']} {t['cell']} U={t['U_V']:+.2f} n={t['n_atoms']} walltime {hhmmss(t['walltime_min'])}")
    if not args.confirm or not man.get("frozen"):
        print("Nothing submitted" + (" (no --confirm)." if man.get("frozen") else " (candidate list not frozen).")); return
    for t in chosen:
        r = subprocess.run(["sbatch", "job-run"], cwd=t["dir"], capture_output=True, text=True)
        m = re.search(r"Submitted batch job (\d+)", r.stdout)
        if m:
            t["job_id"], t["status"] = int(m.group(1)), "submitted"; t.setdefault("history", []).append(["submitted", time.strftime("%Y-%m-%d %H:%M")])
            print(f"  -> job {t['job_id']}")
        else: print(f"  sbatch error: {r.stderr.strip()[:200]}"); break
    save_queue(q)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--dftdir", default="dft", help="folder under rough_sampling_v1 (test runs use another)")
    sub = ap.add_subparsers(dest="cmd")
    p1 = sub.add_parser("prepare"); p1.add_argument("--manifest", default="rough200/rough200_manifest.json")
    p2 = sub.add_parser("status"); p2.add_argument("--recheck", action="store_true")
    p3 = sub.add_parser("submit"); p3.add_argument("--first", type=int, default=8); p3.add_argument("--confirm", action="store_true")
    p3.add_argument("--manifest", default="rough200/rough200_manifest.json")
    a = ap.parse_args()
    DFT = f"{R}/{a.dftdir}"; QUEUE = f"{DFT}/queue.json"
    {"prepare": prepare, "status": status, "submit": submit}.get(a.cmd, lambda _: ap.print_help())(a)
