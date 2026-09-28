#!/usr/bin/env python3
"""Dataset production driver (frozen plan dataset_plan_v1.csv rev 2), 2026-09-27.

Rules implemented here, from the user's decisions:
  * Tasks are generated ONLY from explicit rows of the frozen CSV (never by scanning the structure library);
    files with role retired_mislabeled are never touched.
  * Production numerical standard (parameter_map.md K.8), unchanged: PREC=Normal, ALGO=Fast, NPAR=16, 16 MPI x 8 OpenMP,
    FERMICONVERGE=0.01, EDIFF=1e-7, ENCUT=500, ISMEAR=1/SIGMA=0.2, default mixing; single-sided window SOL_Z0=9.801,
    SOL_Z1=34.603, Lz=44.603, DIPOL 0.5 0.5 0.5; k-mesh N_i = ceil(35 A / |a_i|) (Gamma-centred, 1 along z).
  * Order: main-row relaxations at TARGETMU=-4.9071 first; ideal single points of main and reference rows follow and
    are interleaved by the feeder; children (perturbations/collective/paths) of a parent are generated as soon as ITS
    relaxation is accepted (scripts/production_children.py), not when the whole batch is done.
  * Execution: ONE multi-node 'farm' job (wholenode) queues once; inside it, production.py farm assigns pending tasks to
    free nodes back-to-back (one VASP instance per node, 16 MPI x 8 OpenMP) until the queue is empty or the walltime nears.
    The single-node feeder remains available as a fallback.
  * A relaxation is accepted only if VASP reports 'reached required accuracy' (EDIFFG) and the last CP round closed;
    timed-out / unconverged runs are recorded as failed and are NOT used as relaxed references.
  * Reused states (same production standard) are recorded, not recomputed.

Layout: 05_production/<structure_id>/<config>__mu<value>/{POSCAR,INCAR,KPOINTS,POTCAR,job-run}
        05_production/queue.json (task list + status), 05_production/feeder.log

Usage (from Au_Cl/, via scripts/pyrun.sh):
  production.py prepare-batchA      # write inputs + queue entries for Batch A (idempotent: skips existing tasks)
  production.py feed [--once]       # submit pending tasks up to the limit, refresh statuses (loop unless --once)
  production.py status              # summary
"""
import csv
import fcntl
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time

import numpy as np
from ase.io import read, write

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
LIB = f"{ROOT}/03_pilot/all_defect_structures"
PROD = f"{ROOT}/05_production"
QUEUE = f"{PROD}/queue.json"
LOG = f"{PROD}/feeder.log"
POTCAR_SRC = f"{ROOT}/03_pilot/Flat16x1_muref/POTCAR"
LZ, SOL_Z0, SOL_Z1 = 44.603, 9.801, 34.603
ZMIN = 5.0
PARTITION = "wholenode"      # 2026-09-27: switched from highmem (MaxJobsPU=2/MaxSubmitPU=4, 4x billing). Measured need:
                             # 22-23 GB per 72-atom production run (sacct MaxRSS of the hydra step) vs 257 GB on a standard node.
MAX_ACTIVE = 8               # our submitted+running production jobs across partitions (wholenode allows 64 running per user)
MU_RELAX = "-4.9071"

# library file per structure_id (explicit; nothing is scanned)
LIBFILE = {
    "T-4x4": "T.poscar", "Flat-16x1": "Flat-16x1.poscar", "Flat-8x2": "Flat-8x2.poscar",
    "V1": "V1.poscar", "A1-fcc": "A1_fcc.poscar", "A1-hcp": "A1_hcp.poscar", "V2": "V2.poscar", "V3": "V3.poscar", "A3": "A3.poscar",
    "Step-8x2": "Step-8x2.poscar", "Step-16x2": "Step-16x2.poscar", "Step-8x1": "Step-8x1.poscar", "Step-16x1": "Step-16x1.poscar",
    "Step-24x1": "Step-24x1.poscar", "Au211": "Au211.poscar", "Au221": "Au221.poscar", "Au332": "Au332.poscar", "Au554": "Au554.poscar",
    "Step-8x2_edge-vacancy_plus_foot-adatom": "Step-8x2_edge-vacancy_plus_foot-adatom.poscar",
    "R1-hcp-terminated": "R1-hcp-terminated.poscar", "R2-stripe-wall": "R2-stripe-wall.poscar",
    "Island-7-compact": "Island-7-6x6.poscar", "Pit-7-compact": "V7.poscar",
    "Kink-edge1": "Kink-edge1.poscar", "Kink-edge2": "Kink-edge2.poscar", "Island-7-elongated": "Island-7-elongated.poscar",
    "Pit-7-trench": "Pit-7-trench.poscar", "C1-island-near-step": "C1-island-near-step.poscar", "C2-island+pit": "C2-island+pit.poscar",
    "Island-19-8x8": "Island-19-8x8.poscar", "Island-7-8x8": "Island-7-8x8.poscar", "Pit-19-8x8": "Pit-19-8x8.poscar", "Pit-7-8x8": "Pit-7-8x8.poscar",
}
# states already computed in the production standard (structure_id -> {mu_slot: run dir})
REUSED = {
    "Flat-16x1": {"-4.9071": "03_pilot/Flat16x1_muref", "-5.1071": "03_pilot/Flat16x1_dUp02 (reached mu=-5.098393)"},
    "Step-8x1": {"-4.9071": "03_pilot/Step-8x1_muref_fastcfg"},
    "Step-16x1": {"-4.9071": "03_pilot/Step-16x1_muref_fastcfg", "-5.1071": "03_pilot/Step-16x1_dUp02_fastcfg"},
}
# Batch B/C rows (73-151 atoms; > 200-atom references). Rows whose structure file exists are enqueued at LOWER priority
# than Batch A (relax 40, main ideal 50, children 55, references 60, > 200-atom references 70) so farms take them once
# Batch A is drained. Rows still to be built (kinks, elongated island, trench pit, C1, C2) are skipped until their file exists.
BATCH_B = ["Step-24x1", "Step-16x2", "Island-7-compact", "Pit-7-compact", "Kink-edge1", "Kink-edge2", "Island-7-elongated",
           "Pit-7-trench", "C1-island-near-step", "C2-island+pit"]
BATCH_C = ["Island-19-8x8", "Island-7-8x8", "Pit-19-8x8", "Pit-7-8x8"]
BATCH_A = ["T-4x4", "V1", "A1-fcc", "A1-hcp", "Au221", "Au211", "Step-8x2", "R1-hcp-terminated", "R2-stripe-wall",
           "Flat-16x1", "Flat-8x2", "Step-8x1", "Step-16x1", "V2", "V3", "A3", "Au332", "Au554", "Step-8x2_edge-vacancy_plus_foot-adatom"]

INCAR_SP = """# {task}: dataset_plan_v1 rev 2 production single point (standard K.8). TARGETMU is the internal mu0 reference.
SYSTEM = {task}

NPAR  = 16

PREC   = Normal
ENCUT  = 500
ISPIN  = 1
ISYM   = 0
LREAL  = Auto
ALGO   = Fast
EDIFF  = 1E-7
NELM   = 200
ISMEAR = 1
SIGMA  = 0.2
IBRION = -1
NSW    = 0
LCHARG = .TRUE.
LWAVE  = .FALSE.
LVHAR  = .TRUE.

LSOL    = .TRUE.
ISOL    = 2
C_MOLAR = 1.0
R_ION   = 4.0

LCEP          = .TRUE.
NESCHEME      = 3
TARGETMU      = {mu}
FERMICONVERGE = 0.01
CAP_MAX       = 2.0
NEADJUST      = 1

IDIPOL = 3
LDIPOL = .TRUE.
DIPOL  = 0.5 0.5 0.5

LVAC      = .TRUE.
SOL_Z0    = {sol_z0}
SOL_Z1    = {sol_z1}
SOL_SIGMA = 0.8
"""
INCAR_RELAX = INCAR_SP.replace("production single point", "CP relaxation (movable atoms free, bottom layers fixed)").replace(
    "IBRION = -1\nNSW    = 0", "IBRION = 2\nNSW    = 80\nISIF   = 2\nEDIFFG = -0.02\nPOTIM  = 0.3")

JOBRUN = """#!/bin/bash
#SBATCH --job-name={task}
#SBATCH -o myjob.o%j
#SBATCH -e myjob.e%j
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=8
#SBATCH --time={walltime}
#SBATCH --partition={partition}
#SBATCH --account=CHE190065

module purge
module load intel/19.0.5.281 impi/2019.5.281
module load intel-mkl hdf5 python

ulimit -s unlimited
export OMP_STACKSIZE=512m
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
PW=/anvil/projects/x-che190065/rywang/CEP-HALF/bin
mpirun -np $SLURM_NTASKS $PW/vasp_std > log.out
"""


def log(msg):
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    print(line)
    os.makedirs(PROD, exist_ok=True)
    open(LOG, "a").write(line + "\n")


def load_plan():
    rows = {r["structure_id"]: r for r in csv.DictReader(open(f"{ROOT}/dataset_plan_v1.csv"))}
    for r in rows.values():
        for k in ("n_atoms", "n_ref", "n_relax", "n_perturb", "n_collective", "n_path", "n_reusable_states"):
            r[k] = int(r[k])
        r["mu"] = [m for m in r["mu_list"].split(";") if m]
    return rows


class queue_lock:
    """exclusive lock for read-modify-write of queue.json shared by several farm jobs on DIFFERENT nodes.
    flock() proved ineffective across nodes on this filesystem (2026-09-27: two farms launched the same task
    in the same second), so the lock is an atomic mkdir() of a lock directory, with a stale-lock timeout."""
    LOCKDIR = f"{PROD}/queue.lockdir"
    STALE_S = 900

    def __enter__(self):
        os.makedirs(PROD, exist_ok=True)
        t0 = time.time()
        while True:
            try:
                os.mkdir(self.LOCKDIR)
                open(f"{self.LOCKDIR}/owner", "w").write(f"{os.uname()[1]} pid {os.getpid()} {time.strftime('%F %T')}")
                return self
            except FileExistsError:
                try:
                    age = time.time() - os.stat(self.LOCKDIR).st_mtime
                except FileNotFoundError:
                    continue
                if age > self.STALE_S:
                    log(f"[lock] breaking stale queue lock (age {age:.0f} s)")
                    shutil.rmtree(self.LOCKDIR, ignore_errors=True); continue
                if time.time() - t0 > 2 * self.STALE_S:
                    raise RuntimeError("queue lock not acquired in 30 min")
                time.sleep(1.0 + 2.0 * (os.getpid() % 5) / 5.0)

    def __exit__(self, *a):
        shutil.rmtree(self.LOCKDIR, ignore_errors=True)


def load_queue():
    return json.load(open(QUEUE)) if os.path.exists(QUEUE) else []


def save_queue(q):
    os.makedirs(PROD, exist_ok=True)
    tmp = f"{QUEUE}.{os.uname()[1]}.{os.getpid()}.tmp"        # unique per writer: a shared .tmp name raced between farms on 2026-09-27
    json.dump(q, open(tmp, "w"), indent=1); os.replace(tmp, QUEUE)


def run_ready_atoms(structure_id, row, template_path=None):
    """library template -> run-ready Atoms: Lz, metal bottom at ZMIN, movable/fixed flags."""
    at = read(template_path or f"{LIB}/{LIBFILE[structure_id]}")
    pos = at.get_positions()
    pos[:, 2] += ZMIN - pos[:, 2].min()
    cell = at.get_cell().array.copy(); cell[2] = [0, 0, LZ]
    at.set_cell(cell, scale_atoms=False); at.set_positions(pos)
    assert np.linalg.det(cell) > 0, "left-handed cell"
    zt = pos[:, 2].max()
    assert zt < SOL_Z1 - 15.0, f"metal top {zt:.2f} too high for the fixed window"
    if row["family"] == "vicinal step face":
        movable = pos[:, 2] > zt - 4.8
    else:
        movable = pos[:, 2] > 8.0                       # (111) family: layers at 5.0 and 7.4 fixed
    assert 0 < movable.sum() < len(at)
    return at, movable


def write_poscar(path, at, movable):
    lines = ["Au", "1.0"]
    for v in at.get_cell().array:
        lines.append("  %.16f  %.16f  %.16f" % tuple(v))
    lines += [" Au", f"  {len(at)}", "Selective dynamics", "Cartesian"]
    for p, m in zip(at.get_positions(), movable):
        f = "T" if m else "F"
        lines.append("  %.16f  %.16f  %.16f   %s   %s   %s" % (*p, f, f, f))
    open(path, "w").write("\n".join(lines) + "\n")


def kmesh(at):
    n = [max(1, math.ceil(35.0 / np.linalg.norm(at.get_cell().array[i]))) for i in range(2)]
    return f"Automatic mesh\n0\nGamma\n{n[0]} {n[1]} 1\n0 0 0\n", n


def make_task(q, structure_id, row, config, mu, kind, at, movable, priority, note=""):
    task = f"{structure_id}__{config}__mu{mu}"
    if any(t["task_id"] == task for t in q):
        return None
    d = f"{PROD}/{structure_id}/{config}__mu{mu}"
    os.makedirs(d, exist_ok=True)
    write_poscar(f"{d}/POSCAR", at, movable)
    kp, n = kmesh(at)
    open(f"{d}/KPOINTS", "w").write(kp)
    tmpl = INCAR_RELAX if kind == "relax" else INCAR_SP
    open(f"{d}/INCAR", "w").write(tmpl.format(task=task, mu=mu, sol_z0=SOL_Z0, sol_z1=SOL_Z1))
    shutil.copy(POTCAR_SRC, f"{d}/POTCAR")
    wall = "24:00:00" if kind == "relax" else "04:00:00"
    open(f"{d}/job-run", "w").write(JOBRUN.format(task=task, walltime=wall, partition=PARTITION))
    t = dict(task_id=task, structure_id=structure_id, family=row["family"], tier=row["tier"], config=config, mu=mu, kind=kind,
             dir=d, n_atoms=len(at), n_movable=int(movable.sum()), kpoints=f"{n[0]}x{n[1]}x1", priority=priority,
             status="pending", job_id=None, note=note)
    q.append(t)
    return t


def prepare_batch(batch, prio_relax, prio_main, prio_ref):
    rows = load_plan(); q = load_queue(); made = 0
    for sid in batch:
        row = rows[sid]
        if sid not in LIBFILE or not os.path.exists(f"{LIB}/{LIBFILE[sid]}"):
            log(f"[prepare] {sid}: structure file not available yet -> skipped (will be added when built)"); continue
        at, movable = run_ready_atoms(sid, row)
        assert len(at) == row["n_atoms"], (sid, len(at), row["n_atoms"])
        if row["n_relax"] == 1:
            if make_task(q, sid, row, "relax", MU_RELAX, "relax", at, movable, priority=prio_relax): made += 1
        for mu in row["mu"]:
            if mu in REUSED.get(sid, {}): continue
            prio = prio_main if row["tier"] == "main" else prio_ref
            if make_task(q, sid, row, "ideal", mu, "sp", at, movable, priority=prio): made += 1
    q.sort(key=lambda t: (t["priority"], t["n_atoms"], t["task_id"]))
    save_queue(q)
    log(f"[prepare] {made} new tasks written; queue now {len(q)} entries ({sum(t['status']=='pending' for t in q)} pending)")


def prepare_batchA():
    rows = load_plan(); q = load_queue(); made = 0
    for sid in BATCH_A:
        row = rows[sid]
        if sid not in LIBFILE or not os.path.exists(f"{LIB}/{LIBFILE[sid]}"):
            log(f"[prepare] {sid}: structure file not available yet -> skipped (will be added when built)"); continue
        at, movable = run_ready_atoms(sid, row)
        assert len(at) == row["n_atoms"], (sid, len(at), row["n_atoms"])
        # relaxation (main rows)
        if row["n_relax"] == 1:
            if make_task(q, sid, row, "relax", MU_RELAX, "relax", at, movable, priority=10): made += 1
        # ideal single points, skipping reused slots
        for i, mu in enumerate(row["mu"]):
            if mu in REUSED.get(sid, {}):
                task = f"{sid}__ideal__mu{mu}"
                if not any(t["task_id"] == task for t in q):
                    q.append(dict(task_id=task, structure_id=sid, family=row["family"], tier=row["tier"], config="ideal", mu=mu, kind="sp",
                                  dir=REUSED[sid][mu], n_atoms=row["n_atoms"], priority=0, status="reused", job_id=None,
                                  note="computed earlier in the production standard; not resubmitted"))
                continue
            prio = 20 if row["tier"] == "main" else 30
            if make_task(q, sid, row, "ideal", mu, "sp", at, movable, priority=prio): made += 1
    q.sort(key=lambda t: (t["priority"], t["n_atoms"], t["task_id"]))
    save_queue(q)
    log(f"[prepare] Batch A: {made} new tasks written; queue now {len(q)} entries "
        f"({sum(t['status']=='pending' for t in q)} pending, {sum(t['status']=='reused' for t in q)} reused)")


def squeue_mine():
    """our production jobs (any partition): job id, name, state -- only names that are queue task ids"""
    out = subprocess.run(["squeue", "-u", os.environ["USER"], "-h", "-o", "%i %j %T"], capture_output=True, text=True).stdout
    ids = {t["task_id"] for t in load_queue()}
    return [l.split() for l in out.splitlines() if l.strip() and l.split()[1] in ids]


def sacct_state(job_id):
    out = subprocess.run(["sacct", "-j", str(job_id), "-n", "-X", "-o", "State,Elapsed"], capture_output=True, text=True).stdout.split()
    return (out[0], out[1]) if len(out) >= 2 else (None, None)


def evaluate(t):
    """final status of an ended task from its files."""
    d = t["dir"]; log_out = f"{d}/log.out"; outcar = f"{d}/OUTCAR"
    if not os.path.exists(outcar):
        return "failed", "no OUTCAR"
    txt = open(log_out).read() if os.path.exists(log_out) else ""
    m = re.findall(r"CPM-ion:.*?N_ele=\s*([-\d.]+)\s+mu_e=\s*([-\d.]+)", txt)
    if not m:
        return "failed", "no CPM-ion closure"
    ne, mu = float(m[-1][0]), float(m[-1][1])
    if abs(mu - float(t["mu"])) > 0.011:
        return "failed", f"mu_e {mu:.4f} outside FERMICONVERGE of target {t['mu']}"
    if t["kind"] == "relax":
        oc = open(outcar).read()
        if "reached required accuracy" not in oc:
            return "failed", f"relaxation did not reach EDIFFG (ionic steps: {oc.count('FORCES:') or 'n/a'})"
        nsteps = len(re.findall(r"CPM-ion:", txt))
        return "completed", f"relaxed: {nsteps} ionic steps, final N_e={ne:.6f} mu_e={mu:.6f}"
    return "completed", f"N_e={ne:.6f} mu_e={mu:.6f}"


def feed(once=False):
    while True:
        q = load_queue(); mine = squeue_mine(); active_ids = {j[0] for j in mine}
        # refresh statuses of submitted tasks
        for t in q:
            if t["status"] in ("submitted", "running") and t["job_id"]:
                if str(t["job_id"]) in active_ids:
                    st = [j[2] for j in mine if j[0] == str(t["job_id"])][0]
                    t["status"] = "running" if st == "RUNNING" else "submitted"
                else:
                    state, elapsed = sacct_state(t["job_id"])
                    status, note = evaluate(t)
                    t["status"], t["slurm_state"], t["elapsed"], t["result"] = status, state, elapsed, note
                    log(f"[{status}] {t['task_id']} slurm={state} elapsed={elapsed} :: {note}")
                    if status == "completed" and t["kind"] == "relax":
                        try:
                            subprocess.run([f"{ROOT}/scripts/pyrun.sh", f"{ROOT}/scripts/production_children.py", t["structure_id"]], check=False, timeout=600)
                        except Exception as e:
                            log(f"[children] {t['structure_id']}: generator error {e}")
        # submit
        n_in = len(mine)
        for t in sorted([x for x in q if x["status"] == "pending"], key=lambda x: (x["priority"], x["n_atoms"], x["task_id"])):
            if n_in >= MAX_ACTIVE: break
            r = subprocess.run(["sbatch", "job-run"], cwd=t["dir"], capture_output=True, text=True)
            mj = re.search(r"Submitted batch job (\d+)", r.stdout)
            if mj:
                t["job_id"], t["status"], t["submitted_at"] = int(mj.group(1)), "submitted", time.strftime("%Y-%m-%d %H:%M")
                n_in += 1; log(f"[submit] {t['task_id']} -> job {t['job_id']}")
            else:
                log(f"[submit-error] {t['task_id']}: {r.stderr.strip()[:200]}"); break
        save_queue(q)
        if once: break
        if not any(t["status"] in ("pending", "submitted", "running") for t in q):
            log("[feed] nothing pending or active -- feeder exits"); break
        time.sleep(300)


def retarget(partition):
    """rewrite job-run of every PENDING task to the given partition (idempotent)."""
    q = load_queue(); n = 0
    for t in q:
        if t["status"] == "pending":
            jr = f"{t['dir']}/job-run"; txt = open(jr).read()
            new = re.sub(r"#SBATCH --partition=\S+", f"#SBATCH --partition={partition}", txt)
            if new != txt: open(jr, "w").write(new); n += 1
    log(f"[retarget] {n} pending job-run files now use partition={partition}")


def expected_minutes(t):
    n = t["n_atoms"]; c = 0.55 * n * max(1.0, n / 72.0) ** 0.5
    return c * (10.0 if t["kind"] == "relax" else 1.5)


def farm():
    """In-allocation task farm: run inside a multi-node SLURM job; one VASP instance per node."""
    jobid = os.environ["SLURM_JOB_ID"]
    nodes = subprocess.run(["scontrol", "show", "hostnames", os.environ["SLURM_JOB_NODELIST"]], capture_output=True, text=True).stdout.split()
    end = subprocess.run(["squeue", "-j", jobid, "-h", "-o", "%e"], capture_output=True, text=True).stdout.strip()
    end_ts = time.mktime(time.strptime(end, "%Y-%m-%dT%H:%M:%S")) if end and end != "N/A" else time.time() + 4 * 86400
    log(f"[farm {jobid}] {len(nodes)} nodes {nodes}; allocation ends {end}")
    slots = {n: None for n in nodes}          # node -> (task_id, Popen, t0)
    pw = "/anvil/projects/x-che190065/rywang/CEP-HALF/bin/vasp_std"
    env = dict(os.environ, OMP_NUM_THREADS="8", OMP_STACKSIZE="512m", I_MPI_JOB_RESPECT_PROCESS_PLACEMENT="0", I_MPI_PIN_DOMAIN="omp")
    while True:
        children_for = []
        with queue_lock():
            q = load_queue(); byid = {t["task_id"]: t for t in q}
            # 1) finished farm instances (only the ones this farm launched)
            for node, run in list(slots.items()):
                if run and run[1].poll() is not None:
                    t = byid[run[0]]; status, note = evaluate(t)
                    t["status"], t["result"], t["elapsed_farm_min"], t["node"] = status, note, round((time.time() - run[2]) / 60, 1), node
                    log(f"[{status}] {t['task_id']} on {node} in {t['elapsed_farm_min']} min :: {note}")
                    slots[node] = None
                    if status == "completed" and t["kind"] == "relax": children_for.append(t["structure_id"])
            # 2) tasks submitted as individual SLURM jobs: finished -> evaluate; still PENDING -> adopt into the farm
            mine = squeue_mine(); active = {j[0] for j in mine}
            for j in mine:
                if j[2] == "PENDING" and j[1] in byid and byid[j[1]].get("launcher") != "farm" and byid[j[1]]["status"] == "submitted":
                    subprocess.run(["scancel", j[0]], check=False)
                    byid[j[1]].update(status="pending", job_id=None, note=(byid[j[1]].get("note", "") + f" | single-node job {j[0]} cancelled while PENDING, adopted by farm {jobid}").strip(" |"))
                    log(f"[adopt] {j[1]}: cancelled pending job {j[0]}, will run inside a farm")
                    active.discard(j[0])
            for t in q:
                if t["status"] in ("submitted", "running") and t.get("job_id") and t.get("launcher") != "farm" and str(t["job_id"]) not in active:
                    state, elapsed = sacct_state(t["job_id"]); status, note = evaluate(t)
                    t["status"], t["slurm_state"], t["elapsed"], t["result"] = status, state, elapsed, note
                    log(f"[{status}] {t['task_id']} slurm={state} elapsed={elapsed} :: {note}")
                    if status == "completed" and t["kind"] == "relax": children_for.append(t["structure_id"])
            # 3) launch on free nodes
            remaining_min = (end_ts - time.time()) / 60
            pending = sorted([x for x in q if x["status"] == "pending"], key=lambda x: (x["priority"], x["n_atoms"], x["task_id"]))
            for node in [n for n, r in slots.items() if r is None]:
                cand = next((x for x in pending if expected_minutes(x) * 1.3 + 10 < remaining_min), None)
                if cand is None: break
                pending.remove(cand)
                lo = f"{cand['dir']}/log.out"
                if os.path.exists(lo) and time.time() - os.path.getmtime(lo) < 180:
                    log(f"[skip] {cand['task_id']}: log.out modified {time.time()-os.path.getmtime(lo):.0f} s ago -- another instance may be running; left pending")
                    continue
                cmd = ["mpirun", "-np", "16", "-ppn", "16", "-hosts", node, pw]
                proc = subprocess.Popen(cmd, cwd=cand["dir"], stdout=open(f"{cand['dir']}/log.out", "w"), stderr=open(f"{cand['dir']}/farm.err", "w"), env=env)
                slots[node] = (cand["task_id"], proc, time.time())
                cand.update(status="running", job_id=int(jobid), node=node, started_at=time.strftime("%Y-%m-%d %H:%M"), launcher="farm")
                log(f"[launch] {cand['task_id']} -> {node} (farm {jobid}, expected ~{expected_minutes(cand):.0f} min, {remaining_min/60:.1f} h left)")
            save_queue(q)
            n_pending = len([x for x in q if x["status"] == "pending"])
        for sid in children_for:
            r = subprocess.run([f"{ROOT}/scripts/pyrun.sh", f"{ROOT}/scripts/production_children.py", sid], capture_output=True, text=True, timeout=900)
            if r.returncode != 0:
                log(f"[children-error] {sid}: rc={r.returncode} :: {(r.stderr.strip().splitlines() or ['no stderr'])[-1][:200]}")
        busy = any(r is not None for r in slots.values())
        if not busy and n_pending == 0:
            log(f"[farm {jobid}] queue drained -- exiting"); break
        if not busy and remaining_min < 30:
            log(f"[farm {jobid}] walltime nearly over, nothing launchable -- exiting"); break
        time.sleep(60)


def refresh():
    """refresh statuses of tasks running as individual SLURM jobs and generate children of accepted relaxations; submits nothing."""
    children_for = []
    with queue_lock():
        q = load_queue(); mine = squeue_mine(); active = {j[0] for j in mine}
        for t in q:
            if t["status"] in ("submitted", "running") and t.get("job_id") and t.get("launcher") != "farm" and str(t["job_id"]) not in active:
                state, elapsed = sacct_state(t["job_id"]); st, note = evaluate(t)
                t["status"], t["slurm_state"], t["elapsed"], t["result"] = st, state, elapsed, note
                log(f"[{st}] {t['task_id']} slurm={state} elapsed={elapsed} :: {note}")
                if st == "completed" and t["kind"] == "relax": children_for.append(t["structure_id"])
        save_queue(q)
    for sid in children_for:
        subprocess.run([f"{ROOT}/scripts/pyrun.sh", f"{ROOT}/scripts/production_children.py", sid], check=False, timeout=900)


def status():
    q = load_queue()
    from collections import Counter
    c = Counter(t["status"] for t in q)
    print("queue:", dict(c))
    for t in q:
        if t["status"] in ("submitted", "running", "failed"):
            print(f"  {t['status']:9s} {t['task_id']:55s} job={t['job_id']} {t.get('result','')}")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "status"
    if cmd == "prepare-batchA": prepare_batchA()
    elif cmd == "feed": feed(once="--once" in sys.argv)
    elif cmd == "retarget": retarget(sys.argv[2] if len(sys.argv) > 2 else PARTITION)
    elif cmd == "farm": farm()
    elif cmd == "prepare-batchB": prepare_batch(BATCH_B, 40, 50, 60)
    elif cmd == "prepare-batchC": prepare_batch(BATCH_C, 70, 70, 70)
    elif cmd == "refresh": refresh()
    else: status()
