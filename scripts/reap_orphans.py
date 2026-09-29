#!/usr/bin/env python3
"""Re-queue tasks orphaned by a farm job that ended while they were still running.

Gap this closes (found 2026-09-29): `production.refresh()` and the farm's own status sweep both skip farm-launched tasks
(`t.get("launcher") != "farm"`), and the farm evaluates only the tasks IT dispatched. So when a farm SLURM job hits its
walltime in the middle of a VASP run, that task keeps status "running" for ever: no farm will pick it up again (only
"pending" tasks are dispatched) and nothing reports it. The queue then never drains.

This reaper, run periodically, finds every task with status "running" and launcher "farm" whose SLURM job is no longer
active, and either
  * accepts it, if the run actually finished (production.evaluate reads the CPM-ion closure), or
  * resets it to pending so a farm takes it again.
On a reset it attaches a warm start when the SAME geometry has a completed neighbour potential with a usable CHGCAR
(ICHARG=1 with a COPY of that CHGCAR; never ICHARG=11/12) -- the converged CP state is independent of the starting
density, so this only removes the cold-start cost of the re-run. Every reset is logged and recorded on the task
(`reset_count`, `reset_history`) so the re-run is never silently confused with a first attempt.

Usage:  scripts/pyrun.sh scripts/reap_orphans.py [--dry-run]
"""
import os
import shutil
import subprocess
import sys
import time

sys.path.insert(0, "/anvil/scratch/x-rywang/Au_Cl/scripts")
import production as P  # noqa: E402

DRY = "--dry-run" in sys.argv
NEIGHBOUR = {"-5.4071": ["-5.1071", "-4.9071"], "-5.1071": ["-4.9071", "-4.7071"], "-4.9071": ["-4.7071", "-5.1071"],
             "-4.7071": ["-4.9071", "-5.1071"], "-4.4071": ["-4.7071", "-4.9071"]}


def active_job_ids():
    out = subprocess.run(["squeue", "-u", os.environ.get("USER", ""), "-h", "-o", "%i"], capture_output=True, text=True).stdout
    return {line.split("_")[0].strip() for line in out.splitlines() if line.strip()}


def warm_source(t, q):
    """(mu, CHGCAR path) of a completed state of the SAME geometry at a neighbouring potential, else None."""
    for mu in NEIGHBOUR.get(t["mu"], []):
        s = next((x for x in q if x["structure_id"] == t["structure_id"] and x["config"] == t["config"] and x["mu"] == mu
                  and x["status"] in ("completed", "reused")), None)
        if not s: continue
        chg = f"{s['dir'] if s['dir'].startswith('/') else P.ROOT + '/' + s['dir'].split(' ')[0]}/CHGCAR"
        if os.path.exists(chg) and os.path.getsize(chg) > 1e6: return mu, chg
    return None


def main():
    active = active_job_ids()
    with P.queue_lock():
        q = P.load_queue(); n_acc = n_reset = 0
        for t in q:
            if t["status"] != "running" or t.get("launcher") != "farm": continue
            if str(t.get("job_id")) in active: continue                      # its farm is still alive
            lo = f"{t['dir']}/log.out"
            if os.path.exists(lo) and time.time() - os.path.getmtime(lo) < 180: continue   # still being written: leave it
            st, note = P.evaluate(t)
            if st == "completed":
                if DRY: print(f"  would accept {t['task_id']}: {note}"); n_acc += 1; continue
                t["status"], t["result"] = "completed", note
                P.log(f"[reaped-completed] {t['task_id']} (farm {t.get('job_id')} ended) :: {note}")
                n_acc += 1
                continue
            w = warm_source(t, q)
            if DRY:
                print(f"  would reset {t['task_id']} ({note}); warm start from {w[0] if w else 'none available'}"); n_reset += 1; continue
            if w and "ICHARG" not in open(f"{t['dir']}/INCAR").read():
                shutil.copy(w[1], f"{t['dir']}/CHGCAR")
                inc = open(f"{t['dir']}/INCAR").read()
                open(f"{t['dir']}/INCAR", "w").write(inc.replace("IBRION = -1", "ICHARG = 1\nIBRION = -1"))
                t["warm_start"] = f"ICHARG=1 from {w[1]} (added on re-queue after a walltime kill)"; t["warm_start_mu"] = w[0]
            for f in ("log.out", "OUTCAR", "OSZICAR"):                       # keep the killed run out of the way of evaluate()
                p = f"{t['dir']}/{f}"
                if os.path.exists(p): os.rename(p, f"{p}.killed-{t.get('job_id')}")
            t["reset_count"] = t.get("reset_count", 0) + 1
            t.setdefault("reset_history", []).append(dict(at=time.strftime("%Y-%m-%d %H:%M"), farm=t.get("job_id"), node=t.get("node"),
                                                          ran_from=t.get("started_at"), reason=note,
                                                          warm_start=(w[0] if w else None)))
            t.update(status="pending", job_id=None, node=None, started_at=None, launcher=None, elapsed_farm_min=None, result=None)
            P.log(f"[reaped-requeued] {t['task_id']} (farm ended mid-run: {note}); warm start {w[0] if w else 'none'}; attempt #{t['reset_count']+1}")
            n_reset += 1
        if not DRY:
            q.sort(key=lambda x: (x["priority"], x["n_atoms"], x["task_id"]))
            P.save_queue(q)
    print(f"orphans accepted {n_acc}, re-queued {n_reset}")


main()
