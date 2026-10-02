#!/usr/bin/env python3
"""Recalibrate the single-point cost model from the FINISHED rough runs (the two references and, later, the pairs),
instead of the dataset_v1 extrapolation. Writes rough_sampling_v1/dft/cost_model.json, which rough_dft.walltime_minutes
and the budget tables use when present.

Model: minutes = k * (N / 275)^1.5, k fitted per run; the file records every measured point (task, atoms, elapsed,
SCF time, electronic steps, CP rounds, MaxRSS), k_median and k_max over the finished runs, and the walltime policy
(k_max * 1.5 margin). With fewer than 3 points the policy is stated as provisional.

Usage (from Au_Cl/):  rough_sampling_v1/pyrun_rs.sh size_complex_test_v1/recalibrate_cost.py
"""
import glob
import json
import os
import re
import subprocess
import sys
import time

import numpy as np

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
R = f"{ROOT}/rough_sampling_v1"; T = f"{ROOT}/size_complex_test_v1"


def measure(d, job_id=None):
    oc = open(f"{d}/OUTCAR").read()
    el = re.search(r"Elapsed time \(sec\):\s+([\d.]+)", oc); loops = [float(x) for x in re.findall(r"LOOP:\s+cpu time\s+[\d.]+: real time\s+([\d.]+)", oc)]
    n = len(open(f"{d}/POSCAR").read().splitlines()) - 9
    n = int(open(f"{d}/POSCAR").read().splitlines()[6].split()[0])
    rss = ""
    if job_id:
        out = subprocess.run(["sacct", "-j", str(job_id), "-n", "-o", "MaxRSS"], capture_output=True, text=True).stdout.split()
        vals = [float(s[:-1]) * {"K": 1e-6, "M": 1e-3, "G": 1}[s[-1]] for s in out if s and s[-1] in "KMG"]
        rss = round(max(vals), 1) if vals else ""
    return dict(dir=d, n_atoms=n, elapsed_min=round(float(el.group(1)) / 60, 1) if el else None, scf_min=round(sum(loops) / 60, 1), electronic_steps=len(loops),
                cp_closures=len(re.findall(r"CPM-ion", open(f"{d}/log.out").read())), maxrss_GB=rss)


def main():
    pts = []
    for t in json.load(open(f"{R}/dft/queue.json")):
        if t["status"] == "complete": pts.append(dict(task=t["task_id"], seeded=True, **measure(t["dir"], t.get("job_id"))))
    pq = f"{T}/pair_queue.json"
    if os.path.exists(pq):
        for t in json.load(open(pq))["tasks"]:
            if "dir" in t and os.path.exists(f"{t['dir']}/OUTCAR") and "General timing" in open(f"{t['dir']}/OUTCAR").read():
                pts.append(dict(task=t["task_id"], seeded=True, **measure(t["dir"], t.get("job_id"))))
    pts = [p for p in pts if p["elapsed_min"]]
    if not pts: sys.exit("no finished rough run yet")
    ks = np.array([p["elapsed_min"] / (p["n_atoms"] / 275.0) ** 1.5 for p in pts])
    model = dict(updated=time.strftime("%Y-%m-%d %H:%M"), n_points=len(pts), exponent=1.5, n_ref=275.0, k_median_min=float(np.median(ks)), k_max_min=float(ks.max()),
                 walltime_policy="k_max * (N/275)^1.5 * 1.5, rounded up to the hour" + (" (PROVISIONAL: fewer than 3 measured points)" if len(pts) < 3 else ""),
                 margin=1.5, points=pts)
    json.dump(model, open(f"{R}/dft/cost_model.json", "w"), indent=1)
    print(f"{len(pts)} measured runs; k (min at 275 atoms): median {model['k_median_min']:.0f}, max {model['k_max_min']:.0f} "
          f"(dataset_v1 extrapolation was seeded median 399 / max 563, cold median 1011 / max 1346)")
    for p in pts: print(f"  {p['task']:50s} {p['n_atoms']} Au  elapsed {p['elapsed_min']} min  SCF {p['scf_min']} min  {p['electronic_steps']} steps  {p['cp_closures']} closures  MaxRSS {p['maxrss_GB']} GB")


if __name__ == "__main__":
    main()
