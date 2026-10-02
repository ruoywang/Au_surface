#!/bin/bash
# After the 32 MD jobs: completeness check (STOPS if any run is incomplete), per-run summary, then the two
# CANDIDATE-generation steps only (centres, cells). finalize_200 and rough_dft prepare are NOT run here: they are
# run by hand, their output is a candidate list, and nothing is submitted (review of 2026-10-02).
# Usage: nohup rough_sampling_v1/run_after_md.sh > rough_sampling_v1/logs/pipeline.out 2>&1 &
R=/anvil/scratch/x-rywang/Au_Cl/rough_sampling_v1; cd /anvil/scratch/x-rywang/Au_Cl || exit 1
mkdir -p $R/logs; LOG=$R/logs/pipeline.log
say() { echo "$(date '+%Y-%m-%d %H:%M:%S') $*" | tee -a $LOG; }
say "waiting for the MD jobs (md_P*) to leave the queue"
while squeue -u $USER -h -o "%j" | grep -q "^md_P"; do sleep 120; done
say "queue empty of md_ jobs"
# completeness by trajio.run_complete: DONE line for the full protocol AND the dump ends at the last step
$R/pyrun_rs.sh - <<'EOF' 2>&1 | tee -a $LOG
import json, sys; sys.path.insert(0, "/anvil/scratch/x-rywang/Au_Cl/rough_sampling_v1")
import trajio
M = json.load(open("/anvil/scratch/x-rywang/Au_Cl/rough_sampling_v1/parents/parents_manifest.json"))
bad = []
for m in M:
    ok, why = trajio.run_complete(m["parent_id"])
    if not ok: bad.append((m["parent_id"], why))
print(f"MD complete: {len(M) - len(bad)} of {len(M)}" + ("; INCOMPLETE: " + "; ".join(f"{p} ({w})" for p, w in bad) if bad else ""))
sys.exit(1 if bad else 0)
EOF
if [ "${PIPESTATUS[0]}" -ne 0 ]; then say "STOP: incomplete MD runs listed above; no sampling (no fallback in production mode)"; exit 2; fi
# per-run summary: wall time, final temperature, layer changes
$R/pyrun_rs.sh - <<'EOF' 2>&1 | tee -a $LOG
import glob, os, re, json, numpy as np, sys
sys.path.insert(0, "/anvil/scratch/x-rywang/Au_Cl/rough_sampling_v1")
from ase.io import read
from extract_cells import layers_from_z
import trajio
R = "/anvil/scratch/x-rywang/Au_Cl/rough_sampling_v1"; rows = []
for d in sorted(glob.glob(f"{R}/md/P*_s*")):
    pid = os.path.basename(d); log = open(f"{d}/log.lammps").read()
    loops = [float(x) for x in re.findall(r"Loop time of ([\d.]+) on", log)]
    steps = trajio.dump_steps(f"{d}/traj.lammpstrj"); th = trajio.thermo_table(f"{d}/log.lammps")
    p0 = read(f"{R}/parents/{pid}.extxyz"); fl = trajio.frame_by_step(pid, steps[-1])
    l0, l1 = layers_from_z(p0), layers_from_z(fl)
    T300 = [th[int(round(t / trajio.DT_PS))][1] for t in (10, 20, 165, 172, 179, 186, 193, 199) if int(round(t / trajio.DT_PS)) in th]
    rows.append(dict(parent=pid, wall_min=round(sum(loops) / 60, 1), n_frames=len(steps), last_step=steps[-1],
                     T_at_sampled_frames_K=[round(x) for x in T300], layers_initial=np.bincount(l0, minlength=6).tolist(),
                     layers_final=np.bincount(l1, minlength=6).tolist(), n_changed_layer=int((l0 != l1).sum())))
json.dump(rows, open(f"{R}/md/md_summary.json", "w"), indent=1)
for r in rows: print("  ", r["parent"], r["wall_min"], "min", r["n_frames"], "frames", "T@300K-frames", r["T_at_sampled_frames_K"], "layers", r["layers_initial"], "->", r["layers_final"], "changed", r["n_changed_layer"])
EOF
say "step 6: select_centres (production mode, 300 K frames by timestep)"
$R/pyrun_rs.sh $R/select_centres.py --n 1000 --out centres >> $LOG 2>&1 || { say "select_centres FAILED (see log)"; exit 3; }
say "step 7: extract_cells"
$R/pyrun_rs.sh $R/extract_cells.py --centres $R/centres/centres.jsonl --out cells >> $LOG 2>&1 || { say "extract_cells FAILED (see log)"; exit 4; }
say "CANDIDATES READY: centres/ and cells/. finalize_200.py and rough_dft.py prepare are run by hand; nothing is submitted."
