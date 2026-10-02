#!/bin/bash
# Chain steps 6-10 once the 32 MD jobs have finished. Submits NO DFT (rough_dft.py prepare only writes inputs).
# Usage: nohup rough_sampling_v1/run_after_md.sh > rough_sampling_v1/logs/pipeline.out 2>&1 &
R=/anvil/scratch/x-rywang/Au_Cl/rough_sampling_v1; cd /anvil/scratch/x-rywang/Au_Cl || exit 1
mkdir -p $R/logs; LOG=$R/logs/pipeline.log
say() { echo "$(date '+%Y-%m-%d %H:%M:%S') $*" | tee -a $LOG; }
say "waiting for the MD jobs (md_P*) to leave the queue"
while squeue -u $USER -h -o "%j" | grep -q "^md_P"; do sleep 120; done
say "queue empty of md_ jobs"
ok=0; bad=""
for d in $R/md/P*_s*/; do p=$(basename $d)
  if grep -q "^DONE $p total_ps 200" $d/stdout.txt 2>/dev/null && [ "$(grep -c '^ITEM: TIMESTEP' $d/traj.lammpstrj)" -ge 200 ]; then ok=$((ok+1)); else bad="$bad $p"; fi
done
say "MD finished cleanly: $ok of 32; incomplete:${bad:- none}"
# per-run summary: wall time, final temperature, how many atoms changed layer
$R/pyrun_rs.sh - <<'EOF' 2>&1 | tee -a $LOG
import glob, os, re, json, numpy as np, sys
sys.path.insert(0, "/anvil/scratch/x-rywang/Au_Cl/rough_sampling_v1")
from ase.io import read
from extract_cells import layers_from_z
R = "/anvil/scratch/x-rywang/Au_Cl/rough_sampling_v1"; rows = []
for d in sorted(glob.glob(f"{R}/md/P*_s*")):
    pid = os.path.basename(d); log = open(f"{d}/log.lammps").read() if os.path.exists(f"{d}/log.lammps") else ""
    loops = [float(x) for x in re.findall(r"Loop time of ([\d.]+) on", log)]
    try:
        p0 = read(f"{R}/parents/{pid}.extxyz"); fl = read(f"{d}/traj.lammpstrj", index=-1, format="lammps-dump-text")
        l0, l1 = layers_from_z(p0), layers_from_z(fl)
        rows.append(dict(parent=pid, wall_min=round(sum(loops) / 60, 1), n_frames=len(read(f"{d}/traj.lammpstrj", index=":", format="lammps-dump-text")),
                         layers_initial=np.bincount(l0, minlength=6).tolist(), layers_final=np.bincount(l1, minlength=6).tolist(),
                         n_changed_layer=int((l0 != l1).sum()), z_max_final=float(fl.positions[:, 2].max() - fl.positions[:, 2].min())))
    except Exception as e:
        rows.append(dict(parent=pid, error=str(e)[:120]))
json.dump(rows, open(f"{R}/md/md_summary.json", "w"), indent=1)
for r in rows: print("  ", r.get("parent"), r.get("wall_min"), "min", r.get("n_frames"), "frames", "layers", r.get("layers_initial"), "->", r.get("layers_final"), "changed", r.get("n_changed_layer"), r.get("error", ""))
EOF
say "step 6: select_centres"; $R/pyrun_rs.sh $R/select_centres.py --n 1000 --frames-per-parent 8 --out centres >> $LOG 2>&1 || { say "select_centres FAILED"; exit 2; }
say "step 7: extract_cells"; $R/pyrun_rs.sh $R/extract_cells.py --centres $R/centres/centres.jsonl --out cells >> $LOG 2>&1 || { say "extract_cells FAILED"; exit 3; }
say "step 8: finalize_200"; $R/pyrun_rs.sh $R/finalize_200.py --cells cells --out rough200 --n 200 --seed 20261002 >> $LOG 2>&1 || { say "finalize_200 FAILED"; exit 4; }
say "step 9: rough_dft prepare (inputs + queue + budget; NOTHING SUBMITTED)"; $R/pyrun_rs.sh $R/rough_dft.py --dftdir dft prepare --manifest rough200/rough200_manifest.json >> $LOG 2>&1 || { say "rough_dft prepare FAILED"; exit 5; }
say "step 10: lineage figures"; $R/pyrun_rs.sh $R/render_rough_lineage.py --cells cells --manifest rough200/rough200_manifest.json --out figures >> $LOG 2>&1 || { say "render FAILED"; exit 6; }
say "PIPELINE DONE"
