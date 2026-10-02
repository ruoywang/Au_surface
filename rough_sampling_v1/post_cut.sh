#!/bin/bash
# Run BY HAND after the 8 cut chunks have finished: merge -> candidate list -> DFT inputs + budget -> figures -> report.
# Everything it produces is a CANDIDATE deliverable; nothing is frozen and nothing is submitted.
R=/anvil/scratch/x-rywang/Au_Cl/rough_sampling_v1; cd /anvil/scratch/x-rywang/Au_Cl || exit 1
LOG=$R/logs/pipeline.log; say() { echo "$(date '+%Y-%m-%d %H:%M:%S') $*" | tee -a $LOG; }
n=$(grep -l "^chunk " $R/logs/extract_part*.log 2>/dev/null | wc -l)
[ "$n" -eq 8 ] || { say "STOP: only $n of 8 chunks finished"; exit 2; }
say "merge chunks"; $R/pyrun_rs.sh $R/extract_cells.py --centres $R/centres/centres.jsonl --out cells --merge >> $LOG 2>&1 || { say "merge FAILED"; exit 3; }
say "step 8: finalize_200 (CANDIDATE list, frozen=false)"; $R/pyrun_rs.sh $R/finalize_200.py --cells cells --out rough200 --n 200 --seed 20261002 >> $LOG 2>&1 || { say "finalize_200 FAILED"; exit 4; }
say "step 9: rough_dft prepare (inputs + queue + budget; NOTHING SUBMITTED)"; $R/pyrun_rs.sh $R/rough_dft.py --dftdir dft prepare --manifest rough200/rough200_manifest.json >> $LOG 2>&1 || { say "rough_dft prepare FAILED"; exit 5; }
say "step 9b: first-batch dry run"; $R/pyrun_rs.sh $R/rough_dft.py --dftdir dft submit --first 8 --manifest rough200/rough200_manifest.json >> $LOG 2>&1
say "step 10: lineage figures"; $R/pyrun_rs.sh $R/render_rough_lineage.py --cells cells --manifest rough200/rough200_manifest.json --out figures >> $LOG 2>&1 || { say "render FAILED"; exit 6; }
say "report"; $R/pyrun_rs.sh $R/make_candidate_report.py --out rough200 --dftdir dft --cells cells >> $LOG 2>&1 || { say "report FAILED"; exit 7; }
say "CANDIDATE DELIVERABLES READY (frozen=false, nothing submitted)"
