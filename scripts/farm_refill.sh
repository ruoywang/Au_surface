#!/bin/bash
# Keep short farm jobs in the backfill window (2026-09-27): our account priority is 75 (FAIRSHARE 0), so jobs start only
# via backfill, which on this cluster takes 4-10 h single-node jobs within ~4-6 h. Maintain TARGET 6-hour farms on
# standard/wholenode (alternating) and 2 farms on highmem (starts immediately, MaxJobsPU=2); resubmit as they finish.
# Stops when the production queue has nothing pending or running. Log: 05_production/farm/refill.log
cd /anvil/scratch/x-rywang/Au_Cl/05_production/farm || exit 1
TARGET=${TARGET:-8}; HM=2; i=0
log(){ echo "$(date '+%F %T') $*" | tee -a refill.log; }
while :; do
  n_short=$(squeue -u "$USER" -h -n AuCl_farm6h_standard -o %i | wc -l)
  n_hm=$(squeue -u "$USER" -h -n AuCl_farm_highmem -o %i | wc -l)
  pend=$(/anvil/scratch/x-rywang/Au_Cl/scripts/pyrun.sh -c "import json;q=json.load(open('../queue.json'));print(sum(t['status'] in ('pending','running') for t in q))" 2>/dev/null)
  if [ "${pend:-1}" -eq 0 ]; then log "queue drained (no pending/running tasks) -- refill loop exits"; break; fi
  while [ "$n_short" -lt "$TARGET" ]; do
    part=standard   # 2026-09-27: standard only -- our jobs are 5th of 12 pending there vs 371st of 1129 in wholenode (same nodes, separate queues)
    out=$(sbatch job-farm6h-$part 2>&1) && log "submitted 6h farm on $part: $out" || log "sbatch failed ($part): $out"
    n_short=$((n_short+1))
  done
  while [ "$n_hm" -lt "$HM" ]; do out=$(sbatch job-farm-highmem 2>&1) && log "submitted highmem farm: $out" || { log "highmem sbatch failed: $out"; break; }; n_hm=$((n_hm+1)); done
  sleep 600
done
