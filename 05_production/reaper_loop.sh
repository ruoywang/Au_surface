#!/bin/bash
# Every 10 min: (a) re-queue tasks orphaned by a farm job that ended mid-run (scripts/reap_orphans.py), and
# (b) run extend_mu05.py --refresh so any newly pending task gets its warm start and its predicted NELECT start.
# Exits when the whole queue (both campaigns) has drained, or after 48 h.
cd /anvil/scratch/x-rywang/Au_Cl; t0=$(date +%s)
while :; do
  scripts/pyrun.sh scripts/reap_orphans.py >> 05_production/reaper.log 2>&1
  scripts/pyrun.sh scripts/extend_mu05.py --refresh >> 05_production/reaper.log 2>&1
  left=$(scripts/pyrun.sh -c "import json;q=json.load(open('05_production/queue.json'));print(sum(t['status'] in ('pending','running') for t in q))" 2>/dev/null)
  if [ "${left:-1}" -eq 0 ] || [ $(( $(date +%s) - t0 )) -ge 172800 ]; then break; fi
  sleep 600
done
echo "$(date '+%m-%d %H:%M') reaper loop exiting; pending+running=$left" >> 05_production/reaper.log
