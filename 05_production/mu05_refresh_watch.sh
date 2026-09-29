#!/bin/bash
# Waits until every NON-extension task is finished, then applies warm starts to the still-cold extension tasks
# (scripts/extend_mu05.py --refresh: pending tasks of the mu05_extension campaign only) and stops. Also stops after 12 h.
cd /anvil/scratch/x-rywang/Au_Cl; t0=$(date +%s)
while :; do
  left=$(scripts/pyrun.sh -c "import json;q=json.load(open('05_production/queue.json'));print(sum(t['status'] in ('pending','running') and t.get('campaign')!='mu05_extension' for t in q))" 2>/dev/null)
  el=$(( $(date +%s) - t0 ))
  if [ "${left:-1}" -eq 0 ] || [ $el -ge 43200 ]; then break; fi
  sleep 900
done
echo "=== $(date '+%m-%d %H:%M'): original tasks left = $left; running refresh"
scripts/pyrun.sh scripts/extend_mu05.py --refresh 2>&1 | tail -3
scripts/pyrun.sh -c "import json;q=json.load(open('05_production/queue.json'));e=[t for t in q if t.get('campaign')=='mu05_extension'];print('extension:',{s:sum(t['status']==s for t in e) for s in set(t['status'] for t in e)},'| cold starts:',sum(t['warm_start']=='cold start' for t in e))"
grep -cE "\[failed\]|\[children-error\]" 05_production/feeder.log | sed 's/^/failed+children-error lines in feeder.log: /'
