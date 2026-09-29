#!/bin/bash
# Waits until every ORIGINAL (dataset_plan_v1 rev 2, +-0.2 V / reference) task is finished -- i.e. no non-extension task is
# pending or running -- then applies the pending warm starts and prints the final tally of that campaign. Also breaks on the
# first [failed]/[children-error] line, or after 24 h. The mu05_extension tasks keep running either way.
cd /anvil/scratch/x-rywang/Au_Cl; t0=$(date +%s)
while :; do
  left=$(scripts/pyrun.sh -c "import json;q=json.load(open('05_production/queue.json'));print(sum(t['status'] in ('pending','running') and t.get('campaign')!='mu05_extension' for t in q))" 2>/dev/null)
  nf=$(grep -cE "\[failed\]|\[children-error\]" 05_production/feeder.log)
  el=$(( $(date +%s) - t0 ))
  if [ "${left:-1}" -eq 0 ] || [ "$nf" -gt 0 ] || [ $el -ge 86400 ]; then break; fi
  sleep 300
done
echo "=== $(date '+%m-%d %H:%M'): original (+-0.2 V) tasks still open = $left ; failed/children-error lines = $nf"
echo "--- applying warm starts now that the neighbours are done:"
scripts/pyrun.sh scripts/extend_mu05.py --refresh 2>&1 | tail -3
echo "--- original campaign tally:"
scripts/pyrun.sh -c "
import json, re, collections
q=json.load(open('05_production/queue.json'))
o=[t for t in q if t.get('campaign')!='mu05_extension']; e=[t for t in q if t.get('campaign')=='mu05_extension']
print('original :',{s:sum(t['status']==s for t in o) for s in ('completed','reused','running','pending','failed')})
print('extension:',{s:sum(t['status']==s for t in e) for s in ('completed','running','pending','failed')})
cost=sum(t.get('elapsed_farm_min') or 0 for t in o)/60
print('original measured cost: %.0f node-h over %d computed tasks'%(cost,sum(1 for t in o if t.get('elapsed_farm_min'))))
fam=collections.Counter(t['family'] for t in o if t['status'] in ('completed','reused'))
print('per family:',dict(fam))
off=[]
for t in o:
    m=re.search(r'mu_e=\s*([-\d.]+)',str(t.get('result') or ''))
    if m: off.append(abs(float(m.group(1))-float(t['mu'])))
print('max |mu_e - TARGETMU| over original states: %.4f eV (n=%d)'%(max(off),len(off)))
"
echo "--- last farm lines:"; grep -E "\[(completed|launch|failed)\]" 05_production/feeder.log | tail -n 3
