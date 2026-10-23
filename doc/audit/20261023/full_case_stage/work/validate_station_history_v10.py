"""Three-station actual-history checks, using the same fixed real omega/beta."""
import datetime,json
from validate_station_predictor_v10 import ROOT,V8,V9,V8HASH,V9HASH,CASES,path,sha,snapshot

snapshot();assert sha(V8)==V8HASH and sha(V9)==V9HASH
records={}
for base in CASES:
    case=dict(base,target=base['start']-2,outer_points=3)
    result=path(case,V9,V9HASH,True,'history')
    neutral=(ROOT/'runs'/result['spec']['name']/'out/neutral_run.log').read_text()
    rescue_lines=[line for line in neutral.splitlines() if 'station_shift_rescue stage=' in line]
    result['rescue_lines']=rescue_lines
    result['history_exercised']=bool(result['predictor_history_lines'])
    records[base['label']]=result
output=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),cases=records,
            passed=all(v['accepted'] and v['history_exercised'] for v in records.values()))
p=ROOT/'outputs/station_predictor_v10_history_validation.json';assert not p.exists();p.write_text(json.dumps(output,indent=2)+'\n')
snapshot();assert sha(V8)==V8HASH and sha(V9)==V9HASH
print('STATION HISTORY '+json.dumps(dict(passed=output['passed'],cases={k:dict(accepted=v['accepted'],point=v['point'],
      history_lines=v['predictor_history_lines'],rescue_lines=v['rescue_lines'],malloc_dump_bytes=v['malloc_dump_bytes']) for k,v in records.items()})),flush=True)
