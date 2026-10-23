"""Compare the two existing neutral beta maxima at the native switch bracket."""
import concurrent.futures,json
from full_driver import ROOT,sha
from trace_envelope_x import transfer
from fd_envelope import refine
solver=ROOT/'outputs/bin/START_LST_wall_local_v10'
paths={'tail_beta4':ROOT/'outputs/neutral_tail_beta4_v10_tail_beta4_upstream_progress.json',
       'tail_beta6':ROOT/'outputs/neutral_x_cont_v7a_M2p_downstream_progress.json'}

def solve(source):
    p=paths[source];data=json.loads(p.read_text());old=next(r['point'] for r in data['records'] if r['station']==772)
    actual,transport=transfer(772,773,old,solver,'tail_switch_'+source,True)
    result=refine(773,actual,actual['omega_r'],actual['beta_r'],solver,label='tail_switch_'+source,
        omega_bounds=(4.5,5.5),beta_bounds=(old['beta_r']-.5,old['beta_r']+.5),h=.0025,max_iterations=25)
    return dict(job=dict(station=773,source=source,name='tail_switch_'+source,
        origin_path=str(p),origin_sha256=sha(p),origin_point=old,transport=transport),result=result)

with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(solve,paths))
output=dict(results=results,qualification='Two actual beta stationary upper roots; compare their frequencies, not fixed-frequency growth.')
(ROOT/'outputs/tail_crossover_native_i773.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
print(json.dumps([(r['job']['source'],r['result']['status'],r['result'].get('point')) for r in results]),flush=True)
