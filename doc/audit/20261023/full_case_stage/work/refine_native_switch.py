"""Refine both actual upper neutral roots at the native primary switch bracket."""
import concurrent.futures,json
from full_driver import ROOT,sha
from trace_envelope_x import transfer
from fd_envelope import refine
solver=ROOT/'outputs/bin/START_LST_wall_local_v10'

def solve(source):
    p=ROOT/f'outputs/neutral_x_cont_v7a_{source}_downstream_progress.json'
    data=json.loads(p.read_text())
    old=next(r['point'] for r in data['records'] if r['station']==371)
    actual,transport=transfer(371,372,old,solver,'primary_switch_'+source,True)
    result=refine(372,actual,actual['omega_r'],actual['beta_r'],solver,label='primary_switch_'+source,
        omega_bounds=(13,16),beta_bounds=(old['beta_r']-.5,old['beta_r']+.5),
        h=.0025,max_iterations=25)
    return dict(job=dict(station=372,source=source,name='primary_switch_'+source,
        origin_path=str(p),origin_sha256=sha(p),origin_point=old,transport=transport),result=result)

if __name__=='__main__':
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        results=list(pool.map(solve,['M1p','M2p']))
    output=dict(results=results,
        qualification='Two actual upper neutral beta maxima at original native station i372; compare omega.')
    (ROOT/'outputs/primary_crossover_native_i372.json').write_text(
        json.dumps(output,indent=2,allow_nan=False)+'\n')
    print(json.dumps([(r['job']['source'],r['result']['status'],r['result'].get('point')) for r in results]),flush=True)

