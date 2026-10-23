"""Recompute two stopped local candidates, retaining full Newton traces."""
import concurrent.futures,json
from pathlib import Path
from full_driver import ROOT
from trace_envelope_x import transfer
from fd_envelope import refine

def solve(label):
    progress=ROOT/f'outputs/neutral_x_cont_v7a_{label}_progress.json'
    d=json.loads(progress.read_text());f=d['failures'][-1];solver=Path(d['spec']['solver'])
    actual,transport=transfer(f['from_station'],f['station'],f['seed'],solver,'root_gap_diagnostic_'+label)
    result=refine(f['station'],actual,actual['omega_r'],actual['beta_r'],solver,
        label='root_gap_diagnostic_'+label,omega_bounds=(0,400),beta_bounds=(-200,200),max_iterations=25)
    out=dict(source_progress=str(progress),source_failure=f,transport=transport,result=result)
    (ROOT/f'outputs/root_gap_diagnostic_{label}.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(label=label,status=result['status'],reason=result.get('reason'),
        trace=[dict(omega=r['center']['omega_r'],beta=r['center']['beta_r'],F=r['F'],J=r['J'],h=r['h']) for r in result['trace']])),flush=True)

with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    list(pool.map(solve,['M2p_upstream','M1c_downstream']))
