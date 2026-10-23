"""Test stronger eigensolver tolerances on previously rejected evaluations."""
import concurrent.futures, json
from full_driver import ROOT, PREVIOUS, run, snapshot
checks=snapshot();(ROOT/'outputs/input_integrity_start.json').write_text(json.dumps(checks,indent=2)+'\n')
seeds={v['station']:v for v in json.loads((PREVIOUS/'outputs/seeds.json').read_text())}
specs=[dict(name='tight_i%d_w%.10f'%(i,w),station=i,omega_min=w,omega_max=w,omega_points=1,beta_min=-12.,beta_max=12.,beta_points=49,step=.08,seed=seeds[i]) for i,w in [(351,8.),(381,9.762603680823121),(381,15.05840923439856)]]
results=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
    for job in concurrent.futures.as_completed([pool.submit(run,s) for s in specs]):
        results.append(job.result());(ROOT/'outputs/tight_probe_summary.json').write_text(json.dumps(results,indent=2)+'\n')
assert snapshot()==checks
