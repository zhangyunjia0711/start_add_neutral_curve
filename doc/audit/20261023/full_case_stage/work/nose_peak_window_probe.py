#!/usr/bin/env python3
"""Inspect the actual nearby beta ridge before declaring a predicted peak lost."""
import argparse,json
from pathlib import Path
from types import SimpleNamespace
from full_driver import ROOT,sha,snapshot
from local_peak_x_continue import execute,seed,stamp

def main():
    p=argparse.ArgumentParser();p.add_argument('--source-file',type=Path,required=True)
    p.add_argument('--radius',type=float,default=.08);p.add_argument('--intervals',type=int,default=32)
    p.add_argument('--label',required=True);a=p.parse_args()
    out=ROOT/f'outputs/{a.label}.json'
    if out.exists():raise RuntimeError('Existing probe is preserved')
    d=json.loads(a.source_file.read_text());failure=d['failures'][-1]
    point=failure['predicted_actual'];station=failure['station'];w=point['omega_r'];b=point['beta_r']
    solver=ROOT/'outputs/bin/START_LST_wall_local_v10'
    args=SimpleNamespace(solver=solver,solver_sha=sha(solver),label=a.label,eigen_jump=8.)
    report=dict(created=stamp(),source=str(a.source_file),source_sha256=sha(a.source_file),station=station,
                center=point,radius=a.radius,anchored_window=failure['beta_window'],arms=[],peaks=[])
    def save():out.write_text(json.dumps(report,indent=2)+'\n')
    snapshot();save()
    rows={}
    for side in [-1,1]:
        bounds=[b-a.radius,b] if side<0 else [b,b+a.radius]
        arm=dict(side=side,bounds=bounds)
        try:
            values,name=execute(args,dict(station=station,seed=seed(point),
                omega_min=w,omega_max=w,omega_points=1,beta_min=bounds[0],beta_max=bounds[1],
                beta_points=a.intervals+1,reverse_inner=side<0),f'i{station}_arm{side}')
            arm.update(status='sampled',run=name,samples=values)
            for v in values:rows[round(v['beta_r'],11)]=v
        except Exception as e:arm.update(status='unknown',reason=repr(e))
        report['arms'].append(arm);save()
    ordered=sorted(rows.values(),key=lambda v:v['beta_r'])
    report['peaks']=[dict(point=v,left=ordered[j-1],right=ordered[j+1]) for j,v in enumerate(ordered)
                     if 0<j<len(ordered)-1 and v['growth']>ordered[j-1]['growth'] and v['growth']>ordered[j+1]['growth']]
    report['finished']=stamp();save();snapshot();print(json.dumps({k:v for k,v in report.items() if k!='arms'}),flush=True)

if __name__=='__main__':main()
