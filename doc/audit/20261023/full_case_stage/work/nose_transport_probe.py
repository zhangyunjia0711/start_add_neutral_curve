#!/usr/bin/env python3
"""Independent fixed-parameter adjacent-station jump-gate diagnostic."""
import argparse, json
from pathlib import Path
from types import SimpleNamespace
from full_driver import ROOT, sha, snapshot
from local_peak_x_continue import transport, eigen, stamp

def main():
    p=argparse.ArgumentParser();p.add_argument('--jump',type=float,default=8)
    p.add_argument('--label',default='nose_M1p_transport_j8_20261023')
    p.add_argument('--source-file',type=Path,default=ROOT/'outputs/local_peak_nose_M1p_v10_progress.json')
    p.add_argument('--origin',type=int,default=19);p.add_argument('--target',type=int,default=18)
    a=p.parse_args()
    out=ROOT/f'outputs/{a.label}.json'
    if out.exists():raise RuntimeError('Existing diagnostic is preserved')
    source=a.source_file
    point=next(v['point'] for v in json.loads(source.read_text())['records'] if v['station']==a.origin)
    solver=ROOT/'outputs/bin/START_LST_wall_local_v10'
    args=SimpleNamespace(solver=solver,solver_sha=sha(solver),label=a.label,eigen_jump=a.jump)
    report=dict(started=stamp(),source=str(source),source_sha256=sha(source),source_point=point,
                solver_sha256=args.solver_sha,jump=a.jump,minimum_overlap=.95,legs=[])
    def save():out.write_text(json.dumps(report,indent=2)+'\n')
    snapshot();save()
    try:
        target,name=transport(args,a.origin,a.target,point)
        report['legs'].append(dict(origin=a.origin,target=a.target,run=name,point=target,delta_alpha=abs(eigen(target)-eigen(point))))
        save()
        back,name=transport(args,a.target,a.origin,target)
        report['legs'].append(dict(origin=a.target,target=a.origin,run=name,point=back,delta_alpha=abs(eigen(back)-eigen(target)),
                                  reference_alpha_difference=abs(eigen(back)-eigen(point))))
        report['status']='both_legs_accepted'
    except Exception as e:report.update(status='unresolved',reason=repr(e))
    report['finished']=stamp();save();snapshot();print(json.dumps(report),flush=True)

if __name__=='__main__':main()
