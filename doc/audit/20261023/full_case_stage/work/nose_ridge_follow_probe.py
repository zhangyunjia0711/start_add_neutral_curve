#!/usr/bin/env python3
"""Actual local-ridge path diagnostic at a stalled nose stationary root."""
import argparse,json
from pathlib import Path
from types import SimpleNamespace
from full_driver import ROOT,sha,snapshot
from ridge_path_survey import row_peaks,stamp

def main():
    p=argparse.ArgumentParser();p.add_argument('--source-file',type=Path,required=True)
    p.add_argument('--dw',type=float,required=True);p.add_argument('--steps',type=int,default=20)
    p.add_argument('--halfwidth',type=float,default=.004);p.add_argument('--db',type=float,default=.0000625)
    p.add_argument('--maximum-halfwidth',type=float,default=None,
        help='Explicit finite expansion cap for an unresolved edge maximum; every attempted window is retained')
    p.add_argument('--label',required=True);a=p.parse_args()
    target=ROOT/f'outputs/{a.label}.json'
    if target.exists():raise RuntimeError('Existing evidence is preserved')
    source=json.loads(a.source_file.read_text());failure=source['failures'][-1]
    point=failure['result']['trace'][-1]['center'];station=failure['station']
    solver=ROOT/'outputs/bin/START_LST_wall_local_v10'
    args=SimpleNamespace(station=station,solver=solver,solver_sha=sha(solver),timeout=90,
        halfwidth=a.halfwidth,db=a.db,bmax=80,label=a.label)
    report=dict(created=stamp(),source=str(a.source_file),source_sha256=sha(a.source_file),
        initial_actual=point,station=station,frequency_step=a.dw,steps=a.steps,
        local_beta_halfwidth=a.halfwidth,maximum_halfwidth=a.maximum_halfwidth,beta_step_max=a.db,rows=[],finished=False,
        qualification='Finite actual peak paths; not a neutral root until independently refined and verified')
    active=[dict(point=point,lineage='stalled_actual')]
    def save():target.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    snapshot();save()
    for j in range(a.steps+1):
        omega=point['omega_r']+j*a.dw;next_active=[]
        for origin in active:
            args.halfwidth=a.halfwidth
            while True:
                row=row_peaks(args,origin['point'],omega,origin['lineage']);report['rows'].append(row);save()
                if row['peaks'] or not row['edge_candidates'] or not all(v['complete'] for v in row['arms']):break
                if a.maximum_halfwidth is None or args.halfwidth>=a.maximum_halfwidth:break
                args.halfwidth=min(args.halfwidth*2,a.maximum_halfwidth)
            next_active.extend(dict(point=v['point'],lineage=v['lineage']) for v in row['peaks'])
        active=next_active;report.update(active=active,last_omega=omega,updated=stamp());save()
        print('NOSE RIDGE '+json.dumps(dict(omega=omega,peaks=[v['point'] for v in active])),flush=True)
        if not active or len(active)>4:break
    report['finished']=True;report['updated']=stamp();save();snapshot()

if __name__=='__main__':main()
