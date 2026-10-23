#!/usr/bin/env python3
"""Correct an actual near-neutral ridge sample with independently fine FD."""
import argparse,json
from pathlib import Path
from types import SimpleNamespace
from full_driver import ROOT,sha,snapshot
from local_peak_x_continue import correct,stamp

def main():
    p=argparse.ArgumentParser();p.add_argument('--ridge-file',type=Path,required=True)
    p.add_argument('--row',type=int,required=True);p.add_argument('--source',required=True)
    p.add_argument('--label',required=True);p.add_argument('--h',type=float,default=1e-6)
    a=p.parse_args();out=ROOT/f'outputs/{a.label}.json'
    if out.exists():raise RuntimeError('Existing result is preserved')
    archived=json.loads(a.ridge_file.read_text());peaks=archived['rows'][a.row]['peaks']
    if len(peaks)!=1:raise RuntimeError('Actual ridge seed is ambiguous')
    point=peaks[0]['point'];station=archived['station'];solver=ROOT/'outputs/bin/START_LST_wall_local_v10'
    args=SimpleNamespace(solver=solver,solver_sha=sha(solver),label=a.label,eigen_jump=16.,
        h=a.h,max_iterations=32,strict_gradient=False,beta_location_tol=1e-5,
        omega_min=-80,omega_max=80,omega_trust=.05,beta_trust=.05,
        peak_reseed_radius=0.,minimum_stencil_h=1e-7)
    snapshot();result=correct(args,station,point,point['omega_r'],point['beta_r'],
        [point['beta_r']-.05,point['beta_r']+.05],a.label)
    report=dict(created=stamp(),ridge_source=str(a.ridge_file),ridge_sha256=sha(a.ridge_file),
        actual_seed=point,results=[dict(job=dict(station=station,source=a.source,name=a.label),result=result)])
    out.write_text(json.dumps(report,indent=2)+'\n');snapshot()
    print(json.dumps(dict(station=station,status=result['status'],point=result.get('point'),reason=result.get('reason'))),flush=True)

if __name__=='__main__':main()
