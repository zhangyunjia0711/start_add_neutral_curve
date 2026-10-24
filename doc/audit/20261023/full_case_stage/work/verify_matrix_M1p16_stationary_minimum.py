#!/usr/bin/env python3
"""Native spatial verification of a matrix-found beta minimum, excluded from Gmax."""
import json,sys
from pathlib import Path
import numpy as np
from full_driver import ROOT,CASE,sha,snapshot
from local_peak_x_continue import stamp
from fd_envelope import stencil

def main():
    sys.path.insert(0,str(ROOT/'work/flat_neutral_normal_validation_20261023'))
    from validate_flat import solve
    source=ROOT/'outputs/matrix_neutral_M1p16_20261023.json';d=json.loads(source.read_text())
    c=d['native16_correction']['point'];out=ROOT/'outputs/native16_stationary_beta_minimum_20261023.json'
    if out.exists():raise RuntimeError('Existing verification preserved')
    expected=dict(omega_r=c['omega']['real'],beta_r=c['beta'],alpha_r=c['alpha'],alpha_i=0.)
    snapshot();fresh=solve('native16_stationary_minimum_fresh',CASE/'M1/in',400,16,expected,
        here=ROOT/'work/native16_stationary_beta_minimum_20261023',binary=Path(d['solver']),sha=d['solver_sha256'])
    point=dict(expected,alpha_r=fresh['alpha_r'],alpha_i=fresh['alpha_i'])
    diagnostics=[]
    for h,tight in [(1e-5,False),(5e-6,True)]:
        s=stencil(16,point,point['omega_r'],point['beta_r'],Path(d['solver']),'wall_local',h=h,tight=tight,
            label='native16_minimum_spatial_'+str(h))
        delta=np.linalg.solve(np.array(s['J']),-np.array(s['F']));s['location_correction']=delta.tolist();diagnostics.append(s)
    delta=abs(complex(fresh['alpha_r']-expected['alpha_r'],fresh['alpha_i']))
    passed=bool(fresh['passed'] and abs(fresh['growth'])<1e-7 and delta<1e-7 and
        all(s['J'][1][1]>0 and abs(s['location_correction'][0])<1e-4 and abs(s['location_correction'][1])<1e-5 for s in diagnostics))
    report=dict(created=stamp(),source_path=str(source),source_sha256=sha(source),solver=d['solver'],solver_sha256=d['solver_sha256'],
        point=dict(point,growth=fresh['growth'],iloc=16,x=float(fresh['metadata']['x_global_N'])),
        actual_native_fresh=fresh,alpha_difference_to_real_candidate=delta,spatial_h_h2=diagnostics,passed=passed,
        matrix_curvatures=[x['spatial_beta_curvature'] for x in d['native16_h_h2']],
        status='verified_stationary_beta_minimum' if passed else 'unresolved_minimum_verification',
        envelope_candidate=False,qualification='Actual stationary neutral point classified as a beta minimum; excluded from the maximum-growth beta envelope. Physical mode direction/mechanism and global spectrum remain unclassified.')
    out.write_text(json.dumps(report,indent=2)+'\n');snapshot()
    print(json.dumps(dict(passed=passed,point=report['point'],alpha_difference=delta,
        spatial_curvature=[s['J'][1][1] for s in diagnostics],spatial_correction=[s['location_correction'] for s in diagnostics])),flush=True)

if __name__=='__main__':main()
