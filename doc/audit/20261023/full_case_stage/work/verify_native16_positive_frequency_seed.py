#!/usr/bin/env python3
"""Independently solve the sign-conjugated candidate on original native16 H5."""
import json,sys
from pathlib import Path
import numpy as np
from full_driver import ROOT,CASE,sha,snapshot
from local_peak_x_continue import stamp

def main():
    sys.path.insert(0,str(ROOT/'work/flat_neutral_normal_validation_20261023'))
    from validate_flat import solve,vectors
    source=ROOT/'outputs/native16_free_beta_positive_seed_20261023.json';d=json.loads(source.read_text())
    p=d['point'];expected=dict(omega_r=-p['omega_r'],beta_r=-p['beta_r'],alpha_r=-p['alpha_r'],alpha_i=p['alpha_i'])
    out=ROOT/'outputs/native16_positive_frequency_seed_20261023.json'
    if out.exists():raise RuntimeError('Prior fresh verification preserved')
    snapshot();r=solve('native16_positive_frequency_fresh',CASE/'M1/in',400,16,expected,
        here=ROOT/'work/native16_positive_frequency_seed_20261023',binary=Path(d['solver']),sha=d['solver_sha256'])
    old=d['native_fixed_parameter_validation']['native']
    a=vectors(Path(old['native_vector_file']),400)[old['raw_index']]
    b=vectors(Path(r['native_vector_file']),400)[r['raw_index']]
    mac=float(abs(np.vdot(np.conj(a),b))**2/(np.vdot(a,a).real*np.vdot(b,b).real))
    delta=abs(complex(r['alpha_r']-expected['alpha_r'],r['alpha_i']-expected['alpha_i']))
    passed=bool(r['passed'] and delta<1e-8 and mac>=1-1e-10 and r['growth']>0)
    point=dict(omega_r=expected['omega_r'],beta_r=expected['beta_r'],alpha_r=r['alpha_r'],alpha_i=r['alpha_i'],growth=r['growth'],iloc=16)
    report=dict(created=stamp(),source_path=str(source),source_sha256=sha(source),family='free_beta_M1p_conjugate',
        station=16,coordinates='wall_local',point=point,solver=d['solver'],solver_sha256=d['solver_sha256'],
        expected_algebraic_candidate=expected,actual_native_validation=r,
        alpha_difference_to_candidate=delta,conjugated_raw_primitive_MAC=mac,passed=passed,
        physical_role='Independently fresh original-H5 positive-growth seed at positive frequency; not a neutral root or a global beta maximum')
    out.write_text(json.dumps(report,indent=2)+'\n');snapshot();print(json.dumps({k:report[k] for k in ['point','alpha_difference_to_candidate','conjugated_raw_primitive_MAC','passed']}),flush=True)

if __name__=='__main__':main()
