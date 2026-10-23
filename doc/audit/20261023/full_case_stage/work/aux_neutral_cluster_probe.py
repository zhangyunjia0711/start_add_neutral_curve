#!/usr/bin/env python3
"""Preserve finite nearby PEP candidates at an accepted auxiliary root."""
import argparse,json,re
from pathlib import Path
from full_driver import ROOT,run,sha,snapshot
from local_peak_x_continue import stamp,seed

def main():
    p=argparse.ArgumentParser();p.add_argument('--source-report',type=Path,required=True)
    p.add_argument('--column',type=int,required=True);p.add_argument('--label',required=True);a=p.parse_args()
    out=ROOT/f'outputs/{a.label}.json'
    if out.exists():raise RuntimeError('Existing evidence preserved')
    d=json.loads(a.source_report.read_text());record=next(x for x in d['records'] if x['auxiliary_column']==a.column)
    point=record['point'];b=point['beta_r'];w=point['omega_r']
    spec=dict(name=a.label,mode=d['template'],auxiliary_input_hashes=d['input_hashes'],station=a.column,
        solver=d['solver'],solver_sha256=d['solver_sha256'],coordinates='wall_local',task=33,seed=seed(point),seed_shift=1e-4,
        omega_min=w,omega_max=w,omega_points=1,beta_min=b,beta_max=b+1e-10,beta_points=2,
        anchor_row_seed=True,sample_only=True,adaptive=True,vector_tracking=True,minimum_overlap=.95,
        retry_eigenpairs=8,minimum_step=1e-6,step=.001,jump=8,eigen_tol=1e-7,timeout=180,
        options=['-pep_tol','1e-13','-pep_max_it','1000','-pep_nev','8','-pep_ncv','24',
                 '-st_pc_factor_mat_solver_type','superlu','-pep_view_vectors','binary:out/all_vectors.petsc'])
    snapshot();summary=run(spec)
    log=ROOT/'runs'/a.label/'solver.log';blocks=[];block=[]
    pat=r'^\s*(\d+)\s+\(\s*([^,]+),\s*([^\)]+)\)\s+The error is\s+(\S+)'
    for m in re.finditer(pat,log.read_text(),re.M):
        if int(m[1])==0 and block:blocks.append(block);block=[]
        block.append(dict(index=int(m[1]),alpha_r=float(m[2]),alpha_i=float(m[3]),residual=float(m[4])))
    if block:blocks.append(block)
    pairs=blocks[-1] if blocks else []
    alpha=complex(point['alpha_r'],point['alpha_i'])
    for v in pairs:v['distance_to_accepted_aux_root']=abs(complex(v['alpha_r'],v['alpha_i'])-alpha)
    report=dict(created=stamp(),source_report=str(a.source_report),source_sha256=sha(a.source_report),
        auxiliary_column=a.column,theta=record['theta'],actual_root=point,last_PEP_candidates=pairs,summary=summary,
        native_vector_file=str(ROOT/'runs'/a.label/'out/all_vectors.petsc'),
        qualification='Finite nearby eigenpairs at an auxiliary numerical profile, not complete spectrum or physical x candidates')
    out.write_text(json.dumps(report,indent=2)+'\n');snapshot();print(json.dumps(pairs),flush=True)

if __name__=='__main__':main()
