#!/usr/bin/env python3
"""Track every observed local beta peak through small frequency increments.

Each new row starts from the previous row's actual peak beta/eigenpair. Separate
observed peaks retain separate lineages. This is a finite local-window path
audit, not exhaustive spectral enumeration or a certified global envelope.
"""
import argparse
import datetime
import hashlib
import json
import math
from pathlib import Path

import full_driver as driver

ROOT=driver.ROOT
OPTIONS=['-pep_tol','1e-12','-pep_max_it','1000','-pep_nev','1','-pep_ncv','8',
         '-st_pc_factor_mat_solver_type','superlu']


def stamp():
    return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()


def actual_seed(value):
    return {key:float(value[key]) for key in ['omega_r','beta_r','alpha_r','alpha_i']}


def point_run(args,origin,omega,blo,bhi,points,reverse=False):
    spec=dict(station=args.station,mode='M1',coordinates='wall_local',seed=actual_seed(origin),
        solver=str(args.solver),solver_sha256=args.solver_sha,seed_shift=1e-4,task=33,
        omega_min=float(omega),omega_max=float(omega),omega_points=1,
        beta_min=float(blo),beta_max=float(bhi),beta_points=points,
        reverse_inner=reverse,anchor_row_seed=True,sample_only=True,adaptive=True,
        vector_tracking=True,minimum_overlap=.95,retry_eigenpairs=8,minimum_step=1e-6,
        step=.08,jump=2,eigen_tol=1e-7,parameter_tol=1e-8,growth_tol=1e-9,
        options=OPTIONS,timeout=args.timeout)
    digest=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()[:12]
    spec['name']=args.label+'_'+digest
    summary=driver.run(spec)
    values=driver.table(ROOT/'runs'/spec['name']/'out/neutral_field.plt')
    fatal=bool(summary['exit'] not in [0,1] or summary['petsc_errors'] or summary['unused_options']
               or summary['coordinate_metadata_missing'])
    good=[dict(v,source_run=spec['name']) for v in values if v['status']==0 and
          all(math.isfinite(v[k]) for k in ['alpha_r','alpha_i','growth']) and not fatal]
    return dict(name=spec['name'],summary_path=str(ROOT/'runs'/spec['name']/'summary.json'),
        seconds=summary['seconds'],failed=summary['failed_field_records'],fatal=fatal,
        accepted_samples=good,complete=len(good)==points and summary['rejected']==0)


def row_peaks(args,origin,omega,lineage):
    center=origin['beta_r'];lo=max(-args.bmax,center-args.halfwidth)
    hi=min(args.bmax,center+args.halfwidth)
    arms=[]
    for blo,bhi,reverse in [(lo,center,True),(center,hi,False)]:
        arms.append(point_run(args,origin,omega,blo,bhi,
            max(2,math.ceil((bhi-blo)/args.db-1e-12)+1),reverse))
    unique={}
    for arm in arms:
        for value in arm['accepted_samples']:
            key=round(value['beta_r'],10)
            if key not in unique or value['growth']>unique[key]['growth']:unique[key]=value
    values=sorted(unique.values(),key=lambda v:v['beta_r'])
    peaks=[];edges=[]
    for j,value in enumerate(values):
        if j==0 or j==len(values)-1:
            if len(values)>1 and value['growth']>=values[1 if j==0 else -2]['growth']:
                edges.append(value)
            continue
        left,right=values[j-1],values[j+1]
        if value['growth']>=left['growth'] and value['growth']>=right['growth'] and (
            value['growth']>left['growth'] or value['growth']>right['growth']):
            chosen=value;refinement=None
            # Quadratic vertex is only a proposed parameter. The retained peak
            # is always an actual newly solved sample, never an interpolation.
            x0,x1,x2=left['beta_r'],value['beta_r'],right['beta_r']
            y0,y1,y2=left['growth'],value['growth'],right['growth']
            den=(x0-x1)*(x0-x2)*(x1-x2)
            if abs(den)>1e-15:
                qa=(x2*(y1-y0)+x1*(y0-y2)+x0*(y2-y1))/den
                qb=(x2*x2*(y0-y1)+x1*x1*(y2-y0)+x0*x0*(y1-y2))/den
                proposed=-qb/(2*qa) if qa<0 else x1
                if x0<proposed<x2 and abs(proposed-x1)>1e-9:
                    refinement=point_run(args,value,omega,proposed-1e-9,proposed+1e-9,3)
                    if refinement['accepted_samples']:
                        fresh=max(refinement['accepted_samples'],key=lambda v:v['growth'])
                        if fresh['growth']>=chosen['growth']:chosen=fresh
            peaks.append(dict(point=chosen,coarse_point=value,refinement=refinement,
                lineage=lineage+f'/w{omega:.6g}/peak{len(peaks)}',
                local_window_complete=all(arm['complete'] for arm in arms)))
    return dict(omega=omega,origin=origin,origin_lineage=lineage,beta_window=[lo,hi],
        arms=arms,peaks=peaks,edge_candidates=edges,
        status='observed_local_peaks' if peaks else 'unknown_no_resolved_interior_peak',
        note='Edge maxima and incomplete windows remain unresolved; no other peak is silently substituted.')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--solver',type=Path,required=True);p.add_argument('--solver-sha',required=True)
    p.add_argument('--station',type=int,default=312);p.add_argument('--wstart',type=float,default=4)
    p.add_argument('--wend',type=float,default=12);p.add_argument('--dw',type=float,default=.5)
    p.add_argument('--halfwidth',type=float,default=2);p.add_argument('--db',type=float,default=.25)
    p.add_argument('--bmax',type=float,default=30);p.add_argument('--timeout',type=float,default=90)
    p.add_argument('--seed-file',type=Path,
        help='Actual same-station seed record or wrapper with point; source file is hashed and retained')
    p.add_argument('--resume-from',type=Path,
        help='Start from every actual retained peak at --wstart in an archived ridge-path result')
    p.add_argument('--label',default='ridge_path_i312_v7');args=p.parse_args()
    args.solver=args.solver.resolve();assert driver.sha(args.solver)==args.solver_sha
    driver.snapshot()
    proof_path=ROOT/'outputs/path_competition_v7/summary.json'
    proof=json.loads(proof_path.read_text())
    if args.resume_from:
        archived=json.loads(args.resume_from.read_text())
        assert archived['station']==args.station and archived['coordinates']=='wall_local'
        initial=[dict(point=peak['point'],lineage=peak['lineage']+'/resumed')
            for row in archived['rows'] if abs(row['omega']-args.wstart)<1e-10 for peak in row['peaks']]
        assert initial,'No actual local peak at the requested restart frequency'
        provenance=dict(path=str(args.resume_from),sha256=driver.sha(args.resume_from),
            kind='actual retained local peaks at restart frequency')
    elif args.seed_file:
        archived=json.loads(args.seed_file.read_text());record=archived.get('point',archived)
        assert int(round(record.get('iloc',archived.get('station',-1))))==args.station
        assert abs(record['omega_r']-args.wstart)<1e-10
        assert record.get('status',0)==0
        initial=[dict(point=record,lineage=archived.get('lineage','actual_seed'))]
        provenance=dict(path=str(args.seed_file),sha256=driver.sha(args.seed_file),
            kind='actual independently documented same-station seed',record=archived)
    else:
        seed=actual_seed(proof['origin'])
        assert args.station==312 and abs(args.wstart-seed['omega_r'])<1e-12
        initial=[dict(point=seed,lineage='actual_M1_conjugate_peak')]
        provenance=dict(path=str(proof_path),sha256=driver.sha(proof_path),seed=seed)
    output=dict(created=stamp(),station=args.station,coordinates='wall_local',
        solver=str(args.solver),solver_sha256=args.solver_sha,
        seed_provenance=provenance,
        domain=dict(omega=[args.wstart,args.wend],beta=[-args.bmax,args.bmax]),
        local_window_halfwidth=args.halfwidth,frequency_step=args.dw,beta_step_max=args.db,
        qualification='Observed local beta-peak lineages; all separated peaks retained; incomplete windows and edge candidates unknown, not globally stable',
        rows=[],active=initial,comparisons=[])
    target=ROOT/f'outputs/{args.label}.json'
    n=round((args.wend-args.wstart)/args.dw)
    assert abs(args.wstart+n*args.dw-args.wend)<1e-10
    for step in range(n+1):
        omega=args.wstart+step*args.dw;active=[]
        for parent in output['active']:
            row=row_peaks(args,parent['point'],omega,parent['lineage'])
            output['rows'].append(row)
            for peak in row['peaks']:
                active.append(dict(point=peak['point'],lineage=peak['lineage']))
                if args.station==312 and abs(omega-8)<1e-10:
                    reference=proof['paths']['diagonal_beta14']['point'];beta=reference['beta_r']
                    comparison=point_run(args,peak['point'],omega,beta-1e-9,beta+1e-9,3)
                    if comparison['accepted_samples']:
                        point=min(comparison['accepted_samples'],key=lambda v:abs(v['beta_r']-beta))
                        comparison['same_beta_alpha_difference']=abs(complex(
                            point['alpha_r']-reference['alpha_r'],point['alpha_i']-reference['alpha_i']))
                        comparison['point']=point
                    output['comparisons'].append(comparison)
        output['active']=active;output['updated']=stamp();output['last_omega']=omega
        temp=target.with_suffix('.json.tmp');temp.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n');temp.replace(target)
        print('RIDGE PATH '+json.dumps(dict(omega=omega,active_peaks=len(active),
            peaks=[{k:v['point'][k] for k in ['beta_r','alpha_r','alpha_i','growth']} for v in active])),flush=True)
        if not active:break
    driver.snapshot()


if __name__=='__main__':main()
