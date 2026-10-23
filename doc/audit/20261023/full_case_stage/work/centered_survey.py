#!/usr/bin/env python3
"""Finite real-beta survey launched from actual same-station local-wall seeds.

Frequency proceeds outward from each seed; beta proceeds in two arms outward
from the seed beta. The first (central-beta) eigenpair anchors each next omega
row. Four search families and their joint conjugates are provenance labels;
missing families/cells remain explicit unknowns.
"""
import argparse
import concurrent.futures
import configparser
import csv
import datetime
import json
import math
from pathlib import Path
import re
import time

import full_driver as driver

ROOT=driver.ROOT
FAMILIES=['M1','M2','M3_local_a','M3_local_b']
OPTIONS=['-pep_tol','1e-12','-pep_max_it','1000','-pep_nev','1','-pep_ncv','8',
         '-st_pc_factor_mat_solver_type','superlu']


def stamp():
    return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()


def coarse_stations(grid):
    """139 original stations prescribed by physical-distance scope audit."""
    targets=list(range(-35,0,2))+[0]+list(range(0,161,2))+list(range(200,1081,40))+[1100]
    stations={min(range(len(grid)),key=lambda j:abs(35*grid[j]-x))+1 for x in targets}
    stations.update([1,2,21,61,119,240,241,312,313,317,321,331,341,351,361,371,
                     381,418,419,464,481,599,711,841])
    return sorted(stations)


def omega_segments(seed_omega,wmax,dw):
    """Common grid 0:dw:wmax, in outward traversal order from seed omega."""
    n=round(wmax/dw)
    assert abs(n*dw-wmax)<1e-10
    values=[j*dw for j in range(n+1)]
    lower=[v for v in values if v<=seed_omega+1e-12]
    upper=[v for v in values if v>=seed_omega-1e-12]
    result=[]
    if lower and len(lower)>1: result.append(('wdn',lower[0],lower[-1],len(lower),True))
    if upper: result.append(('wup',upper[0],upper[-1],len(upper),False))
    # A seed above the window has no upper segment. A singleton lower endpoint
    # still needs to be evaluated, e.g. a zero-size requested frequency window.
    if lower and not upper and len(lower)==1: result.append(('wdn',lower[0],lower[0],1,True))
    return result


def build_plan(args,axial,catalogue):
    grid=catalogue['grids']['main_841']['x_nondim']
    stations=args.stations or coarse_stations(grid)
    lookup={(v['family'],v['iloc']):v for v in axial.get('records',[])
            if v.get('status')==0 and v.get('coordinates')=='wall_local'}
    jobs=[];unknown=[]
    binary_hash=driver.sha(args.solver)
    for station in stations:
        assert 1<=station<=841
        for family in args.families:
            actual=lookup.get((family,station))
            if actual is None:
                unknown.append(dict(station=station,x=grid[station-1],x_mm=35*grid[station-1],
                    family=family,status='unknown',reason='no accepted same-station local-wall axial seed'))
                continue
            assert abs(actual['x']-grid[station-1])<1e-10
            for conjugate in [False,True]:
                seed={key:actual[key] for key in ['omega_r','beta_r','alpha_r','alpha_i']}
                if conjugate:
                    for key in ['omega_r','beta_r','alpha_r']:seed[key]=-seed[key]
                beta0=seed['beta_r']
                assert abs(beta0)<args.bmax,'Search domain must contain every seed beta'
                source=dict(kind='same_station_actual_local_axial_seed',record=actual,
                    joint_conjugate_derived=conjugate,
                    conjugacy='(omega,beta,alpha)->(-omega,-beta,-conj(alpha))',
                    physical_mode_identity='unclassified continuation of named search start')
                for wlabel,wlo,whi,nw,reverse in omega_segments(seed['omega_r'],args.wmax,args.dw):
                    for blabel,blo,bhi,breverse in [('bleft',-args.bmax,beta0,True),
                                                    ('bright',beta0,args.bmax,False)]:
                        nb=math.ceil((bhi-blo)/args.db-1e-12)+1
                        label=f'{family}{"c" if conjugate else "p"}'
                        jobs.append(dict(name=f'{args.phase}_{label}_i{station:03d}_{wlabel}_{blabel}_{binary_hash[:8]}',
                            family=family,conjugate=conjugate,station=station,mode='M1',
                            seed=seed,seed_shift=1e-4,seed_provenance=source,
                            solver=str(args.solver),solver_sha256=binary_hash,coordinates='wall_local',
                            task=33,sample_only=True,omega_min=wlo,omega_max=whi,omega_points=nw,
                            beta_min=blo,beta_max=bhi,beta_points=nb,beta_center=beta0,
                            frequency_direction=wlabel,beta_arm=blabel,reverse_outer=reverse,
                            reverse_inner=breverse,anchor_row_seed=True,adaptive=True,vector_tracking=True,
                            minimum_overlap=.95,minimum_step=1e-5,retry_eigenpairs=8,step=args.step,
                            jump=2,options=OPTIONS,timeout=args.timeout))
    return dict(created=stamp(),phase=args.phase,solver=str(args.solver),solver_sha256=binary_hash,
        coordinates='wall_local',grid='M1 main_841',stations=stations,families=args.families,
        domain=dict(omega=[0,args.wmax],beta=[-args.bmax,args.bmax],dw=args.dw,db_max=args.db),
        seed_source=str(args.axial),seed_source_sha256=driver.sha(args.axial),
        seed_source_solver_sha256=axial['solver_sha256'],
        extra_seed_sources=[dict(path=str(p),sha256=driver.sha(p)) for p in args.extra_axial],
        source_sha256={str(p):driver.sha(p) for p in [Path(__file__),Path(driver.__file__)]},
        missing_same_station_seeds=unknown,jobs=jobs,
        qualification='Finite declared-family real-beta search; unresolved starts/cells remain unknown, not stable')


def summarize(summary):
    spec=summary['spec'];directory=ROOT/'runs'/spec['name']
    cfg=configparser.ConfigParser();cfg.read(directory/'config.ini');neutral=cfg['neutral']
    istrue=lambda v:v.strip().lower() in ('t','true','yes','1')
    configuration_ok=all(istrue(neutral.get(k,''))==bool(spec[k]) for k in
        ['sample_only','vector_tracking','anchor_row_seed','reverse_outer','reverse_inner'])
    configuration_ok=configuration_ok and istrue(neutral.get('adaptive_steps',''))
    configuration_ok=configuration_ok and int(neutral.get('retry_eigenpairs','0'))==8
    configuration_ok=configuration_ok and abs(float(neutral.get('minimum_overlap','0'))-.95)<1e-12
    field=driver.table(directory/'out/neutral_field.plt');rows=[]
    run_bad=bool(summary.get('exit') not in [0,1] or not configuration_ok or summary.get('petsc_errors') or summary.get('unused_options')
        or summary.get('coordinate_metadata_missing') or summary['exit']==124)
    for row in range(1,spec['omega_points']+1):
        values=[v for v in field if v['row']==row]
        frac=(row-1)/max(1,spec['omega_points']-1)
        if spec['reverse_outer']:frac=1-frac
        omega=spec['omega_min']+frac*(spec['omega_max']-spec['omega_min'])
        valid=[v for v in values if v['status']==0 and all(math.isfinite(v[k]) for k in
            ['omega_r','beta_r','alpha_r','alpha_i','growth'])]
        reasons=[]
        if run_bad:reasons.append('run-level PETSc/options/metadata/runtime failure')
        if len(values)!=spec['beta_points']:reasons.append('unexpected number of sampled beta points')
        if any(abs(v['omega_r']-omega)>1e-9 for v in values):reasons.append('omega row ordering mismatch')
        if values and abs(values[0]['beta_r']-spec['beta_center'])>1e-9:
            reasons.append('row did not start at central beta')
        central=next((v for v in values if abs(v['beta_r']-spec['beta_center'])<1e-9),None)
        if central is None or central['status']!=0:reasons.append('central-beta row anchor unresolved')
        failed=sum(v['status']!=0 for v in values)
        if failed:reasons.append('one or more beta samples unresolved')
        winner=max(valid,key=lambda v:v['growth']) if valid else None
        rows.append(dict(row=row,omega=omega,status='sampled' if valid and not reasons else 'unknown',
            reasons=reasons,samples=len(values),accepted_samples=len(valid),failed=failed,
            winner=winner,central=central,
            beta_edge=bool(winner and abs(abs(winner['beta_r'])-max(abs(spec['beta_min']),abs(spec['beta_max'])))<1e-8)))
    return dict(name=spec['name'],family=spec['family'],conjugate=spec['conjugate'],
        station=spec['station'],frequency_direction=spec['frequency_direction'],beta_arm=spec['beta_arm'],
        beta_center=spec['beta_center'],coordinates='wall_local',summary_path=str(directory/'summary.json'),
        seconds=summary['seconds'],terminal=summary['terminal'],exit=summary['exit'],
        actual_configuration_verified=configuration_ok,
        rejected=summary.get('rejected'),petsc_errors=summary.get('petsc_errors'),rows=rows)


def save_progress(args,plan,results,started):
    output=dict(created=plan['created'],updated=stamp(),phase=args.phase,
        coordinates='wall_local',elapsed_seconds=round(time.time()-started,3),completed=len(results),
        planned=len(plan['jobs']),stations=plan['stations'],domain=plan['domain'],families=plan['families'],
        missing_same_station_seeds=plan['missing_same_station_seeds'],results=results,
        plan_path=str(ROOT/f'outputs/{args.phase}_plan.json'))
    tmp=ROOT/f'outputs/{args.phase}_progress.tmp'
    tmp.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    tmp.replace(ROOT/f'outputs/{args.phase}_progress.json')
    print('CENTERED PROGRESS '+json.dumps({k:output[k] for k in ['phase','completed','planned','elapsed_seconds']}),flush=True)


def export_boundary_report(phase):
    """Observed finite-grid boundary values; unresolved points remain explicit."""
    plan=json.loads((ROOT/f'outputs/{phase}_plan.json').read_text())
    progress=json.loads((ROOT/f'outputs/{phase}_progress.json').read_text())
    finished={v['name']:v for v in progress['results']}
    groups={}
    for spec in plan['jobs']:
        key=spec['station'],spec['family'],spec['conjugate']
        groups.setdefault(key,[]).append(spec)
    rows=[]
    bmax=plan['domain']['beta'][1];wmax=plan['domain']['omega'][1]
    for (station,family,conjugate),specs in sorted(groups.items()):
        samples=[];edge_complete={side:True for side in ['omega_max','beta_min','beta_max']}
        for spec in specs:
            field=driver.table(ROOT/'runs'/spec['name']/'out/neutral_field.plt')
            done=finished.get(spec['name'])
            valid_run=done and done.get('actual_configuration_verified',False) and not done.get('petsc_errors')
            for value in field:
                samples.append(dict(value,source_run=spec['name']))
            for side in edge_complete:
                if side=='omega_max':
                    participates=spec['omega_min']-1e-9<=wmax<=spec['omega_max']+1e-9
                    values=[v for v in field if abs(v['omega_r']-wmax)<1e-9]
                    expected=spec['beta_points']
                else:
                    edge=-bmax if side=='beta_min' else bmax
                    participates=spec['beta_min']-1e-9<=edge<=spec['beta_max']+1e-9
                    values=[v for v in field if abs(v['beta_r']-edge)<1e-9]
                    expected=spec['omega_points']
                if participates and (not valid_run or len(values)!=expected or any(v['status']!=0 for v in values)):
                    edge_complete[side]=False
        valid=[v for v in samples if v['status']==0 and math.isfinite(v['growth'])]
        for side in ['all','omega_max','beta_min','beta_max']:
            values=valid if side=='all' else [v for v in valid if
                (abs(v['omega_r']-wmax)<1e-9 if side=='omega_max' else
                 abs(v['beta_r']-(-bmax if side=='beta_min' else bmax))<1e-9)]
            winner=max(values,key=lambda v:v['growth']) if values else None
            complete=(all(spec['name'] in finished for spec in specs) and all(v['status']==0 for v in samples)) if side=='all' else edge_complete[side]
            rows.append(dict(station=station,family=family,conjugate=conjugate,boundary=side,
                status='sampled' if complete and values else 'unknown',samples=len(values),
                max_growth=winner['growth'] if winner else None,omega=winner['omega_r'] if winner else None,
                beta=winner['beta_r'] if winner else None,alpha_r=winner['alpha_r'] if winner else None,
                alpha_i=winner['alpha_i'] if winner else None,source_run=winner['source_run'] if winner else None,
                observed_positive=bool(winner and winner['growth']>1e-7)))
    result=dict(updated=stamp(),phase=phase,completed=progress['completed'],planned=progress['planned'],
        full_short_finished=progress['completed']==progress['planned'],domain=plan['domain'],rows=rows,
        arm_seconds=[dict(name=r['name'],seconds=r.get('seconds'),rejected=r.get('rejected'))
                     for r in sorted(progress['results'],key=lambda r:r.get('seconds',0),reverse=True)],
        qualification='Observed sampled lower bounds of the named finite search; edge positivity signals a need to test domain extension, unknown is not stable')
    path=ROOT/f'outputs/{phase}_boundary_report.json'
    path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    with path.with_suffix('.csv').open('w',newline='') as out:
        writer=csv.DictWriter(out,fieldnames=list(rows[0]) if rows else ['station'])
        writer.writeheader();writer.writerows(rows)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--solver',type=Path,required=True);p.add_argument('--solver-sha',required=True)
    p.add_argument('--axial',type=Path,default=ROOT/'outputs/axial_wall_local.json')
    p.add_argument('--extra-axial',type=Path,nargs='*',default=[])
    p.add_argument('--seed-solver-sha',help='Explicit historical actual-seed solver hash; seed provenance is retained')
    p.add_argument('--families',nargs='+',default=FAMILIES,
        choices=FAMILIES+['M3near'])
    p.add_argument('--phase',default='centered_coarse_v6');p.add_argument('--stations',nargs='*',type=int)
    p.add_argument('--plan-from',type=Path,
        help='Reprioritize an archived plan; reuse completed runs and retain explicitly interrupted arms as unknown')
    p.add_argument('--bmax',type=float,default=20);p.add_argument('--wmax',type=float,default=40)
    p.add_argument('--dw',type=float,default=4);p.add_argument('--db',type=float,default=1)
    p.add_argument('--step',type=float,default=.4);p.add_argument('--workers',type=int,default=6)
    p.add_argument('--timeout',type=float,default=7200);p.add_argument('--execute',action='store_true')
    args=p.parse_args();args.solver=args.solver.resolve()
    assert args.solver.is_file() and driver.sha(args.solver)==args.solver_sha
    version=re.search(r'_v(\d+)',args.solver.name)
    assert version and int(version[1])>=6,'Centered traversal requires reviewed v6-or-later API'
    driver.snapshot()
    axial=json.loads(args.axial.read_text())
    assert axial['solver_sha256']==(args.seed_solver_sha or args.solver_sha), 'Historical seed solver must be explicitly confirmed'
    for path in args.extra_axial:
        extra=json.loads(path.read_text())
        assert extra['solver_sha256'] in [args.solver_sha,args.seed_solver_sha], 'Unconfirmed extra seed executor hash'
        axial['records'].extend(extra.get('records',[]))
    catalogue=json.loads((ROOT/'outputs/seed_catalog.json').read_text())
    if args.plan_from:
        plan=json.loads(args.plan_from.read_text())
        assert plan['solver_sha256']==args.solver_sha
        plan['replanned_from']=dict(path=str(args.plan_from),sha256=driver.sha(args.plan_from))
        plan['phase']=args.phase;plan['created']=stamp()
        priority={119:0,241:1,841:2,599:3,312:4,317:5,381:6}
        plan['jobs'].sort(key=lambda s:(0 if s['family']=='M3near' else 1,
            priority.get(s['station'],100),s['family'],s['name']))
        for spec in plan['jobs']:
            if not (ROOT/'runs'/spec['name']).exists():spec['timeout']=args.timeout
        plan['replanning_notes']='M3near and flat representative boundaries first; interrupted original arms stay unknown; unstarted arm timeouts reduced without rewriting historical runs.'
    else:plan=build_plan(args,axial,catalogue)
    plan_path=ROOT/f'outputs/{args.phase}_plan.json'
    if plan_path.exists():
        old=json.loads(plan_path.read_text())
        assert old['jobs']==plan['jobs'],'Plan changed under existing run label; choose a new phase'
    else:plan_path.write_text(json.dumps(plan,indent=2,allow_nan=False)+'\n')
    print('CENTERED PLAN '+json.dumps(dict(jobs=len(plan['jobs']),stations=len(plan['stations']),
        missing_family_stations=len(plan['missing_same_station_seeds']),path=str(plan_path))),flush=True)
    if not args.execute:return
    started=time.time();results=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        pending={pool.submit(driver.run,spec):spec for spec in plan['jobs']}
        for future in concurrent.futures.as_completed(pending):
            spec=pending[future]
            try:result=summarize(future.result())
            except Exception as error:
                result=dict(name=spec['name'],station=spec['station'],family=spec['family'],
                    conjugate=spec['conjugate'],status='exception',exception=repr(error))
            results.append(result);save_progress(args,plan,results,started)
    driver.snapshot()


if __name__=='__main__':main()
