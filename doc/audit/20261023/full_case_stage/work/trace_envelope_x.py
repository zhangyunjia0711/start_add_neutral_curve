#!/usr/bin/env python3
"""Follow known neutral stationary candidates along original wall x.

Every new station first receives a true fixed-(omega,beta) iloc continuation
with vector matching, then a local two-equation refinement. This is a set of
individual candidate branches; a separate multi-path/domain competition audit
is required to turn it into an envelope. Failures never imply stability.
"""
import argparse
import concurrent.futures
import datetime
import hashlib
import json
from pathlib import Path
import time

from full_driver import ROOT,run,table,sha,snapshot
from centered_survey import coarse_stations,OPTIONS
from fd_envelope import refine


def stamp():return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()


def transfer(origin,target,seed,solver,label,station_predictor=False):
    if target==origin:return dict(seed),None
    spec=dict(station=target,seed_station=origin,mode='M1',
        seed={k:float(seed[k]) for k in ('omega_r','beta_r','alpha_r','alpha_i')},seed_shift=1e-4,
        solver=str(solver),solver_sha256=sha(solver),coordinates='wall_local',
        task=33,outer_parameter='iloc',inner_parameter='beta',
        outer_min=min(origin,target),outer_max=max(origin,target),outer_points=2,
        inner_min=seed['beta_r'],inner_max=seed['beta_r']+1e-8,inner_points=2,
        reverse_outer=target<origin,anchor_row_seed=True,sample_only=True,
        step=.2,jump=4,adaptive=True,vector_tracking=True,minimum_overlap=.95,
        minimum_step=1e-6,retry_eigenpairs=8,options=OPTIONS,timeout=600)
    if station_predictor:spec['station_predictor']=True
    token=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()[:12]
    spec['name']=f'{label}_transport_i{origin:03d}_{target:03d}_{token}'
    summary=run(spec)
    if summary['petsc_errors'] or summary['unused_options'] or summary['coordinate_metadata_missing'] or summary['exit']==124:
        raise RuntimeError('Transport runtime/configuration error: '+spec['name'])
    rows=table(ROOT/'runs'/spec['name']/'out/neutral_field.plt')
    candidates=[v for v in rows if abs(v['iloc']-target)<1e-7 and abs(v['beta_r']-seed['beta_r'])<1e-10 and v['status']==0]
    if not candidates:raise RuntimeError('No matched target station in '+spec['name'])
    return candidates[0],spec['name']


def chain(spec):
    path=ROOT/f'outputs/{spec["name"]}_progress.json'
    if path.exists():
        data=json.loads(path.read_text());assert data['spec']==spec
        if data.get('finished'):return data
    else:data=dict(created=stamp(),spec=spec,records=[],failures=[],finished=False,
        qualification='Individual continued neutral candidates; no global envelope completeness claim.')
    accepted={r['station']:r for r in data['records'] if r.get('status')=='candidate_converged'}
    current=dict(spec['seed']);station=spec['seed_station']
    targets=spec['stations'];direction=1 if targets[-1]>station else -1

    def save():
        data['updated']=stamp();tmp=path.with_suffix('.tmp')
        tmp.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');tmp.replace(path)

    def attempt(target,depth=0):
        nonlocal current,station
        if target in accepted:
            record=accepted[target];current=record['point'];station=target;return True
        start=station;origin=dict(current);transport=None;result=None
        try:
            actual,transport=transfer(start,target,origin,Path(spec['solver']),spec['name'],spec.get('station_predictor',False))
            result=refine(target,actual,actual['omega_r'],actual['beta_r'],Path(spec['solver']),
                'wall_local',label=spec['name']+f'_i{target:03d}',
                omega_bounds=(spec.get('omega_min',0),spec['omega_max']),beta_bounds=(-spec['beta_max'],spec['beta_max']),
                h=spec.get('fd_h',.02),max_iterations=spec.get('max_iterations',15))
            if result['status']!='candidate_converged':raise RuntimeError(result.get('reason','root unresolved'))
            record=dict(station=target,status=result['status'],point=result['point'],
                transport_from=start,transport_run=transport,result=result,requested_coarse=target in targets)
            data['records'].append(record);accepted[target]=record
            current=result['point'];station=target;save()
            print(json.dumps(dict(chain=spec['name'],station=target,omega=current['omega_r'],beta=current['beta_r'],
                growth=current['growth'],coarse=target in targets)),flush=True)
            return True
        except Exception as error:
            failure=dict(station=target,from_station=start,seed=origin,transport_run=transport,
                depth=depth,reason=repr(error),status='unknown',refinement=result)
            data['failures'].append(failure);save()
            if abs(target-start)<=1:return False
            middle=(target+start)//2
            if middle in (target,start):return False
            if not attempt(middle,depth+1):return False
            return attempt(target,depth+1)

    for index,target in enumerate(targets):
        if not attempt(target):
            data['unresolved_requested_stations']=targets[index:];break
    data['finished']=True;save();return data


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed-results',type=Path,required=True);p.add_argument('--solver',type=Path,required=True)
    p.add_argument('--label',required=True);p.add_argument('--station',type=int,default=312)
    p.add_argument('--workers',type=int,default=2);p.add_argument('--omega-max',type=float,default=400)
    p.add_argument('--signed-omega',action='store_true',help='Continue through zero; negative-frequency points remain actual solves')
    p.add_argument('--station-predictor',action='store_true',help='Use reviewed optional station predictor; solver v10 or later')
    p.add_argument('--fd-h',type=float,default=.02);p.add_argument('--max-iterations',type=int,default=15)
    p.add_argument('--beta-max',type=float,default=200);p.add_argument('--execute',action='store_true')
    p.add_argument('--directions',nargs='+',choices=['upstream','downstream'],default=['upstream','downstream'])
    p.add_argument('--source',nargs='+');a=p.parse_args();a.solver=a.solver.resolve();snapshot()
    d=json.loads(a.seed_results.read_text());seed_records=[]
    for record in d['results']:
        result=record['result'];source=record['job']['source']
        if result.get('status')!='candidate_converged' or record['job']['station']!=a.station:continue
        if a.source and source not in a.source:continue
        point=result['point']
        if any(abs(point['omega_r']-r['point']['omega_r'])<1e-5 and abs(point['beta_r']-r['point']['beta_r'])<1e-4
               and abs(point['alpha_r']-r['point']['alpha_r'])<1e-4 for r in seed_records):continue
        seed_records.append(dict(source=source,point=point,origin_job=record['job']['name']))
    grid=json.loads((ROOT/'outputs/seed_catalog.json').read_text())['grids']['main_841']['x_nondim']
    stations=coarse_stations(grid);jobs=[]
    for record in seed_records:
        for direction in a.directions:
            targets=sorted([i for i in stations if i<a.station],reverse=True) if direction=='upstream' else sorted(i for i in stations if i>a.station)
            jobs.append(dict(name=f'{a.label}_{record["source"]}_{direction}',seed_station=a.station,
                seed=record['point'],seed_provenance=record,stations=targets,solver=str(a.solver),
                solver_sha256=sha(a.solver),station_predictor=a.station_predictor,
                fd_h=a.fd_h,max_iterations=a.max_iterations,
                omega_min=-a.omega_max if a.signed_omega else 0,
                omega_max=a.omega_max,beta_max=a.beta_max))
    plan=ROOT/f'outputs/{a.label}_plan.json';content=dict(jobs=jobs,stations=stations,
        seed_input=str(a.seed_results),seed_input_sha256=sha(a.seed_results))
    if plan.exists():assert json.loads(plan.read_text())==content
    else:plan.write_text(json.dumps(content,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(chains=len(jobs),stations=len(stations),plan=str(plan))),flush=True)
    if a.execute:
        with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
            for result in pool.map(chain,jobs):
                print(json.dumps(dict(chain=result['spec']['name'],accepted=len(result['records']),failures=len(result['failures']),
                    unresolved=len(result.get('unresolved_requested_stations',[])))),flush=True)
        snapshot()


if __name__=='__main__':main()
