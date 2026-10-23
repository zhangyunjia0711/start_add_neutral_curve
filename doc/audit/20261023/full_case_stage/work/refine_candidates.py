#!/usr/bin/env python3
"""Refine all observed source sign changes; candidates require later union audit.

Input is the read-only envelope_reduce schema. Every origin is an actual field
sample. Domain endpoints are retained as unresolved endpoint candidates instead
of being submitted to the interior two-equation Newton method.
"""
import argparse
import concurrent.futures
import hashlib
import json
from pathlib import Path
import time

from full_driver import ROOT,sha,snapshot
from fd_envelope import refine


def build_jobs(data,solver,label):
    peaks={v['peak_id']:v for v in data['beta_peak_candidates']}
    domains={v['phase']:v for v in data.get('inputs',[])}
    jobs=[];deferred=[]
    for bracket in data['neutral_brackets']:
        if bracket['scope']!='source':continue
        left,right=bracket['left'],bracket['right']
        if left['omega']==right['omega']:continue
        choices=[]
        for endpoint in (left,right):
            candidates=[peaks[i] for i in endpoint['beta_peak_candidate_ids']
                        if i in peaks and peaks[i]['is_local_max_on_available_neighbors']]
            interior=[v for v in candidates if v['kind']=='interior_local_maximum']
            # Retain distinct local maxima, including weaker positive lobes.
            interior.sort(key=lambda v:v['growth'],reverse=True)
            if interior:
                top=interior[0]['growth']
                choices.extend(v for v in interior if v['growth']>0 or v['growth']>=top-.1)
            for v in candidates:
                if v['is_domain_endpoint'] and (v['growth']>0 or v['seed']==endpoint.get('seed')):
                    deferred.append(dict(bracket_id=bracket['bracket_id'],station=bracket['station'],
                        source=bracket['source'],kind='beta_endpoint_needs_separate_treatment',peak=v))
        chosen=[]
        for peak in sorted(choices,key=lambda v:abs(v['growth'])):
            # Only deduplicate nearly identical beta guesses; both distant peaks
            # are retained even when one is weaker than the sampled winner.
            if any(abs(peak['beta']-q['beta'])<1.1 for q in chosen):continue
            chosen.append(peak)
        for peak in chosen:
            origin=peak['seed']
            fraction=left['known_G']/(left['known_G']-right['known_G'])
            w_linear=left['omega']+fraction*(right['omega']-left['omega'])
            job=dict(station=bracket['station'],source=bracket['source'],coordinates=data['coordinates'],
                bracket=bracket,origin_peak=peak,seed=origin,omega=origin['omega_r'],beta=peak['beta'],
                omega_linear_estimate=w_linear,
                solver=str(solver),solver_sha256=sha(solver),
                omega_bounds=(max(0,left['omega']-4),right['omega']+4),beta_bounds=(-80,80),
                qualification='Single continued mode interior candidate; no global union certification.')
            token=hashlib.sha256(json.dumps(job,sort_keys=True).encode()).hexdigest()[:12]
            job['name']=f'{label}_{bracket["source"]}_i{bracket["station"]:03d}_{token}'
            jobs.append(job)
    return jobs,deferred


def execute(job):
    path=ROOT/'outputs/refined_candidates'/f'{job["name"]}.json'
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():
        old=json.loads(path.read_text());assert old['job']==job;return old
    result=refine(job['station'],job['seed'],job['omega'],job['beta'],Path(job['solver']),
        job['coordinates'],label=job['name'],omega_bounds=job['omega_bounds'],beta_bounds=job['beta_bounds'])
    if result['status']=='candidate_converged':
        w=result['point']['omega_r'];b=job['bracket']
        result['inside_observed_bracket']=b['left']['omega']-1e-6<=w<=b['right']['omega']+1e-6
    output=dict(job=job,result=result)
    path.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(name=job['name'],status=result['status'],reason=result.get('reason'),
        point=result.get('point'),inside_observed_bracket=result.get('inside_observed_bracket'))),flush=True)
    return output


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('input',type=Path);p.add_argument('--solver',type=Path,required=True)
    p.add_argument('--label',required=True);p.add_argument('--workers',type=int,default=2)
    p.add_argument('--max-jobs',type=int);p.add_argument('--stations',type=int,nargs='+')
    p.add_argument('--execute',action='store_true');a=p.parse_args()
    snapshot();d=json.loads(a.input.read_text());jobs,deferred=build_jobs(d,a.solver,a.label)
    if a.stations:jobs=[j for j in jobs if j['station'] in a.stations]
    if a.max_jobs:jobs=jobs[:a.max_jobs]
    output=dict(input=str(a.input),input_sha256=sha(a.input),jobs=jobs,deferred_endpoints=deferred,
        qualification='Refinement only; missing competitors, domains and resolution checks remain separate.')
    plan=ROOT/f'outputs/{a.label}_plan.json'
    if plan.exists():assert json.loads(plan.read_text())==output,'Use a new label for changed inputs'
    else:plan.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(jobs=len(jobs),deferred_endpoints=len(deferred),plan=str(plan))),flush=True)
    if not a.execute:return
    results=[];start=time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
        pending={pool.submit(execute,j):j for j in jobs}
        for f in concurrent.futures.as_completed(pending):
            j=pending[f]
            try:result=f.result()
            except Exception as e:result=dict(job=j,result=dict(status='exception',reason=repr(e)))
            results.append(result)
            progress=dict(completed=len(results),planned=len(jobs),elapsed_seconds=time.time()-start,results=results)
            temp=ROOT/f'outputs/{a.label}_progress.tmp'
            temp.write_text(json.dumps(progress,indent=2,allow_nan=False)+'\n')
            temp.replace(ROOT/f'outputs/{a.label}_progress.json')
    snapshot()


if __name__=='__main__':main()
