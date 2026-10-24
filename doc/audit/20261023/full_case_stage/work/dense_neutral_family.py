"""Densify one verified neutral family on original x stations.

Frozen actual origins, true iloc/vector transport, bounded local correction,
and h/2 validation are required. An unresolved station never means stable.
This does not enumerate the complete spectrum or establish a causal family.
"""
import argparse, concurrent.futures, datetime, json, shutil, time
from pathlib import Path
from full_driver import ROOT, sha, snapshot
from trace_envelope_x import transfer
from fd_envelope import refine


def stamp():
    return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()


def write(path, data):
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');temp.replace(path)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',required=True)
    p.add_argument('--boundary',required=True,choices=['upper_local','lower_local'])
    p.add_argument('--label',required=True)
    p.add_argument('--workers',type=int,default=6)
    p.add_argument('--limit',type=int)
    a=p.parse_args();snapshot()
    area=ROOT/'outputs'/a.label;area.mkdir(exist_ok=True)
    (area/'stations').mkdir(exist_ok=True)
    frozen=area/'candidate_inputs.json'
    if not frozen.exists():shutil.copyfile(ROOT/'outputs/neutral_variable_beta_candidates.json',frozen)
    known={r['station']:r for r in json.loads(frozen.read_text())['unique_candidates_data']
           if r['source']==a.source and r['frequency_boundary']==a.boundary and r['omega_r']>=0}
    assert len(known)>=2,'Need independently verified neighbouring origins'
    solver=ROOT/'outputs/bin/START_LST_wall_local_v10'
    jobs=[]
    for station in range(min(known),max(known)+1):
        if station in known:continue
        neighbours=sorted(known,key=lambda i:abs(i-station))[:2]
        origins=[known[i] for i in neighbours if abs(i-station)<=20]
        jobs.append(dict(station=station,origins=origins))
    plan=dict(created=stamp(),source=a.source,boundary=a.boundary,
        solver=str(solver),solver_sha256=sha(solver),input_sha256=sha(frozen),
        station_range=[min(known),max(known)],jobs=jobs,
        policy=dict(fd_h=.00125,beta_halfwidth=.5,omega_tolerance=1e-4,
            beta_tolerance=1e-5,minimum_MAC=.95,maximum_origin_station_distance=20),
        qualification='Verified local family only; global competition and causal classification are separate.')
    planpath=area/'plan.json'
    if planpath.exists():
        old=json.loads(planpath.read_text())
        for k in plan:
            if k!='created':assert old[k]==plan[k],k
        plan=old
    else:write(planpath,plan)

    def solve(job):
        station=job['station'];path=area/'stations'/('i%03d.json'%station)
        if path.exists():return json.loads(path.read_text())
        results=[];failures=[];start=time.monotonic()
        for number,old in enumerate(job['origins']):
            label=a.label+'_i%03d_from%03d_p%d'%(station,old['station'],number)
            entry=dict(station=station,source=a.source,name=label,
                origin_station=old['station'],origin_result_file=old['result_file'],
                origin_verification_run=old['verification_run'])
            try:
                actual,name=transfer(old['station'],station,old,solver,label,True)
                entry['transport_run']=name
                result=refine(station,actual,actual['omega_r'],actual['beta_r'],solver,
                    label=label,h=.00125,max_iterations=25,omega_bounds=(0,80),
                    beta_bounds=(old['beta_r']-.5,old['beta_r']+.5))
                if result['status']=='candidate_converged':
                    actual_kind='upper_local' if result['verification']['J'][0][0]<0 else 'lower_local'
                    if actual_kind!=a.boundary:
                        # Preserve the actual other root, but never count it as this family.
                        failures.append(dict(job=entry,reason='Frequency-boundary orientation changed',actual_boundary=actual_kind))
                    else:
                        results.append(dict(job=entry,result=result));break
                else:failures.append(dict(job=entry,reason=result.get('reason','unresolved')))
                results.append(dict(job=entry,result=result))
            except Exception as error:failures.append(dict(job=entry,reason=repr(error)))
        accepted=[r for r in results if r['result'].get('status')=='candidate_converged'
            and (r['result']['verification']['J'][0][0]<0)==(a.boundary=='upper_local')]
        data=dict(station=station,finished=True,created=stamp(),seconds=time.monotonic()-start,
            status='actual_family_candidate' if accepted else 'unknown',results=results,failures=failures)
        write(path,data);return data

    active=jobs[:a.limit] if a.limit else jobs
    progress=dict(created=stamp(),plan_sha256=sha(planpath),total_jobs=len(jobs),
        scheduled_jobs=len(active),results=[],finished=False)
    write(area/'progress.json',progress)
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
        pending={pool.submit(solve,j):j for j in active}
        for f in concurrent.futures.as_completed(pending):
            d=f.result();short={k:d[k] for k in ['station','status','seconds']}
            progress['results'].append(short);progress['updated']=stamp()
            write(area/'progress.json',progress);print(json.dumps(short),flush=True)
    progress.update(finished=True,updated=stamp());write(area/'progress.json',progress);snapshot()


if __name__=='__main__':main()
