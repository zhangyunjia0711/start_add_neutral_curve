"""Actual clustered-spectrum regressions for opt-in discrete station prediction."""
import concurrent.futures, datetime, json, math, pathlib, re
from full_driver import ROOT, run, table, sha, snapshot

V8=ROOT/'outputs/bin/START_LST_wall_local_v8'
V9=ROOT/'outputs/bin/START_LST_wall_local_v9'
V8HASH='4efa0d3531d55a20060c91da6310bcd4e971eddc98af263947e3d0f54c0a2a8a'
V9HASH='31e52b504fecb0a5ca75d5fe6eec3764d8d149b83f63e80763d6e2bb4b692cea'
CASES=[dict(label='M1p',start=204,target=203,seed=dict(omega_r=35.08638891897888,beta_r=20.84934718914317,
       alpha_r=59.52502181299014,alpha_i=1.379992911708211e-9)),
       dict(label='M3near',start=56,target=55,seed=dict(omega_r=0.,beta_r=-7.,
       alpha_r=76.57982977151184,alpha_i=.1189407621523704))]
FLOAT=r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?'
MAC=re.compile(r'MAC best index=\s*(\d+)\s+eigen=\s*\(\s*('+FLOAT+r')\s*,\s*('+FLOAT+r')\s*\)\s+overlap=\s*('+FLOAT+r')\s+residual=\s*('+FLOAT+r')\s+jump=\s*('+FLOAT+r')')

def path(case,binary,digest,enabled,suffix='forward'):
    seed=case['seed'];b=seed['beta_r'];start=case['start'];target=case['target']
    spec=dict(name=f'station_predictor_{binary.name[-2:]}_{case["label"]}_{suffix}_i{start}_{target}_{"on" if enabled else "off"}',
        station=target,seed_station=start,mode='M1',seed=seed,seed_shift=1e-4,
        solver=str(binary),solver_sha256=digest,coordinates='wall_local',task=33,
        outer_parameter='iloc',inner_parameter='beta',outer_min=min(start,target),outer_max=max(start,target),
        outer_points=2,inner_min=b,inner_max=b+1e-8,inner_points=2,reverse_outer=target<start,
        anchor_row_seed=True,sample_only=True,step=.2,jump=4,adaptive=True,vector_tracking=True,
        minimum_overlap=.95,minimum_step=1e-6,retry_eigenpairs=8,station_predictor=enabled,
        options=['-pep_tol','1e-12','-pep_max_it','1000','-pep_nev','1','-pep_ncv','8',
                 '-st_pc_factor_mat_solver_type','superlu','-malloc_debug','-malloc_dump','-options_left'],timeout=240)
    summary=run(spec);directory=ROOT/'runs'/spec['name']
    field=table(directory/'out/neutral_field.plt');log=(directory/'solver.log').read_text(errors='replace')
    neutral=(directory/'out/neutral_run.log').read_text(errors='replace')
    matches=[];current=None
    for line in log.splitlines():
        match=MAC.search(line)
        if match:
            v=[float(x.replace('D','E')) for x in match.groups()[1:]]
            current=dict(index=int(match[1]),alpha_r=v[0],alpha_i=v[1],mac=v[2],residual=v[3],jump=v[4])
        accepted=re.search(r'MAC accepted=\s*([TF])',line)
        if accepted and current is not None:
            current['accepted']=accepted[1]=='T';matches.append(current);current=None
    accepted=[v for v in matches if v['accepted']]
    target_points=[v for v in field if int(v['iloc'])==target and v['status']==0]
    source_points=[v for v in field if int(v['iloc'])==start and v['status']==0]
    point=min(target_points,key=lambda v:abs(v['beta_r']-b)) if target_points else None
    result=dict(spec=spec,summary=summary,field=field,target_points=target_points,source_points=source_points,point=point,
        accepted=len(field)==4 and len(target_points)==2 and summary['rejected']==0 and all(v['status']==0 for v in field),
        min_accepted_mac=min((v['mac'] for v in accepted),default=None),
        max_accepted_residual=max((v['residual'] for v in accepted),default=None),
        max_accepted_jump=max((v['jump'] for v in accepted),default=None),mac_records=matches,
        malloc_dump_bytes=sum(map(int,re.findall(r'\[ 0\] (\d+) bytes',log))),
        predictor_probe_lines=[line for line in neutral.splitlines() if 'station_predictor_probe ' in line],
        predictor_history_lines=[line for line in neutral.splitlines() if 'station_predictor_history' in line],
        predictor_target_lines=[line for line in neutral.splitlines() if 'station_predictor_target' in line],
        fallback_lines=[line for line in neutral.splitlines() if 'station_predictor_fallback' in line])
    assert summary['petsc_errors']==0 and summary['unused_options']==0 and not summary['coordinate_metadata_missing']
    # Full START retains problem objects at process exit even in the frozen v8
    # baseline. Record this separately; it is not a probe-copy leak assertion.
    (directory/'station_predictor_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

def compare_default(a,b):
    assert len(a['field'])==len(b['field'])
    max_delta=0.
    for x,y in zip(a['field'],b['field']):
        assert x['status']==y['status'] and x['iloc']==y['iloc'] and x['beta_r']==y['beta_r']
        if x['status']==0:
            max_delta=max(max_delta,abs(complex(x['alpha_r'],x['alpha_i'])-complex(y['alpha_r'],y['alpha_i'])))
    assert max_delta<1e-8 and a['summary']['solves']==b['summary']['solves']
    return max_delta

def main():
    assert sha(V8)==V8HASH and sha(V9)==V9HASH;snapshot()
    results={};jobs=[(case,binary,digest,enabled) for case in CASES for binary,digest,enabled in
                    [(V8,V8HASH,False),(V9,V9HASH,False),(V9,V9HASH,True)]]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        for args,result in zip(jobs,pool.map(lambda args:path(*args),jobs)):
            case,binary,digest,enabled=args
            results.setdefault(case['label'],{})[f'{binary.name[-2:]}_{"on" if enabled else "off"}']=result
    output=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),cases=results,
                expected_v8_failures=['M1p i204→203','M3near i56→55'],qualification='Finite local spatial/MAC regression, not a full spectrum claim.')
    passes=[]
    for case in CASES:
        records=results[case['label']]
        records['default_max_delta_alpha']=compare_default(records['v8_off'],records['v9_off'])
        forward=records['v9_on']
        if forward['accepted']:
            point=forward['point'];seed={k:point[k] for k in ['omega_r','beta_r','alpha_r','alpha_i']}
            reverse=path(dict(label=case['label'],start=case['target'],target=case['start'],seed=seed),V9,V9HASH,True,'reverse')
            records['reverse']=reverse
            if reverse['accepted']:
                actual=reverse['point'];original=case['seed']
                records['return_delta_alpha']=abs(complex(actual['alpha_r'],actual['alpha_i'])-complex(original['alpha_r'],original['alpha_i']))
        passes.append(not records['v8_off']['accepted'] and not records['v9_off']['accepted'] and forward['accepted'] and
                      records.get('reverse',{}).get('accepted',False) and records.get('return_delta_alpha',math.inf)<1e-7)
    output['passed']=all(passes)
    p=ROOT/'outputs/station_predictor_v9_validation.json';assert not p.exists();p.write_text(json.dumps(output,indent=2)+'\n')
    snapshot();assert sha(V8)==V8HASH and sha(V9)==V9HASH
    print('STATION PREDICTOR '+json.dumps(dict(passed=output['passed'],cases={label:dict(enabled_accepted=v['v9_on']['accepted'],
        point=v['v9_on']['point'],min_mac=v['v9_on']['min_accepted_mac'],max_residual=v['v9_on']['max_accepted_residual'],
        probes=len(v['v9_on']['predictor_probe_lines']),fallbacks=len(v['v9_on']['fallback_lines']),
        return_delta_alpha=v.get('return_delta_alpha'),default_delta_alpha=v['default_max_delta_alpha']) for label,v in results.items()})),flush=True)

if __name__=='__main__':main()
