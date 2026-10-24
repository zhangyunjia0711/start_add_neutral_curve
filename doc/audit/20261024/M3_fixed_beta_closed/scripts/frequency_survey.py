"""Bounded signed-frequency scans from actual fixed-beta M3 native origins."""
import argparse,concurrent.futures,datetime,hashlib,json
from common import ROOT,PREVIOUS,SOLVER,SOLVER_HASH,OPTIONS,run,table,sha,snapshot

def stamp():return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()

def main():
    p=argparse.ArgumentParser();p.add_argument('--stations',type=int,nargs='+',default=[119,321,600,841])
    p.add_argument('--cap',type=float,default=40);p.add_argument('--workers',type=int,default=4)
    p.add_argument('--spacing',type=float,default=.5);p.add_argument('--timeout',type=float,default=180)
    p.add_argument('--label',default='pilot');p.add_argument('--vectors',action='store_true');a=p.parse_args()
    seeds_path=ROOT/'outputs/prior_axis_seeds.json';seeds={int(r['iloc']):r for r in json.loads(seeds_path.read_text())['records']}
    guard=snapshot();jobs=[]
    for i in a.stations:
        for sign in [-1,1]:
            origin={k:float(seeds[i][k]) for k in ['omega_r','beta_r','alpha_r','alpha_i']}
            assert origin['omega_r']==0 and origin['beta_r']==-7
            options=OPTIONS+(['-pep_view_vectors','binary:out/all_vectors.petsc:default:append'] if a.vectors else [])
            spec=dict(station=i,mode='M3',seed=origin,seed_shift=1e-4,solver=str(SOLVER),solver_sha256=SOLVER_HASH,
                coordinates='wall_local',task=33,outer_parameter='iloc',outer_min=i,outer_max=i,outer_points=1,
                inner_parameter='omega',inner_min=-a.cap if sign<0 else 0,inner_max=0 if sign<0 else a.cap,
                inner_points=int(a.cap/a.spacing)+1,reverse_inner=sign<0,anchor_row_seed=True,
                sample_only=True,step=.2,jump=4,adaptive=True,vector_tracking=True,minimum_overlap=.95,
                minimum_step=1e-6,retry_eigenpairs=8,growth_tol=1e-9,parameter_tol=1e-8,
                options=options,timeout=a.timeout)
            token=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()[:10]
            spec['name']=f'M3b7_{a.label}_i{i:03d}_{"neg" if sign<0 else "pos"}_{token}';jobs.append(spec)
    out=ROOT/'outputs'/('frequency_'+a.label+'.json')
    result=dict(created=stamp(),solver=str(SOLVER),solver_sha256=SOLVER_HASH,prior_seeds_sha256=sha(seeds_path),
        jobs=jobs,results=[],original_inputs=guard,finished=False,qualification='Fixed-beta followed M3 sheet; signed frequencies are actual solves.')
    def save():
        temp=out.with_suffix('.tmp');temp.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');temp.replace(out)
    def solve(spec):
        summary=run(spec);d=ROOT/'runs'/spec['name'];roots=table(d/'out/neutral_roots.plt');field=table(d/'out/neutral_field.plt')
        complete=summary['terminal'] in ['converged','no bracketed root'] and not any(summary[k] for k in ['petsc_errors','unused_options','failed_field_records','coordinate_metadata_missing']) and summary['exit']!=124
        assert all(q['beta_r']==-7 and q['beta_i']==0 for q in roots+field)
        ordered=sorted(field,key=lambda q:q['omega_r'])
        brackets=[dict(left=x,right=y) for x,y in zip(ordered,ordered[1:]) if x['growth']*y['growth']<0]
        return dict(name=spec['name'],station=spec['station'],complete=complete,summary=summary,roots=roots,
            sampled_frequency=[min(q['omega_r'] for q in field),max(q['omega_r'] for q in field)] if field else None,
            field_peak=max(field,key=lambda q:q['growth']) if field else None,
            endpoints=[field[0],field[-1]] if field else [],root_file=str(d/'out/neutral_roots.plt'),brackets=brackets)
    save()
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
        for d in pool.map(solve,jobs):
            result['results'].append(d);result['updated']=stamp();save()
            print(json.dumps(dict(name=d['name'],complete=d['complete'],roots=d['roots'],peak=d['field_peak'])),flush=True)
    result['finished']=True;result['updated']=stamp();result['original_inputs_final']=snapshot();save()

if __name__=='__main__':main()
