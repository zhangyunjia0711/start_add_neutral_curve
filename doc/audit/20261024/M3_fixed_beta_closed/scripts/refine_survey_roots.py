"""Native safeguarded fixed-beta roots, seeded by actual frequency brackets."""
import argparse,concurrent.futures,hashlib,json,datetime
from pathlib import Path
from common import ROOT,SOLVER,SOLVER_HASH,OPTIONS,run,table,sha

def solve(station,bracket,label):
    lo,hi=bracket['left'],bracket['right'];assert lo['growth']*hi['growth']<0
    seed={k:lo[k] for k in ['omega_r','beta_r','alpha_r','alpha_i']};assert seed['beta_r']==-7
    spec=dict(station=station,mode='M3',seed=seed,seed_shift=1e-5,solver=str(SOLVER),solver_sha256=SOLVER_HASH,
        coordinates='wall_local',task=33,outer_parameter='iloc',outer_min=station,outer_max=station,outer_points=1,
        inner_parameter='omega',inner_min=lo['omega_r'],inner_max=hi['omega_r'],inner_points=2,
        sample_only=False,step=.1,jump=4,adaptive=True,vector_tracking=True,minimum_overlap=.95,
        minimum_step=1e-6,retry_eigenpairs=8,growth_tol=1e-10,parameter_tol=1e-10,
        options=OPTIONS+['-pep_view_vectors','binary:out/all_vectors.petsc:default:append'],timeout=120)
    token=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()[:10]
    spec['name']=f'M3b7_root_{label}_i{station:03d}_{token}'
    s=run(spec);p=ROOT/'runs'/spec['name'];roots=table(p/'out/neutral_roots.plt')
    passed=s['terminal']=='converged' and not any(s[k] for k in ['petsc_errors','unused_options','failed_field_records','coordinate_metadata_missing']) and len(roots)==1
    if roots:assert all(q['beta_r']==-7 and abs(q['growth'])<1e-8 for q in roots)
    return dict(station=station,bracket=bracket,passed=passed,summary=s,roots=roots,
        boundary='lower' if lo['growth']<0 else 'upper',raw_vector_file=str(p/'out/all_vectors.petsc'))

def main():
    p=argparse.ArgumentParser();p.add_argument('--survey',type=Path,required=True);p.add_argument('--workers',type=int,default=4)
    a=p.parse_args();d=json.loads(a.survey.read_text());jobs=[]
    for r in d['results']:
        if r['complete']:
            for j,b in enumerate(r['brackets']):jobs.append((r['station'],b,Path(a.survey).stem+'_'+str(j)))
    result=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),
        survey_path=str(a.survey),survey_sha256=sha(a.survey),results=[],finished=False)
    out=ROOT/'outputs'/(a.survey.stem+'_roots.json')
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
        for v in pool.map(lambda t:solve(*t),jobs):
            result['results'].append(v);out.write_text(json.dumps(result,indent=2)+'\n')
            print('ROOT '+json.dumps(dict(station=v['station'],passed=v['passed'],boundary=v['boundary'],roots=v['roots'])),flush=True)
    result['finished']=True;out.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
