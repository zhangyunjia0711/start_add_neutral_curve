"""Finite-domain, multi-start survey. Missing/failed cells remain unknown.

The six starts close the three archived seed families under joint conjugacy.
They are search provenance labels, not a claim of complete spectral enumeration.
"""
import argparse, concurrent.futures, datetime, json, math, pathlib, time
from full_driver import ROOT, run, table, snapshot, sha
from seed_catalog import load, get_candidates

FAMILIES=['M1','M2','M3']
COARSE=sorted(set(range(1,842,20))|{119,241,312,317,321,371,381,421,481})
PILOT=[1,61,119,241,312,381,601,841]

def seed_for(catalogue, station, family, conjugate=False,coordinates='global'):
    axial_path=ROOT/'outputs/axial_seeds.json'
    records=json.loads(axial_path.read_text()).get('records',[]) if axial_path.exists() else []
    vector_path=ROOT/'outputs/axial_vectors.json'
    if vector_path.exists():records+=json.loads(vector_path.read_text()).get('records',[])
    if coordinates=='wall_local':
        records=[]
        for path in [ROOT/'outputs/axial_wall_local.json',ROOT/'outputs/local_seed_bootstrap.json']:
            if path.exists():records+=json.loads(path.read_text()).get('records',[])
    usable=[v for v in records if v['family']==family and v.get('status')==0]
    defaults={'M1':(8.17586,-.98903),'M2':(15,15),'M3':(0,-7)}
    om,be=defaults[family]
    archived=get_candidates(catalogue,station,'fixed_'+family,omega=om,beta=be,limit=1)[0]
    target_x=catalogue['grids']['main_841']['x_nondim'][station-1]
    best_axial=min(usable,key=lambda v:abs(v['x']-target_x)) if usable else None
    if best_axial and abs(best_axial['x']-target_x)<=abs(archived['x']-target_x)+1e-12:
        original=best_axial
        source=dict(kind='actual_axial_seed',**original)
    else:
        original=archived
        source=dict(kind='archived_nearest_initial_guess',**original)
    seed={k:original[k] for k in ['omega_r','beta_r','alpha_r','alpha_i']}
    if conjugate:
        for k in ['omega_r','beta_r','alpha_r']:seed[k]=-seed[k]
    source['derived_joint_conjugate']=conjugate
    source['requested_coordinates']=coordinates
    source['source_is_new_model_solution']=coordinates=='global' or source['kind']=='actual_axial_seed'
    return seed,source

def specs(catalogue, stations, solver, label='survey', bmax=20, wmax=40, dw=4, db=1, step=.4,coordinates='wall_local',wmin=0):
    result=[]
    for i in stations:
        for family in FAMILIES:
            for conjugate in [False,True]:
                tag=family+('c' if conjugate else 'p')
                seed,source=seed_for(catalogue,i,family,conjugate,coordinates)
                result.append(dict(name=f'{label}_{tag}_i{i:03d}',family=family,conjugate=conjugate,
                    station=i,mode='M1',seed=seed,seed_provenance=source,solver=str(solver),solver_sha256=sha(solver),
                    task=33,coordinates=coordinates,sample_only=True,omega_min=wmin,omega_max=wmax,omega_points=round((wmax-wmin)/dw)+1,
                    beta_min=-bmax,beta_max=bmax,beta_points=round(2*bmax/db)+1,step=step,
                    adaptive=True,vector_tracking=True,minimum_overlap=.95,minimum_step=1e-5,
                    options=['-pep_tol','1e-12','-pep_max_it','1000','-pep_nev','8','-pep_ncv','24','-st_pc_factor_mat_solver_type','superlu'],
                    timeout=7200))
    return result

def check_continuity(field):
    # Diagnostic heuristic only. A flag is not automatically deleted/relabelled.
    bad=[]
    groups={}
    for v in field:
        if v['status']==0:groups.setdefault(v['row'],[]).append(v)
    for row,values in groups.items():
        values.sort(key=lambda v:v['beta_r'])
        for a,b in zip(values,values[1:]):
            db=abs(a['beta_r']-b['beta_r']); da=abs(complex(a['alpha_r']-b['alpha_r'],a['alpha_i']-b['alpha_i']))
            if db<.02 and da>max(1e-4,30*db):
                bad.append(dict(row=row,beta_a=a['beta_r'],beta_b=b['beta_r'],delta_alpha=da))
    return bad

def summarize(result):
    d=ROOT/'runs'/result['spec']['name'];field=table(d/'out/neutral_field.plt')
    anomalies=check_continuity(field);rows=[]
    for j in range(1,result['spec']['omega_points']+1):
        f=[v for v in field if v['row']==j]
        valid=[v for v in f if v['status']==0 and all(math.isfinite(v[k]) for k in ['alpha_r','alpha_i','growth'])]
        if not valid:
            rows.append(dict(row=j,status='unknown',omega=result['spec']['omega_min']+(j-1)*(result['spec']['omega_max']-result['spec']['omega_min'])/max(1,result['spec']['omega_points']-1),samples=len(f)))
            continue
        winner=max(valid,key=lambda v:v['growth'])
        bad=sum(v['status']!=0 for v in f)
        rows.append(dict(row=j,omega=winner['omega_r'],winner=winner,samples=len(f),failed=bad,
            status='unknown' if bad or result.get('petsc_errors',0) else 'sampled',
            continuity_flags=sum(v['row']==j for v in anomalies),
            beta_edge=abs(winner['beta_r'])>=result['spec']['beta_max']-.0001))
    return dict(name=result['spec']['name'],family=result['spec']['family'],conjugate=result['spec']['conjugate'],
        station=result['spec']['station'],summary_path=str(d/'summary.json'),rows=rows,
        terminal=result['terminal'],seconds=result['seconds'],petsc_errors=result.get('petsc_errors'),
        rejected=result['rejected'],continuity_flags=anomalies)

def main():
    p=argparse.ArgumentParser();p.add_argument('--phase',default='wallpilot')
    p.add_argument('--coordinates',choices=['global','wall_local'],default='wall_local')
    p.add_argument('--wmin',type=float,default=0);p.add_argument('--wmax',type=float,default=40)
    p.add_argument('--bmax',type=float,default=20);p.add_argument('--dw',type=float,default=4)
    p.add_argument('--db',type=float,default=1);p.add_argument('--step',type=float,default=.4)
    p.add_argument('--workers',type=int,default=6);p.add_argument('--solver',type=pathlib.Path,required=True)
    p.add_argument('--stations',nargs='*',type=int);a=p.parse_args()
    assert a.solver.exists();snapshot();catalogue=load()
    stations=a.stations or (PILOT if a.phase.endswith('pilot') else COARSE)
    plan=specs(catalogue,stations,a.solver,label=a.phase,coordinates=a.coordinates,wmin=a.wmin,wmax=a.wmax,bmax=a.bmax,dw=a.dw,db=a.db,step=a.step)
    (ROOT/f'outputs/{a.phase}_plan.json').write_text(json.dumps(plan,indent=2)+'\n')
    results=[];started=time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
        pending={pool.submit(run,s):s for s in plan}
        for f in concurrent.futures.as_completed(pending):
            spec=pending[f]
            try:results.append(summarize(f.result()))
            except Exception as e:
                results.append(dict(name=spec['name'],station=spec['station'],family=spec['family'],conjugate=spec['conjugate'],status='exception',exception=repr(e)))
            output=dict(phase=a.phase,coordinates=a.coordinates,elapsed_seconds=time.time()-started,completed=len(results),planned=len(plan),stations=stations,domain=dict(omega=[a.wmin,a.wmax],beta=[-a.bmax,a.bmax]),results=results)
            temp=ROOT/f'outputs/{a.phase}_progress.tmp';temp.write_text(json.dumps(output,indent=2)+'\n');temp.replace(ROOT/f'outputs/{a.phase}_progress.json')
            print('PROGRESS '+json.dumps({k:output[k] for k in ['phase','elapsed_seconds','completed','planned']}),flush=True)
    snapshot()

if __name__=='__main__':main()
