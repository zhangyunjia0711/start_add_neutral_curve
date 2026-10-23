"""Finite single-family M3near neutral-envelope search, immutable solver only."""
import argparse,datetime,hashlib,json,math
from pathlib import Path
import full_driver as driver
import fd_envelope as fd
ROOT=driver.ROOT
BIN=ROOT/'outputs/bin/START_LST_wall_local_v7'
HASH='7aa6016570f406587156f229107c0e5e4f76c8a85ba89dda267039ff1d688d79'
STATE=ROOT/'outputs/local_mode_audit_M3near_candidate.json'
LOG=ROOT/'outputs/local_mode_audit_工作日志_2026-10-22.md'

def stamp():return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
def log(message):
    with LOG.open('a',encoding='utf-8') as f:f.write('\n'+stamp()+' '+message+'\n')

def scan(seed,omega,lo,hi,db,label,reverse=False):
    spec=dict(station=119,mode='M1',coordinates='wall_local',
              seed={k:seed[k] for k in ('omega_r','beta_r','alpha_r','alpha_i')},seed_shift=1e-4,
              solver=str(BIN),solver_sha256=HASH,task=33,sample_only=True,
              omega_min=omega,omega_max=omega,omega_points=1,beta_min=lo,beta_max=hi,
              beta_points=math.ceil((hi-lo)/db-1e-12)+1,reverse_inner=reverse,
              anchor_row_seed=True,adaptive=True,vector_tracking=True,minimum_overlap=.95,
              minimum_step=1e-6,retry_eigenpairs=8,step=.4,jump=2,
              options=['-pep_tol','1e-12','-pep_max_it','1000','-pep_nev','1','-pep_ncv','8',
                       '-st_pc_factor_mat_solver_type','superlu'],timeout=240)
    digest=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()[:10]
    spec['name']='m3near_candidate_'+label+'_'+digest
    summary=driver.run(spec);path=ROOT/'runs'/spec['name']/'out/neutral_field.plt';rows=driver.table(path)
    valid=[dict(r,source_path=str(path),source_run=spec['name'],family='M3near',physical_role='single-family near-wall candidate, not global envelope') for r in rows if r['status']==0]
    valid.sort(key=lambda r:r['beta_r'])
    local=[]
    for i in range(1,len(valid)-1):
        a,b,c=valid[i-1:i+2]
        # A rejected/absent grid cell must not be bridged in peak detection.
        if max(b['beta_r']-a['beta_r'],c['beta_r']-b['beta_r'])>db*1.01:continue
        if b['growth']>a['growth'] and b['growth']>=c['growth']:
            local.append(dict(point=b,left=a,right=c))
    return dict(label=label,summary=summary,valid=valid,local_maxima=local,
                winner=max(valid,key=lambda r:r['growth']) if valid else None,
                all_grid_cells_accepted=bool(len(rows)==spec['beta_points'] and all(r['status']==0 for r in rows) and not summary['petsc_errors'] and not summary['coordinate_metadata_missing']))

def save(state):
    state['updated']=stamp();STATE.write_text(json.dumps(state,indent=2,allow_nan=False)+'\n')

def stationary(seed,omega,beta,label,h=.02,max_iterations=8):
    result=dict(status='unresolved',station=119,omega=omega,trace=[],solver=str(BIN),solver_sha256=HASH,
                qualification='One tracked M3near local beta maximum, not global competition proof.')
    current=dict(seed);b=beta
    try:
        for iteration in range(max_iterations):
            point=fd.stencil(119,current,omega,b,BIN,'wall_local',h=h,label='m3near_peak_'+label)
            result['trace'].append(point);gb=point['F'][1];curvature=point['J'][1][1]
            if curvature>=0:raise RuntimeError('Not a strict beta maximum')
            current=point['center']
            if abs(gb)<1e-7:
                result.update(status='local_maximum',point=dict(current,source_run=point['name'],source_path=str(ROOT/'runs'/point['name']/'out/neutral_field.plt'),
                              family='M3near',physical_role='single-family local beta maximum, not final global validation'),derivative=gb,curvature=curvature)
                break
            delta=max(-2.,min(2.,-gb/curvature));b+=delta
            if not -200<b<200:raise RuntimeError('stationary beta step leaves declared ±200 domain')
        else:result['reason']='iteration limit'
    except Exception as e:result['reason']=repr(e)
    return result

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--omega',type=float,default=0.);parser.add_argument('--seed',type=Path)
    parser.add_argument('--lo',type=float,default=-30);parser.add_argument('--hi',type=float,default=30)
    parser.add_argument('--db',type=float,default=1);parser.add_argument('--label',default='w0_b30')
    parser.add_argument('--refine',action='store_true');parser.add_argument('--beta',type=float)
    parser.add_argument('--stationary',action='store_true')
    parser.add_argument('--max-iterations',type=int,default=24)
    args=parser.parse_args();assert driver.sha(BIN)==HASH
    seed=json.loads((args.seed or ROOT/'outputs/local_mode_audit_nearwall_seed.json').read_text())
    if 'point' in seed:seed=seed['point']
    state=json.loads(STATE.read_text()) if STATE.exists() else dict(created=stamp(),solver=str(BIN),solver_sha256=HASH,scans=[],refinements=[],qualification='Finite M3near locally tracked candidate only; does not replace multi-family/max-beta envelope auditing.')
    log('M3near可变beta候选阶段：'+json.dumps(vars(args),default=str,ensure_ascii=False))
    if args.stationary:
        assert args.beta is not None
        result=stationary(seed,args.omega,args.beta,args.label)
        state.setdefault('stationary',[]).append(result);save(state)
        out=ROOT/f'outputs/local_mode_audit_M3near_peak_{args.label}.json';out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        print('M3NEAR STATIONARY '+json.dumps({k:v for k,v in result.items() if k!='trace'},indent=2),flush=True)
    elif args.refine:
        assert args.beta is not None
        result=fd.refine(119,seed,args.omega,args.beta,BIN,'wall_local',label='m3near_fd_'+args.label,
                         max_iterations=args.max_iterations,omega_bounds=(0,400),beta_bounds=(-200,200),h=.02)
        state['refinements'].append(result);save(state)
        print('M3NEAR REFINE '+json.dumps({k:v for k,v in result.items() if k not in ['trace','verification']},indent=2),flush=True)
    else:
        beta0=seed['beta_r']
        assert args.lo<=beta0<=args.hi
        for lo,hi,reverse,side in [(args.lo,beta0,True,'left'),(beta0,args.hi,False,'right')]:
            if hi-lo<1e-10:continue
            result=scan(seed,args.omega,lo,hi,args.db,args.label+'_'+side,reverse)
            state['scans'].append(result);save(state)
            print('M3NEAR SCAN '+json.dumps(dict(label=result['label'],winner=result['winner'],local_maxima=result['local_maxima'],all_grid_cells_accepted=result['all_grid_cells_accepted'])),flush=True)
    assert driver.sha(BIN)==HASH

if __name__=='__main__':main()
