"""True-station, fixed-beta neutral continuation with safeguarded frequency correction."""
import argparse,concurrent.futures,datetime,hashlib,json,math
import numpy as np
from common import ROOT,PREVIOUS,SOLVER,SOLVER_HASH,OPTIONS,run,table

def stamp():return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
def seed(q):return {k:float(q[k]) for k in ['omega_r','beta_r','alpha_r','alpha_i']}

def field(station,q,low,high,n,label,origin_station=None):
    source=station if origin_station is None else origin_station
    spec=dict(station=station,seed_station=source,mode='M3',seed=seed(q),seed_shift=1e-5,
        solver=str(SOLVER),solver_sha256=SOLVER_HASH,coordinates='wall_local',task=33,
        outer_parameter='iloc',outer_min=min(source,station),outer_max=max(source,station),outer_points=1 if source==station else 2,
        inner_parameter='omega',inner_min=low,inner_max=high,inner_points=n,
        reverse_outer=station<source,sample_only=True,anchor_row_seed=True,station_predictor=True,
        step=.1,jump=4,adaptive=True,vector_tracking=True,minimum_overlap=.95,
        minimum_step=1e-6,retry_eigenpairs=8,options=OPTIONS,timeout=90)
    token=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()[:10];spec['name']=label+'_'+token
    s=run(spec);values=table(ROOT/'runs'/spec['name']/'out/neutral_field.plt')
    if s['exit']==124 or any(s[k] for k in ['petsc_errors','unused_options','failed_field_records','coordinate_metadata_missing']):
        raise RuntimeError('Unresolved field '+spec['name'])
    values=[v for v in values if int(v['iloc'])==station and v['status']==0]
    if len(values)!=n:raise RuntimeError('Missing target samples '+spec['name'])
    assert all(v['beta_r']==-7 for v in values)
    return sorted(values,key=lambda v:v['omega_r']),spec['name']

def frequency_root(station,actual,guess,orientation,label):
    state=dict(actual);w=float(guess);history=[];h=.0025
    result=dict(status='unknown',station=station,trace=history,expected_orientation=orientation)
    try:
        for iteration in range(24):
            values,name=field(station,state,w-h,w+h,3,label+'_stencil')
            left,center,right=values;g=center['growth'];deriv=(right['growth']-left['growth'])/(2*h)
            row=dict(run=name,point=center,h=h,g=g,g_omega=deriv,samples=values);history.append(row)
            if not math.isfinite(deriv) or abs(deriv)<1e-10:raise RuntimeError('Frequency derivative near zero; endpoint treatment needed')
            if deriv*orientation<=0:raise RuntimeError('Frequency-boundary orientation changed; do not cross a turning point')
            if abs(g)<1e-9:
                fine,nf=field(station,center,w-h/2,w+h/2,3,label+'_verify')
                dw=(fine[2]['growth']-fine[0]['growth'])/h
                if dw*orientation>0 and abs(fine[1]['growth'])<1e-8 and abs(fine[1]['growth']/dw)<1e-6:
                    result.update(status='converged',point=fine[1],verification=dict(run=nf,h=h/2,
                        g_omega=dw,samples=fine,omega_correction=-fine[1]['growth']/dw));return result
                h/=2;state=fine[1];continue
            step=max(-1.,min(1.,-g/deriv));accepted=False;trials=[]
            for k in range(7):
                candidate=w+step*2**(-k)
                if not -100<candidate<40:continue
                vv,nn=field(station,center,candidate-1e-9,candidate+1e-9,3,label+'_line')
                trial=vv[1];trials.append(dict(run=nn,point=trial))
                if abs(trial['growth'])<abs(g):state=trial;w=candidate;accepted=True;break
            row['line_search']=trials
            if not accepted:raise RuntimeError('No decreasing frequency correction; near turning point or wrong prediction')
        raise RuntimeError('Frequency corrector iteration limit')
    except Exception as error:result['reason']=repr(error);return result

def main():
    p=argparse.ArgumentParser();p.add_argument('--label',default='primary');p.add_argument('--seed-file',default='frequency_pilot_fields_roots.json')
    p.add_argument('--seed-station',type=int,default=119);p.add_argument('--workers',type=int,default=4)
    a=p.parse_args();doc=json.loads((ROOT/'outputs'/a.seed_file).read_text())
    candidates=[r for r in doc['results'] if r['station']==a.seed_station and r['passed']]
    grid=json.loads((PREVIOUS/'outputs/seed_catalog.json').read_text())['grids']['main_841']['x_nondim']
    targets=[1,8,16,20,24,28,32,36,40,48,56,64,70,80,90,100,110,119,140,160,190,210,240,241,260,280,300,321,350,380,400,420,440,454,480,520,550,580,600,625,650,680,700,712,730,760,800,841]
    jobs=[]
    for j,c in enumerate(candidates):
        for direction in [-1,1]:
            seq=sorted([i for i in targets if (i-a.seed_station)*direction>0],reverse=direction<0)
            jobs.append(dict(name=f'M3b7_{a.label}_{c["boundary"]}{j}_{"up" if direction<0 else "down"}',
                source_result=c,seed_station=a.seed_station,seed=c['roots'][0],targets=seq,
                orientation=1 if c['boundary']=='lower' else -1))
    def chain(job):
        out=ROOT/'outputs'/(job['name']+'.json')
        if out.exists():return json.loads(out.read_text())
        state=dict(created=stamp(),job=job,records=[],failures=[],finished=False)
        history=[dict(station=job['seed_station'],point=job['seed'])]
        def save():
            temp=out.with_suffix('.tmp');temp.write_text(json.dumps(state,indent=2,allow_nan=False)+'\n');temp.replace(out)
        def attempt(target,depth=0):
            old=history[-1];q=old['point'];w=q['omega_r']
            if len(history)>1:
                p0=history[-2];w+=(grid[target-1]-grid[old['station']-1])/(grid[old['station']-1]-grid[p0['station']-1])*(q['omega_r']-p0['point']['omega_r'])
            record=dict(station=target,origin=old['station'],guess=w,depth=depth)
            try:
                transported,name=field(target,q,q['omega_r'],q['omega_r']+1e-9,2,job['name']+'_transfer',old['station'])
                record['transfer_run']=name
                result=frequency_root(target,transported[0],w,job['orientation'],job['name']+f'_i{target:03d}')
                record['result']=result
                if result['status']!='converged':raise RuntimeError(result.get('reason'))
                record['point']=result['point'];state['records'].append(record);history.append(record)
                save();print('TRACE '+json.dumps(dict(chain=job['name'],station=target,omega=record['point']['omega_r'])),flush=True);return True
            except Exception as error:
                record['reason']=repr(error);state['failures'].append(record);save();start=old['station']
                if abs(target-start)<=1:return False
                mid=(target+start)//2
                return attempt(mid,depth+1) and attempt(target,depth+1)
        save()
        for i in job['targets']:
            if not attempt(i):state['unresolved_from_station']=i;break
        state.update(finished=True,updated=stamp());save();return state
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
        for s in pool.map(chain,jobs):print('CHAIN_END',s['job']['name'],len(s['records']),len(s['failures']),flush=True)

if __name__=='__main__':main()
