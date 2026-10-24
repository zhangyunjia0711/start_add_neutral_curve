"""Bounded signed-frequency exploration of the continued rough-case M3 spatial root."""
import argparse,concurrent.futures,datetime,json
import numpy as np
from common import ROOT,sha,snapshot
from fixed_beta_pencil import OperatorLinearFamily
from spatial_pencil import solve

def arm(direction,seed,label,limit):
    family=OperatorLinearFamily();station=seed['station'];p=family.native(station)
    q=np.load(seed['vector_file']);alpha=complex(seed['point']['alpha_r'],seed['point']['alpha_i']);omega=0.;step=.025
    out=ROOT/'outputs'/('%s_frequency_%s.json'%(label,'positive' if direction>0 else 'negative'))
    assert not out.exists();wave=ROOT/'work'/out.stem;wave.mkdir(exist_ok=True)
    result=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),station=station,
        direction=direction,source=seed['source_path'],records=[],roots=[],failures=[],finished=False)
    last=dict(point=seed['point'],vector_file=seed['vector_file'])
    while abs(omega)<limit-1e-12:
        target=direction*min(limit,abs(omega)+step)
        try:
            r,v=solve(p,target,alpha,q,.98);newalpha=complex(r['alpha_r'],r['alpha_i'])
            assert abs(newalpha-alpha)<1.,'Large alpha move'
            path=wave/('q%04d.npy'%len(result['records']));np.save(path,v)
            item=dict(point=r,vector_file=str(path),vector_sha256=sha(path));result['records'].append(item)
            if r['growth']*last['point']['growth']<0:
                a,b=sorted([last,item],key=lambda a:a['point']['omega_r']);trace=[]
                for k in range(60):
                    wa,wb=a['point']['omega_r'],b['point']['omega_r'];ga,gb=a['point']['growth'],b['point']['growth']
                    t=np.clip(-ga/(gb-ga),.1,.9);w=wa+t*(wb-wa);near=a if t<.5 else b
                    rr,vv=solve(p,w,complex(near['point']['alpha_r'],near['point']['alpha_i']),np.load(near['vector_file']),.97)
                    f=wave/('root_trial_%d_%d.npy'%(len(result['roots']),k));np.save(f,vv);entry=dict(point=rr,vector_file=str(f),vector_sha256=sha(f));trace.append(entry)
                    if abs(rr['growth'])<1e-10:
                        result['roots'].append(dict(point=rr,vector_file=str(f),trace=trace,bracket=[a,b]));print('FREQUENCY_ROOT',station,rr,flush=True);break
                    if rr['growth']*ga>0:a=entry
                    else:b=entry
                else:raise RuntimeError('Neutral root refinement did not converge')
            omega=target;alpha=newalpha;q=v;last=item;step=min(.1,step*1.2)
            if len(result['records'])%50==0:print('FREQUENCY',direction,omega,alpha,r['growth'],flush=True)
            out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        except Exception as error:
            result['failures'].append(dict(omega=omega,target=target,reason=repr(error)));step/=2
            out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
            if step<1e-6:result['unresolved_at_omega']=omega;break
    result.update(finished=True,original_inputs=snapshot());out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return result

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--seed',default='rough_physical_M3_seed.json');parser.add_argument('--label',default='physical_M3');parser.add_argument('--limit',type=float,default=15.);a=parser.parse_args()
    source=ROOT/'outputs'/a.seed;seed=json.loads(source.read_text());assert seed['finished'];seed['source_path']=str(source)
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        for r in pool.map(lambda sign:arm(sign,seed,a.label,a.limit),[-1,1]):print('SURVEY_FINISHED',r['direction'],len(r['records']),len(r['roots']),r.get('unresolved_at_omega'),flush=True)

if __name__=='__main__':main()
