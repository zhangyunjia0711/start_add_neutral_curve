"""First fixed-station frequency folds, with actual full-matrix verification."""
import concurrent.futures,json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/neutral_m3_frequency_limits_20261027')
sys.path.insert(0,str(ROOT/'work'))
from common import CASES,snapshot,sha
from solver import Pencil,qep,mac
from fold_core import Equations,S,W
SOURCE=ROOT.parent/'neutral_xbeta_extended_frequencies_20261027'

def fold(case,i,a,b,f,q,label):
    p=Pencil(case,i,W*f);_,q0,_,p0=qep(p.spatial(b),complex(a),q,True);eq=Equations(p,q0,p0);z0=np.array([a,b,f])/S;scale=eq.scale(z0);trace=[];cache={}
    def ev(z):
        key=tuple(z)
        if key not in cache:cache.clear();cache[key]=eq.evaluate(z)
        return cache[key]
    def fun(z):
        r=ev(z);assert r['MAC']>.9
        trace.append(dict(alpha_r=float(z[0]*20),beta_r=float(z[1]*10),f_kHz=float(z[2]*100),scaled_F=(r['F']/scale).tolist(),right_residual=r['right_residual'],MAC=r['MAC']))
        return r['F']/scale
    lower=np.array([a-40,b-30,f-150])/S;upper=np.array([a+40,b+30,f+150])/S
    opt=least_squares(fun,z0,jac=lambda z:ev(z)['J']/scale[:,None],bounds=(lower,upper),x_scale='jac',ftol=1e-12,xtol=1e-12,gtol=1e-12,max_nfev=100)
    r=ev(opt.x);aa,bb,ff=opt.x*S;status='solved' if max(abs(r['F']/scale))<1e-7 and r['right_residual']<1e-7 else 'unresolved'
    out=dict(case=case,station=i,x_mm=CASES[case]['x_mm'][i-1],alpha_r=float(aa),alpha_i=0.,beta_r=float(bb),f_kHz=float(ff),omega_r=float(ff*W),status=status,scaled_F=(r['F']/scale).tolist(),right_residual=r['right_residual'],scaled_residual=r['scaled_residual'],MAC=r['MAC'],trace=trace,message=opt.message)
    if status=='solved':out['derivative_check']=eq.derivative_check(opt.x)
    file=ROOT/'outputs'/f'{label}_fixed_station_fold.json';vec=file.with_suffix('.npy');np.save(vec,r['q']);out['vector_file']=str(vec);out['vector_sha256']=sha(vec);file.write_text(json.dumps(out,indent=2)+'\n');return out

def one(job):
    case,sign=job;fstart=200 if sign>0 else -112;d=json.loads((SOURCE/'outputs'/f'{case}_f{fstart}kHz_S_final.json').read_text());points=[p for p in d['points'] if p['kind']==0];stations=set(int(p['xi']) for p in points)
    i=max(stations,key=lambda i:np.ptp([p['beta_r'] for p in points if p['xi']==i]));rows=sorted([p for p in points if p['xi']==i],key=lambda p:p['beta_r']);a,b,q=rows[0]['alpha_r'],rows[0]['beta_r'],np.load(rows[0]['vector_file']);f=float(fstart);z=np.array([a,b])/S[:2];step=2. if sign>0 else .2;records=[];fails=[]
    # Neutral continuation stops at a fixed-station fold; this is only a seed.
    for k in range(1000):
        p=Pencil(case,i,f*W);aa,q,r,left=qep(p.spatial(float(z[1]*10)),complex(z[0]*20),q,True);eq=Equations(p,q,left);state=eq.evaluate(np.r_[z,f/100]);J=state['J'][:2,:2];gf=state['J'][:2,2]
        if abs(f)>1500:break
        df=sign*step
        try:
            trial=z+np.linalg.solve(J,-gf*df/100);new=Equations(Pencil(case,i,(f+df)*W),q,left)
            for iteration in range(20):
                rr=new.evaluate(np.r_[trial,(f+df)/100]);F=rr['F'][:2]
                if np.linalg.norm(F)<1e-12 and rr['right_residual']<1e-8:break
                dz=np.linalg.solve(rr['J'][:2,:2],-F)
                if np.linalg.norm(dz)>.15:dz*=.15/np.linalg.norm(dz)
                trial+=dz
            assert np.linalg.norm(F)<1e-12 and rr['right_residual']<1e-8 and rr['MAC']>.97 and np.linalg.norm(trial-z)<.5
        except Exception as error:
            fails.append(dict(f_kHz=f,df=df,reason=repr(error)));step/=2
            if step<.001:break
            continue
        f+=df;z=trial;q=rr['q'];records.append(dict(f_kHz=f,alpha_r=float(z[0]*20),beta_r=float(z[1]*10),MAC=rr['MAC']))
        if iteration<=3:step=min(5. if sign>0 else .5,step*1.2)
    a,b=z*S[:2];label=f'{case}_{"positive" if sign>0 else "negative"}';result=fold(case,i,float(a),float(b),f,q,label);result['continuation']=dict(records=records,failures=fails)
    file=ROOT/'outputs'/f'{label}_pilot.json';file.write_text(json.dumps(result,indent=2)+'\n');print('PILOT',case,sign,result['status'],i,result['f_kHz'],result['x_mm'],result['beta_r'],flush=True);return result
def main():
    snapshot()
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:results=list(pool.map(one,[(case,sign) for case in CASES for sign in [-1,1]]))
    (ROOT/'outputs/pilot_summary.json').write_text(json.dumps(results,indent=2)+'\n')
if __name__=='__main__':main()
