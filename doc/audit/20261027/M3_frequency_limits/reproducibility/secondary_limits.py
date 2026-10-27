"""Continue independently discovered downstream unstable pockets to extinction."""
import concurrent.futures,json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/neutral_m3_frequency_limits_20261027')
sys.path.insert(0,str(ROOT/'work'))
from common import CASES
from solver import Pencil,original,mac
from pilot import fold
from fold_core import W

def peak(case,i,f,b,a,q):
    p=Pencil(case,i,f*W);trace=[]
    for it in range(30):
        a,q,r=p.solve_alpha(float(b),complex(a),q,True);h=.02;am,_,_=p.solve_alpha(b-h,a,q);ap,_,_=p.solve_alpha(b+h,a,q);curvature=(ap.imag-2*a.imag+am.imag)/(h*h);der=r['alpha_beta_i'];trace.append(dict(beta=b,alpha_r=a.real,alpha_i=a.imag,derivative=der,curvature=curvature))
        if abs(der)<1e-9:break
        if curvature>0:delta=float(np.clip(-der/curvature,-.8,.8))
        else:delta=-np.sign(der)*.2
        accepted=False
        for retry in range(15):
            aa,v,rr=p.solve_alpha(b+delta,a,q,True)
            if aa.imag<a.imag+1e-12:b+=delta;a,q=aa,v;accepted=True;break
            delta/=2
        assert accepted,('Peak stalled',case,i,f,trace[-1])
    assert abs(der)<1e-6 and curvature>0,('Unresolved peak',case,i,f,trace[-1])
    return dict(f_kHz=float(f),beta_r=float(b),alpha_r=a.real,alpha_i=a.imag,growth=-a.imag,trace=trace),q

def one(case):
    i=360 if case=='bluntplate' else 90;source=ROOT/'outputs'/f'{case}_outside_negative_i{i:04d}.json';d=json.loads(source.read_text());m=d['maximum_growth'];f=d['target_f_kHz'];b=m['beta'];a=complex(m['alpha_r'],m['alpha_i']);state=dict(alpha_r=a.real,alpha_i=a.imag,beta_r=b,omega_r=f*W);actual,q,receipt=original(case,i,state,b,case+'_secondary_unstable_independent');assert abs(actual['growth']-m['growth'])<1e-7
    row,q=peak(case,i,f,b,a,q);records=[row];cache={float(f):(row,q)};step=2.;errors=[]
    for it in range(500):
        newf=f-step
        try:new,v=peak(case,i,newf,row['beta_r'],complex(row['alpha_r'],row['alpha_i']),q);assert mac(q,v)>.97
        except Exception as e:
            step/=2;errors.append(dict(f_kHz=f,step=step,reason=repr(e)))
            if step<1e-4:raise RuntimeError(('Secondary frequency unresolved',case,f,errors[-1]))
            continue
        records.append(new);cache[newf]=(new,v)
        if new['growth']<0:break
        f,row,q=newf,new,v
        if len(new['trace'])<6:step=min(5.,step*1.2)
    else:raise RuntimeError(('Frequency range exhausted',case,f))
    def value(ff):
        if ff not in cache:
            key=min(cache,key=lambda k:abs(ff-k));r,q=cache[key];cache[ff]=peak(case,i,ff,r['beta_r'],complex(r['alpha_r'],r['alpha_i']),q)
        return cache[ff][0]['growth']
    fc=brentq(value,newf,f,xtol=1e-8);row,q=cache[fc];result=fold(case,i,row['alpha_r'],row['beta_r'],fc,q,case+'_negative_secondary');assert result['status']=='solved';result.update(source=str(source),initial_independent_check=receipt,frequency_peak_path=records,peak_path_failures=errors)
    (ROOT/'outputs'/f'{case}_negative_secondary_pilot.json').write_text(json.dumps(result,indent=2)+'\n');print('SECONDARY_FIXED_FOLD',case,i,fc,result['beta_r'],result['alpha_r'],flush=True);return result
def main():
    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as pool:results=list(pool.map(one,CASES))
    (ROOT/'outputs/secondary_pilot_summary.json').write_text(json.dumps(results,indent=2)+'\n')
if __name__=='__main__':main()
