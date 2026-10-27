"""Optimize actual-profile frequency folds around all leading native extrema."""
import concurrent.futures,json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/neutral_m3_frequency_limits_20261027')
sys.path.insert(0,str(ROOT/'work'))
from common import CASES,sha
from solver import Pencil,qep,mac,original
from physical import Pair
from trace_folds import solve
from fold_core import S,W

def one(job):
    case,sign=job[:2];variant=job[2] if len(job)>2 else '';label=f'{case}_{"positive" if sign>0 else "negative"}'+('_'+variant if variant else '');source=ROOT/'outputs'/f'{label}_native_fold_curve.json';d=json.loads(source.read_text());records={p['station']:p for p in d['records']};peak=(max if sign>0 else min)(records.values(),key=lambda p:p['f_kHz']);i0=peak['station'];pairs={};cache={};trace=[]
    def evaluate(x):
        key=float(x).hex()
        if key in cache:return cache[key]
        i=min(int(np.floor(x)),CASES[case]['stations']-1);t=x-i;assert i in records and i+1 in records
        if i not in pairs:pairs[i]=Pair(case,i,0.,label+'_frequency_limit')
        p=pairs[i].pencil(float(t));a,b=records[i],records[i+1];z=np.array([(1-t)*a[k]+t*b[k] for k in ['alpha_r','beta_r','f_kHz']])/S;q=np.load((a if t<.5 else b)['vector_file']);row,v,zz=solve(p,z,q,max_nfev=100)
        row.update(xi=float(x),native_bracket=[i,i+1],theta=float(t),x_mm=(1-t)*CASES[case]['x_mm'][i-1]+t*CASES[case]['x_mm'][i],generation=str(pairs[i].build(float(t))/'generation.json'))
        cache[key]=(row,v,p);trace.append({k:row[k] for k in ['xi','x_mm','f_kHz','alpha_r','beta_r','right_residual','raw_MAC']});return cache[key]
    choices=[]
    for i in range(max(min(records),i0-2),min(max(records),i0+3)):
        if i not in records or i+1 not in records:continue
        for x in [float(i),float(i+1)]:choices.append(evaluate(x))
        opt=minimize_scalar(lambda t:-sign*evaluate(i+float(t))[0]['f_kHz'],bounds=(1e-6,1-1e-6),method='bounded',options=dict(xatol=2e-6,maxiter=40));choices.append(evaluate(i+float(opt.x)))
    row,q,p=(max if sign>0 else min)(choices,key=lambda item:item[0]['f_kHz']);qfile=ROOT/'outputs'/f'{label}_critical.npy';np.save(qfile,q)
    f,a,b=row['f_kHz'],row['alpha_r'],row['beta_r'];p.omega=f*W;p.C=[list(c) for c in p.baseC];p.C[0][0]=p.C[0][0]+p.omega*p.G
    aa,v,r=p.solve_alpha(b,a,q,True);assert abs(aa.imag)<1e-8 and abs(r['alpha_beta_i'])<1e-7 and mac(q,v)>1-1e-7
    directory=Path(row['generation']).parent;state=dict(alpha_r=a,alpha_i=0.,beta_r=b,omega_r=f*W);actual,vector,receipt=original(directory,5,state,b,label+'_critical_independent');assert abs(actual['growth'])<1e-8 and mac(vector,q)>1-1e-7
    tests=[]
    for delta in [.01,.1,1.]:
        for side in [-1,1]:
            frequency=f+side*sign*delta;p.omega=frequency*W;p.C=[list(c) for c in p.baseC];p.C[0][0]=p.C[0][0]+p.omega*p.G;evals=[]
            def growth(beta):
                alpha,v,r=p.solve_alpha(float(beta),complex(a),q,True);assert mac(q,v)>.95;evals.append(dict(beta=float(beta),alpha_r=alpha.real,alpha_i=alpha.imag,growth=-alpha.imag,MAC=mac(q,v)));return alpha.imag
            opt=minimize_scalar(growth,bounds=(b-2.,b+2.),method='bounded',options=dict(xatol=1e-8));alpha_i=growth(float(opt.x));assert (alpha_i<0 if side<0 else alpha_i>0),(label,delta,side,alpha_i)
            tests.append(dict(f_kHz=frequency,side='inside' if side<0 else 'outside',distance_kHz=delta,beta_at_local_max_growth=float(opt.x),maximum_growth=float(-alpha_i),evaluations=evals))
    result=dict(case=case,sign=sign,status='actual_profile_extremum_verified',critical=row,vector_file=str(qfile),vector_sha256=sha(qfile),native_extreme=peak,native_source=str(source),native_source_sha256=sha(source),native_source_status=d['status'],physical_samples=trace,alpha_beta_i=r['alpha_beta_i'],original_check=receipt,frequency_side_checks=tests,scope='Actual-profile refinement around the largest native extremum of this tracked neutral-merger branch; full branch coverage and outside scans audited separately')
    out=ROOT/'outputs'/f'{label}_refined_limit.json';out.write_text(json.dumps(result,indent=2)+'\n');print('REFINED_LIMIT',label,row['f_kHz'],row['x_mm'],row['beta_r'],row['alpha_r'],flush=True);return result
def main():
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:results=list(pool.map(one,[(case,sign) for case in CASES for sign in [-1,1]]))
    (ROOT/'outputs/refined_limit_summary.json').write_text(json.dumps(results,indent=2)+'\n')
if __name__=='__main__':main()
