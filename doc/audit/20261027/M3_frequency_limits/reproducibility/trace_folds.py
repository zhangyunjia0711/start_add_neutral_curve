"""Track the two neutral-beta merger curves along every original CFD station."""
import concurrent.futures,json,sys,time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/neutral_m3_frequency_limits_20261027')
sys.path.insert(0,str(ROOT/'work'))
from common import CASES,snapshot,sha
from solver import Pencil,qep,mac
from fold_core import Equations,S,W

def solve(p,zguess,qguess,max_nfev=60):
    a,b,f=zguess*S;p.omega=f*W;p.C=[list(row) for row in p.baseC];p.C[0][0]=p.C[0][0]+p.omega*p.G
    aa,q0,record,p0=qep(p.spatial(float(b)),complex(a),qguess,True)
    assert mac(qguess,q0)>.95
    eq=Equations(p,q0,p0);z0=np.array([aa.real,b,f])/S;scale=eq.scale(z0);cache={};trace=[]
    def ev(z):
        key=tuple(z)
        if key not in cache:cache.clear();cache[key]=eq.evaluate(z)
        return cache[key]
    def fun(z):
        r=ev(z);assert r['MAC']>.95
        trace.append(dict(z=z.tolist(),scaled_F=(r['F']/scale).tolist(),residual=r['right_residual'],MAC=r['MAC']))
        return r['F']/scale
    lower=z0-np.array([1.5,2.,1.]);upper=z0+np.array([1.5,2.,1.])
    opt=least_squares(fun,z0,jac=lambda z:ev(z)['J']/scale[:,None],bounds=(lower,upper),x_scale='jac',ftol=2e-12,xtol=2e-12,gtol=2e-12,max_nfev=max_nfev)
    r=ev(opt.x);assert max(abs(r['F']/scale))<1e-7 and r['right_residual']<1e-7,('Fold not solved',opt.message,trace[-1])
    overlap=mac(qguess,r['q']);assert overlap>.95,('Fold continuity',overlap)
    a,b,f=opt.x*S
    return dict(alpha_r=float(a),alpha_i=0.,beta_r=float(b),f_kHz=float(f),omega_r=float(f*W),right_residual=r['right_residual'],coefficient_scaled_residual=r['scaled_residual'],scaled_F=(r['F']/scale).tolist(),raw_MAC=overlap,iterations=len(trace)),r['q'],opt.x

def blended(a,b,t,f):
    p=object.__new__(Pencil);p.baseC=[[(1-t)*a.baseC[k][j]+t*b.baseC[k][j] for j in range(3)] for k in range(3)];p.G=(1-t)*a.G+t*b.G;p.C=[list(row) for row in p.baseC];p.omega=f*W;p.C[0][0]=p.C[0][0]+p.omega*p.G;return p

def one(job):
    case,sign=job[:2];variant=job[2] if len(job)>2 else '';label=f'{case}_{"positive" if sign>0 else "negative"}'+('_'+variant if variant else '');source=ROOT/'outputs'/f'{label}_pilot.json';pilot=json.loads(source.read_text());assert pilot['status']=='solved';seed=pilot['station'];q0=np.load(pilot['vector_file']);z0=np.array([pilot['alpha_r'],pilot['beta_r'],pilot['f_kHz']])/S
    path=ROOT/'outputs'/f'{label}_native_fold_curve.json';wave=ROOT/'work'/f'{label}_fold_vectors';wave.mkdir(exist_ok=True)
    state=dict(case=case,sign=sign,source=str(source),source_sha256=sha(source),records=[],failures=[],directions=[],status='running')
    seed_record=None
    def store(i,row,q):
        qp=wave/f'i{i:04d}.npy';np.save(qp,q);row=dict(row,station=int(i),xi=float(i),x_mm=CASES[case]['x_mm'][i-1],vector_file=str(qp),vector_sha256=sha(qp));state['records'].append(row)
        if len(state['records'])%10==0:path.write_text(json.dumps(state,indent=2)+'\n')
        return row
    p=Pencil(case,seed,z0[2]*100*W);row,q,z=solve(p,z0,q0);seed_record=store(seed,row,q);z0=z.copy();q0=q.copy()
    for direction,end in [(-1,1),(1,CASES[case]['stations'])]:
        z=z0.copy();q=q0.copy();previous=seed;completed=True
        for i in range(seed+direction,end+direction,direction):
            p=Pencil(case,i,z[2]*100*W)
            try:row,v,zz=solve(p,z,q)
            except Exception as e:
                state['failures'].append(dict(from_station=previous,to_station=i,reason=repr(e)));ok=False
                a=Pencil(case,previous,z[2]*100*W);b=p
                for divisions in [4,16]:
                    try:
                        localz=z.copy();localq=q.copy()
                        for theta in np.linspace(0,1,divisions+1)[1:]:row,localq,localz=solve(blended(a,b,float(theta),float(localz[2]*100)),localz,localq)
                        row,v,zz=solve(p,localz,localq);row['operator_substeps']=divisions;ok=True;break
                    except Exception as ee:state['failures'].append(dict(from_station=previous,to_station=i,subdivisions=divisions,reason=repr(ee)))
                if not ok:completed=False;break
            store(i,row,v);q,z=v,zz;previous=i
            if len(state['records'])%100==0:print('FOLD_TRACE',label,len(state['records']),i,row['x_mm'] if 'x_mm' in row else CASES[case]['x_mm'][i-1],row['f_kHz'],flush=True)
        state['directions'].append(dict(direction=direction,last_station=previous,complete=completed))
    state['records'].sort(key=lambda r:r['station']);state['status']='complete_all_native_stations' if all(r['complete'] for r in state['directions']) else 'partial_with_recorded_failures'
    extreme=(max if sign>0 else min)(state['records'],key=lambda r:r['f_kHz']);state['native_extreme']=extreme;state['station_count']=len(state['records']);path.write_text(json.dumps(state,indent=2)+'\n');print('FOLD_CURVE',label,state['status'],len(state['records']),extreme['station'],extreme['x_mm'],extreme['f_kHz'],extreme['beta_r'],flush=True);return dict(case=case,sign=sign,status=state['status'],file=str(path),native_extreme=extreme)
def main():
    snapshot()
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:results=list(pool.map(one,[(case,sign) for case in CASES for sign in [-1,1]]))
    (ROOT/'outputs/native_fold_summary.json').write_text(json.dumps(results,indent=2)+'\n');snapshot()
if __name__=='__main__':main()
