"""Continue beta-neutral mergers through folds of their station projection.

Native-operator guide only. Actual profiles and the frozen solver independently
validate any resulting frequency extremum before it becomes a reported limit.
"""
import concurrent.futures,json,sys
from pathlib import Path
import numpy as np
from scipy.linalg import null_space
ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/neutral_m3_frequency_limits_20261027')
sys.path.insert(0,str(ROOT/'work'))
from common import CASES,sha
from solver import Pencil,qep,mac
from contour import Family
from trace_folds import blended
from fold_core import Equations,S,W
XS=50.
S4=np.r_[XS,S]
def one(job):
    case,variant=job if isinstance(job,tuple) else (job,'secondary');source=ROOT/'outputs'/(f'{case}_negative_secondary_native_fold_curve.json' if variant=='secondary' else f'{case}_negative_tail_fold_curve.json');d=json.loads(source.read_text());seed=min(d['records'],key=lambda r:r['station']);family=Family(case,0.);known={}
    for kind in ['negative_native','negative_tail','negative_secondary_native']:
        file=ROOT/'outputs'/f'{case}_{kind}_fold_curve.json'
        if file.exists():
            for row in json.loads(file.read_text())['records']:known.setdefault(row['station'],[]).append((kind,row))
    p=Pencil(case,seed['station'],seed['omega_r']);_,q0,_,p0=qep(p.spatial(seed['beta_r']),seed['alpha_r'],np.load(seed['vector_file']),True);z0=np.array([seed['station'],seed['alpha_r'],seed['beta_r'],seed['f_kHz']])/S4
    def pencil(x,side):
        C,dC,i,t=family.pencil(float(x),side);a,b=family.native(i),family.native(i+1);return blended(a,b,float(t),0.),i,float(t)
    def evaluate(z,qref,pref,side=1):
        x,a,b,f=z*S4;p,i,t=pencil(x,side);eq=Equations(p,qref,pref);r=eq.evaluate(z[1:]);h=1e-4
        # One-sided station derivative at a native corner, otherwise central.
        left=max(float(i),x-h);right=min(float(i+1),x+h);assert right-left>1e-10
        lp,_,_=pencil(left,1);rp,_,_=pencil(right,-1);fl=Equations(lp,qref,pref).evaluate(z[1:])['F'];fr=Equations(rp,qref,pref).evaluate(z[1:])['F'];J=np.column_stack([(fr-fl)/(right-left)*XS,r['J']]);sc=eq.scale(z[1:]);U,s,V=np.linalg.svd(J/sc[:,None],full_matrices=True);assert s[-1]>s[0]*1e-11,('Critical curve rank loss',s,z*S4);tangent=V[-1]
        point=dict(xi=float(x),x_mm=float((1-t)*family.x[i-1]+t*family.x[i]),alpha_r=float(a),beta_r=float(b),f_kHz=float(f),omega_r=float(W*f),right_residual=r['right_residual'],scaled_F=(r['F']/sc).tolist(),MAC=r['MAC'])
        return r['F'],J,sc,tangent,r['q'],r['p'],point
    _,_,_,t,q,p,row=evaluate(z0,q0,p0,-1)
    if t[0]>0:t=-t
    path=ROOT/'outputs'/f'{case}_negative_{variant}_arc.json';wave=ROOT/'work'/f'{case}_negative_{variant}_arc_vectors';wave.mkdir(exist_ok=True);state=dict(case=case,variant=variant,status='running',records=[],failures=[],source=str(source),source_sha256=sha(source),model='Piecewise native-operator guide; no physical off-grid claims')
    def store(z,t,q,p,row):
        k=len(state['records']);qfile=wave/f'q{k:04d}.npy';pfile=wave/f'p{k:04d}.npy';np.save(qfile,q);np.save(pfile,p);state['records'].append(dict(z=z.tolist(),tangent=t.tolist(),point=row,vector_file=str(qfile),left_vector_file=str(pfile)))
        if k%10==0:path.write_text(json.dumps(state,indent=2)+'\n')
    z=z0.copy();store(z,t,q,p,row);step=.005;arc=0.
    def correct(pred,plane,q,p):
        zz=pred.copy();px=pred[0]*XS;cell=min(int(np.floor(px)),family.count-1)
        if abs(px-round(px))<1e-8 and plane[0]<0:cell=max(1,int(round(px))-1)
        lo,hi=cell/XS,(cell+1)/XS;constrained=np.linalg.norm(plane[1:])>1e-12
        for iteration in range(35):
            side=-1 if constrained and abs(zz[0]*XS-(cell+1))<1e-8 else 1
            F,J,scale,tt,v,w,rr=evaluate(zz,q,p,side);value=np.r_[F/scale,np.dot(zz-pred,plane)];A=np.vstack([J/scale[:,None],plane]);weights=1/np.maximum(np.linalg.norm(A,axis=1),1e-30)
            if max(abs(value))<1e-8 and rr['right_residual']<1e-7:return zz,tt,v,w,rr
            delta=np.linalg.solve(A*weights[:,None],-value*weights)
            if np.linalg.norm(delta)>.08:delta*=.08/np.linalg.norm(delta)
            factor=1.
            if constrained and delta[0]>0 and zz[0]+delta[0]>hi:factor=.8*(hi-zz[0])/delta[0]
            if constrained and delta[0]<0 and zz[0]+delta[0]<lo:factor=.8*(lo-zz[0])/delta[0]
            assert factor>1e-9;zz+=factor*delta
        raise RuntimeError(('Arc corrector limit',rr))
    for loop in range(2500):
        if abs(z[0]*XS-round(z[0]*XS))<1e-8:
            _,_,_,tt,q,p,row=evaluate(z,q,p,int(np.sign(t[0])))
            if np.dot(tt,t)<0:tt=-tt
            t=tt
        accepted=False
        for retry in range(16):
            try:
                pred=z+step*t;x=z[0]*XS;boundary=int(np.floor(x+1e-8))+1 if t[0]>0 else int(np.ceil(x-1e-8))-1
                if not 1<=boundary<=family.count:state['status']='domain_boundary';break
                if (pred[0]*XS-boundary)*t[0]>=0:
                    pred=z+(boundary/XS-z[0])/t[0]*t;plane=np.array([1.,0.,0.,0.])
                else:plane=t
                zz,tt,v,w,rr=correct(pred,plane,q,p)
                if np.dot(tt,t)<0:tt=-tt
                assert np.dot(tt,t)>.3 and mac(q,v)>.95 and np.linalg.norm(zz-z)<2*step
                arc+=np.linalg.norm(zz-z);z,t,q,p,row=zz,tt,v,w,rr;store(z,t,q,p,row);accepted=True;break
            except Exception as e:
                state['failures'].append(dict(index=len(state['records'])-1,step=step,reason=repr(e)));step/=2
                if step<1e-6:break
        if state['status']=='domain_boundary':break
        if not accepted:state['status']='unresolved_step';break
        if len(state['records'])>30 and abs(row['xi']-round(row['xi']))<1e-8:
            i=round(row['xi'])
            for kind,knownrow in known.get(i,[]):
                error=np.linalg.norm((np.array([row[k]-knownrow[k] for k in ['alpha_r','beta_r','f_kHz']]))/S)
                if error<1e-6 and arc>.2 and mac(q,np.load(knownrow['vector_file']))>1-1e-6:state.update(status='connected_to_known_branch',connection=dict(kind=kind,station=i,error=float(error)));break
            if state['status']!='running':break
        step=min(.04,step*1.1)
        if len(state['records'])%50==0:print('ARC_FOLD',case,len(state['records']),row['x_mm'],row['f_kHz'],row['alpha_r'],flush=True)
    else:state['status']='point_limit'
    state['minimum_frequency']=min(state['records'],key=lambda r:r['point']['f_kHz']);path.write_text(json.dumps(state,indent=2)+'\n');print('ARC_FOLD_FINAL',case,state['status'],len(state['records']),state['minimum_frequency']['point'],flush=True);return state
def main():
    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as pool:results=list(pool.map(one,CASES))
    (ROOT/'outputs/arc_fold_summary.json').write_text(json.dumps(results,indent=2)+'\n')
if __name__=='__main__':main()
