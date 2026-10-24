"""Pseudo-arclength neutral contour of the exact fixed-beta frequency pencil.

Real alpha is solved on Im(omega)=0. Native discrete operators are exact;
between stations their declared linear interpolation defines the contour.
Closure requires a Poincare-section return and an independently corrected
same eigenvector/eigenvalue, never a manually added closing endpoint.
"""
import argparse,datetime,json,math
from pathlib import Path
import numpy as np
from common import ROOT,PREVIOUS,snapshot,sha
from profile_interpolation import evaluate,mac
from fixed_beta_pencil import OperatorLinearFamily

SCALE=np.array([100.,20.])
def stamp():return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()

def main():
    p=argparse.ArgumentParser();p.add_argument('--label',default='M3b7_closed_arc_v1')
    p.add_argument('--max-points',type=int,default=700);p.add_argument('--step',type=float,default=.06)
    p.add_argument('--resume',action='store_true');a=p.parse_args()
    family=OperatorLinearFamily();source=ROOT/'outputs/frequency_pilot_fields_roots.json'
    seed=next(v['roots'][0] for v in json.loads(source.read_text())['results'] if v['station']==119 and v['boundary']=='upper')
    fresh=evaluate('M3',119,seed,seed['omega_r'],1e-6,'M3b7_temporal119_source')
    grid=np.array(json.loads((PREVIOUS/'outputs/seed_catalog.json').read_text())['grids']['main_841']['x_nondim'])*35
    out=ROOT/'outputs'/(a.label+'.json');wave=ROOT/'work'/a.label;wave.mkdir(exist_ok=True)
    def point(z,q,wshift):
        xi,alpha=z*SCALE;pencil=family.pencil(float(xi));r,v=pencil.evaluate(float(alpha),q,wshift)
        i=min(int(np.floor(xi)),840);theta=xi-i
        r.update(xi=float(xi),x_mm=float((1-theta)*grid[i-1]+theta*grid[i]),alpha_i=0.,
            native_bracket=[i,i+1],theta=float(theta),model='native' if min(theta,1-theta)<1e-10 else 'linear_native_operator_interpolation')
        return r,v
    def gradient(r):return np.array([r['omega_xi_i'],r['omega_alpha_i']])*SCALE
    def correct(pred,t,q,wshift,maxit=12):
        z=pred.copy();trace=[];v=q;omega=wshift
        for k in range(maxit):
            if not 1<=z[0]*SCALE[0]<=841:raise RuntimeError('Contour reaches original streamwise domain boundary')
            r,v=point(z,v,omega);omega=r['omega_r']
            F=np.array([r['omega_i'],float(np.dot(z-pred,t))]);g=gradient(r)
            trace.append(dict(point=r,F=F.tolist()))
            if abs(F[0])<1e-10 and abs(F[1])<1e-9:
                assert r['original_real_frequency_polynomial_residual']<1e-7
                return z,r,v,trace
            J=np.array([g,t]);delta=np.linalg.solve(J,-F)
            if np.linalg.norm(delta)>.5:delta*=.5/np.linalg.norm(delta)
            z+=delta
        raise RuntimeError('Pseudo-arclength corrector did not converge')
    z=np.array([119.,seed['alpha_r']])/SCALE;r,q=point(z,fresh['q'],seed['omega_r'])
    assert abs(r['omega_i'])<1e-9 and abs(r['omega_r']-seed['omega_r'])<1e-7
    g=gradient(r);t=np.array([g[1],-g[0]]);t/=np.linalg.norm(t)
    if t[0]<0:t=-t
    startz=z.copy();startq=q.copy();startr=r.copy();startt=t.copy();step=a.step
    state=dict(created=stamp(),beta=-7.,source_path=str(source),source_sha256=sha(source),
        model='Piecewise linear interpolation of original native fixed-beta discrete pencils; no extra CFD grid stations',
        seed_spatial_point=seed,records=[],failures=[],finished=False,status='running',closure=None)
    def store(z,r,q,t,ds,iterations):
        index=len(state['records']);path=wave/('q%04d.npy'%index);np.save(path,q)
        item=dict(index=index,z=z.tolist(),point=r,tangent=t.tolist(),step=ds,iterations=iterations,
            vector_file=str(path),vector_sha256=sha(path));state['records'].append(item)
        tmp=out.with_suffix('.tmp');state['updated']=stamp();tmp.write_text(json.dumps(state,indent=2,allow_nan=False)+'\n');tmp.replace(out)
    if a.resume and out.exists():
        state=json.loads(out.read_text());assert not state['finished'];last=state['records'][-1]
        z=np.array(last['z']);q=np.load(last['vector_file']);r=last['point'];t=np.array(last['tangent']);step=last['step']
    else:
        assert not out.exists(),'Existing phase preserved; use --resume or a new label'
        store(z,r,q,t,step,0)
    accumulated=sum(np.linalg.norm(np.array(b['z'])-np.array(c['z'])) for b,c in zip(state['records'][1:],state['records'][:-1]))
    while len(state['records'])<a.max_points:
        oldz=z.copy();oldq=q.copy();oldr=r.copy();oldsection=float(np.dot(z-startz,startt));accepted=False
        for retry in range(10):
            try:
                pred=z+step*t;zz,rr,v,trace=correct(pred,t,q,r['omega_r'])
                overlap=mac(oldq,v)
                if overlap<.95 or np.linalg.norm(zz-oldz)>step*2:raise RuntimeError('Arc step too long for a single followed vector')
                grad=gradient(rr);tt=np.array([grad[1],-grad[0]]);tt/=np.linalg.norm(tt)
                if np.dot(tt,t)<0:tt=-tt
                if np.dot(tt,t)<.35:raise RuntimeError('Large contour turn; reduce arclength step')
                z,r,q,t=zz,rr,v,tt;accepted=True;accumulated+=np.linalg.norm(z-oldz)
                r['adjacent_raw_MAC']=overlap;store(z,r,q,t,step,len(trace));break
            except Exception as error:
                state['failures'].append(dict(previous_index=len(state['records'])-1,step=step,reason=repr(error)))
                step/=2
                if step<1e-5:break
        if not accepted:state.update(status='unresolved_minimum_arc_step',finished=True);break
        section=float(np.dot(z-startz,startt));distance=float(np.linalg.norm(z-startz))
        if len(state['records'])>30 and accumulated>2 and oldsection<0<=section and distance<max(2*step,.02) and np.dot(t,startt)>.7:
            zz,rr,v,trace=correct(startz,startt,q,r['omega_r'])
            closure=dict(section_return_actual=rr,scaled_distance_to_start=float(np.linalg.norm(zz-startz)),
                raw_MAC_to_start=mac(startq,v),omega_difference=abs(rr['omega_r']-startr['omega_r']),
                point_count_before_canonical_start_duplicate=len(state['records']),independent_correction_trace=trace)
            closure['passed']=closure['scaled_distance_to_start']<1e-7 and closure['raw_MAC_to_start']>1-1e-7 and closure['omega_difference']<1e-7
            state['closure']=closure
            if closure['passed']:
                store(zz,rr,v,startt,0,len(trace));state.update(status='closed_contour_verified',finished=True);break
        if len(state['records'])%10==0:
            print('ARC',len(state['records']),r['xi'],r['x_mm'],r['omega_r'],r['alpha_r'],step,flush=True)
        if len(trace)<=3 and overlap>.99:step=min(.15,step*1.15)
    if not state['finished']:state.update(status='point_budget_reached_not_closed',finished=True)
    state['updated']=stamp();state['original_inputs']=snapshot();out.write_text(json.dumps(state,indent=2,allow_nan=False)+'\n')
    print('ARC_FINAL',state['status'],len(state['records']),flush=True)

if __name__=='__main__':main()
