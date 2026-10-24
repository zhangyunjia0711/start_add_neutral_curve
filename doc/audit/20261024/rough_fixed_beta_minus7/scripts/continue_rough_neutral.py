"""Fixed-beta pseudo-arclength contour from the physically aligned rough M3 seed."""
import argparse,datetime,json
import numpy as np
from common import ROOT,sha,snapshot
from profile_interpolation import evaluate,mac
from fixed_beta_pencil import OperatorLinearFamily

SCALE=np.array([50.,20.])
def stamp():return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--label',default='M3b7_closed_arc_v2');parser.add_argument('--max-points',type=int,default=1200)
    parser.add_argument('--resume',action='store_true');parser.add_argument('--resume-from');a=parser.parse_args()
    source=ROOT/'outputs/physical_M3_frequency_positive.json';survey=json.loads(source.read_text());assert survey['finished']
    seed=survey['roots'][0];station=survey['station'];s=dict(seed['point']);s['alpha_i']=0.
    fresh=evaluate('M3',station,s,s['omega_r'],1e-7,'rough_M3_neutral_start_spatial')
    assert abs(fresh['point']['growth'])<1e-8 and mac(fresh['q'],np.load(seed['vector_file']))>1-1e-7
    family=OperatorLinearFamily();grid=np.array(json.loads((ROOT/'reference/outputs/seed_catalog.json').read_text())['grids']['main_841']['x_nondim'])*35
    out=ROOT/'outputs'/(a.label+'.json');wave=ROOT/'work'/a.label;wave.mkdir(exist_ok=True)
    def point(z,q,w):
        xi,alpha=z*SCALE
        if not 1<=xi<=361:raise RuntimeError('Contour reaches the original streamwise domain boundary')
        p=family.pencil(float(xi));r,v=p.evaluate(float(alpha),q,w);i=p.native_bracket[0]
        r.update(xi=float(xi),x_mm=float((1-p.theta)*grid[i-1]+p.theta*grid[i]),alpha_i=0.,native_bracket=p.native_bracket,theta=p.theta,
            model='native' if min(p.theta,1-p.theta)<1e-10 else 'linear_native_operator_interpolation')
        return r,v
    def gradient(r):return np.array([r['omega_xi_i'],r['omega_alpha_i']])*SCALE
    def correct(pred,t,q,w,maxit=12):
        z=pred.copy();trace=[];v=q
        for k in range(maxit):
            r,v=point(z,v,w);w=r['omega_r'];F=np.array([r['omega_i'],np.dot(z-pred,t)]);trace.append(dict(point=r,F=F.tolist()))
            if abs(F[0])<1e-10 and abs(F[1])<1e-9:
                assert r['original_real_frequency_polynomial_residual']<1e-7
                return z,r,v,trace
            delta=np.linalg.solve(np.array([gradient(r),t]),-F)
            if np.linalg.norm(delta)>.3:delta*=.3/np.linalg.norm(delta)
            z+=delta
        raise RuntimeError('Neutral arclength corrector limit')
    z=np.array([station,s['alpha_r']])/SCALE;r,q=point(z,fresh['q'],s['omega_r'])
    assert abs(r['omega_i'])<1e-8 and abs(r['omega_r']-s['omega_r'])<1e-7
    g=gradient(r);t=np.array([g[1],-g[0]]);t/=np.linalg.norm(t)
    if t[0]<0:t=-t
    startz=z.copy();startq=q.copy();startr=r.copy();startt=t.copy();step=.03
    state=dict(created=stamp(),beta=-7.,source_path=str(source),source_sha256=sha(source),scale=SCALE.tolist(),
        model='Guide only: piecewise linear native-operator interpolation; final points will be solved on native stations or actual interpolated local profiles.',
        seed_spatial_point=fresh['point'],seed_spatial_run=fresh['run'],records=[],failures=[],finished=False,status='running')
    def store(z,r,q,t,ds,iterations):
        j=len(state['records']);path=wave/('q%04d.npy'%j);np.save(path,q)
        state['records'].append(dict(index=j,z=z.tolist(),point=r,tangent=t.tolist(),step=ds,iterations=iterations,vector_file=str(path),vector_sha256=sha(path)))
        state['updated']=stamp();tmp=out.with_suffix('.tmp');tmp.write_text(json.dumps(state,indent=2,allow_nan=False)+'\n');tmp.replace(out)
    if a.resume_from:
        previous=ROOT/'outputs'/a.resume_from;state=json.loads(previous.read_text());assert not out.exists()
        state.update(finished=False,status='running',resumed_from=str(previous),resumed_from_sha256=sha(previous));last=state['records'][-1]
        z=np.array(last['z']);q=np.load(last['vector_file']);r=last['point'];t=np.array(last['tangent']);step=max(last['step'],.003)
    elif a.resume and out.exists():
        state=json.loads(out.read_text());assert not state['finished'];last=state['records'][-1]
        z=np.array(last['z']);q=np.load(last['vector_file']);r=last['point'];t=np.array(last['tangent']);step=last['step']
    else:
        assert not out.exists();store(z,r,q,t,step,0)
    accumulated=sum(np.linalg.norm(np.array(b['z'])-np.array(c['z'])) for b,c in zip(state['records'][1:],state['records'][:-1]))
    while len(state['records'])<a.max_points:
        oldz=z.copy();oldq=q.copy();oldsection=float(np.dot(z-startz,startt));accepted=False
        for retry in range(12):
            try:
                pred=z+step*t;xi=float(z[0]*SCALE[0]);px=float(pred[0]*SCALE[0])
                boundary=int(np.floor(xi+1e-8))+1 if t[0]>0 else int(np.ceil(xi-1e-8))-1
                native_bridge=1<=boundary<=361 and (px-boundary)*t[0]>=0
                if native_bridge:
                    plane=np.array([1.,0.]);pred=z+(boundary/SCALE[0]-z[0])/t[0]*t
                    zz,rr,v,trace=correct(pred,plane,q,r['omega_r'],maxit=30)
                    assert rr['omega_alpha_i']*r['omega_alpha_i']>0,'Native bridge would skip a real fold'
                else:zz,rr,v,trace=correct(pred,t,q,r['omega_r'])
                overlap=mac(oldq,v)
                if overlap<.95 or np.linalg.norm(zz-oldz)>step*2:raise RuntimeError('Arc-step continuity gate')
                grad=gradient(rr);tt=np.array([grad[1],-grad[0]]);tt/=np.linalg.norm(tt)
                if native_bridge:
                    if tt[0]*t[0]<0:tt=-tt
                    rr['native_grid_corner_bridge']=True
                else:
                    if np.dot(tt,t)<0:tt=-tt
                    if np.dot(tt,t)<.35:raise RuntimeError('Sharp turn; reduce arclength step')
                z,r,q,t=zz,rr,v,tt;accepted=True;accumulated+=np.linalg.norm(z-oldz);r['adjacent_raw_MAC']=overlap
                store(z,r,q,t,step,len(trace));break
            except Exception as error:
                state['failures'].append(dict(previous_index=len(state['records'])-1,step=step,reason=repr(error)));step/=2
                if step<1e-5:break
        if not accepted:state.update(status='unresolved_minimum_arc_step',finished=True);break
        section=float(np.dot(z-startz,startt));distance=float(np.linalg.norm(z-startz))
        if len(state['records'])>30 and accumulated>2 and oldsection<0<=section and distance<max(2*step,.02) and np.dot(t,startt)>.7:
            pred=oldz+(-oldsection)/(section-oldsection)*(z-oldz)
            zz,rr,v,trace=correct(pred,startt,q,r['omega_r'])
            closure=dict(section_interpolation_predictor=pred.tolist(),initial_distance=float(np.linalg.norm(pred-startz)),
                section_return_actual=rr,scaled_distance_to_start=float(np.linalg.norm(zz-startz)),raw_MAC_to_start=mac(startq,v),
                omega_difference=abs(rr['omega_r']-startr['omega_r']),independent_correction_trace=trace)
            closure['passed']=closure['scaled_distance_to_start']<1e-7 and closure['raw_MAC_to_start']>1-1e-7 and closure['omega_difference']<1e-7
            state['closure']=closure
            if closure['passed']:
                store(zz,rr,v,startt,0,len(trace));state.update(status='closed_contour_verified',finished=True);break
        if len(state['records'])%10==0:print('ARC',len(state['records']),r['xi'],r['x_mm'],r['omega_r'],r['alpha_r'],step,flush=True)
        if len(trace)<=3 and overlap>.99:step=min(.1,step*1.15)
    if not state['finished']:state.update(status='point_budget_reached_not_closed',finished=True)
    state['updated']=stamp();state['original_inputs']=snapshot();out.write_text(json.dumps(state,indent=2,allow_nan=False)+'\n')
    print('ARC_FINAL',state['status'],len(state['records']),flush=True)

if __name__=='__main__':main()
