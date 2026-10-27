"""Stationary x-beta contour: quadratic beta eigenproblem, real alpha.

Native operator interpolation guides topology. Delivery is corrected on native
CFD stations; guide points between them are explicitly marked as guides.
"""
import argparse,json,datetime,time
from pathlib import Path
from functools import lru_cache
import numpy as np
from common import ROOT,CASES,sha,snapshot
from solver import Pencil,mac,qep,original
SCALE=np.array([50.,20.])

class Family:
    def __init__(self,case,omega):self.case=case;self.omega=omega;self.count=CASES[case]['stations'];self.x=np.array(CASES[case]['x_mm']);self.native=lru_cache(maxsize=6)(self._native)
    def _native(self,i):
        path=ROOT/'native_pencils'/self.case/f'i{i:03d}'/'frequency_coefficient.json';start=time.monotonic()
        while not path.exists():
            assert time.monotonic()-start<900,('Waiting for matrix preparation',self.case,i)
            time.sleep(.2)
        return Pencil(self.case,i,self.omega)
    def pencil(self,xi,side=1):
        if abs(xi-round(xi))<1e-9:xi=float(round(xi))
        assert 1<=xi<=self.count,('Domain boundary',xi)
        i=min(int(np.floor(xi)),self.count-1)
        if side<0 and abs(xi-round(xi))<1e-9:i=max(1,int(round(xi))-1)
        theta=xi-i;a,b=self.native(i),self.native(i+1)
        C=[[(1-theta)*a.C[k][j]+theta*b.C[k][j] for j in range(3)] for k in range(3)]
        dC=[[b.C[k][j]-a.C[k][j] for j in range(3)] for k in range(3)]
        return C,dC,i,theta
    def evaluate(self,z,q,target,side=1):
        xi,alpha=z*SCALE;C,dC,i,theta=self.pencil(float(xi),side)
        K=[sum(alpha**k*C[k][j] for k in range(3)) for j in range(3)]
        beta,v,r,p=qep(K,target,q,True)
        Pb=K[1]+2*beta*K[2]
        Pa=sum(beta**j*(C[1][j]+2*alpha*C[2][j]) for j in range(3))
        Px=sum(alpha**k*beta**j*dC[k][j] for k in range(3) for j in range(3))
        den=np.vdot(p,Pb@v);ba=-np.vdot(p,Pa@v)/den;bx=-np.vdot(p,Px@v)/den
        r.update(xi=float(xi),x_mm=float((1-theta)*self.x[i-1]+theta*self.x[i]),alpha_r=float(alpha),alpha_i=0.,omega_r=float(self.omega),omega_i=0.,
            beta_r=float(beta.real),beta_i=float(beta.imag),beta_alpha_i=float(ba.imag),beta_xi_i=float(bx.imag),
            beta_alpha_r=float(ba.real),beta_xi_r=float(bx.real),native_bracket=[i,i+1],theta=float(theta),
            real_beta_polynomial_residual=float(np.linalg.norm(sum(alpha**k*beta.real**j*C[k][j] for k in range(3) for j in range(3))@v)))
        return r,v

def seed_root(case,frequency,component):
    source=ROOT/'outputs/neutral_seeds.json';data=json.loads(source.read_text())
    s=next(r for r in data['records'] if r['case']==case and r['f_kHz']==frequency and r['side']=='lower')
    q=np.load(s['vector_file']);i=s['station']
    if s.get('kind')=='operator_guide':
        return dict(station=i,point=s,receipt=dict(model='Original-operator guide seed; independent physical checks performed after finalization'),component=component,frequency_kHz=frequency,source_record=s),q
    fresh,v,receipt=original(case,i,s,s['beta_r'],f'{case}_f{frequency}kHz_S_neutral_seed')
    assert abs(fresh['growth'])<1e-8 and mac(v,q)>1-1e-7
    return dict(station=i,point=fresh,receipt=receipt,component=component,frequency_kHz=frequency,source_record=s),v


def main():
    parser=argparse.ArgumentParser();parser.add_argument('case');parser.add_argument('frequency',type=int,choices=[50,100]);parser.add_argument('component',choices=['A','B']);parser.add_argument('--resume',action='store_true');a=parser.parse_args();case=a.case;conversion=json.loads((ROOT/'outputs/frequency_conversion.json').read_text());omega=next(r['omega_nd'] for r in conversion['frequencies'] if r['f_kHz']==a.frequency);family=Family(case,omega)
    out=ROOT/'outputs'/f'{case}_f{a.frequency}kHz_{a.component}_guide.json';wave=ROOT/'work'/f'{case}_f{a.frequency}kHz_{a.component}_guide';wave.mkdir(exist_ok=True)
    seed,q=seed_root(case,a.frequency,a.component);s=seed['point'];z=np.array([seed['station'],s['alpha_r']])/SCALE
    r,q=family.evaluate(z,q,s['beta_r'])
    def gradient(r):return np.array([r['beta_xi_i'],r['beta_alpha_i']])*SCALE
    def correct(pred,t,q,b,maxit=14):
        zz=pred.copy();trace=[];v=q
        constrained=abs(t[1])>1e-12
        px=float(pred[0]*SCALE[0]);cell=min(int(np.floor(px)),family.count-1)
        if abs(px-round(px))<1e-9 and t[0]<0:cell=max(1,int(round(px))-1)
        lo,hi=cell/SCALE[0],(cell+1)/SCALE[0]
        for it in range(maxit):
            side=-1 if constrained and abs(zz[0]*SCALE[0]-(cell+1))<1e-8 else 1
            rr,v=family.evaluate(zz,v,b,side=side);b=rr['beta_r'];F=np.array([rr['beta_i'],np.dot(zz-pred,t)])
            trace.append(dict(point=rr,F=F.tolist()))
            if abs(F[0])<1e-11 and abs(F[1])<1e-9:
                assert rr['real_beta_polynomial_residual']<1e-7
                return zz,rr,v,trace
            dz=np.linalg.solve(np.array([gradient(rr),t]),-F)
            if np.linalg.norm(dz)>.2:dz*=.2/np.linalg.norm(dz)
            factor=1.
            if constrained and dz[0]>0 and zz[0]+dz[0]>hi:factor=.8*(hi-zz[0])/dz[0]
            if constrained and dz[0]<0 and zz[0]+dz[0]<lo:factor=.8*(lo-zz[0])/dz[0]
            assert factor>1e-9,'Corrector left its native interpolation cell'
            zz+=factor*dz
        raise RuntimeError('Neutral corrector iteration limit')
    g=gradient(r);t=np.array([g[1],-g[0]]);t/=np.linalg.norm(t)
    if t[0]<0:t=-t
    startz=z.copy();startq=q.copy();startt=t.copy();startr=r.copy()
    state=dict(case=case,omega_r=omega,f_kHz=a.frequency,component=a.component,scale=SCALE.tolist(),seed=seed,records=[],failures=[],finished=False,status='running',model='Guide: piecewise linear native operator interpolation only')
    def store(z,r,q,t,step,it):
        j=len(state['records']);path=wave/f'q{j:04d}.npy';np.save(path,q)
        state['records'].append(dict(index=j,z=z.tolist(),point=r,tangent=t.tolist(),step=step,iterations=it,vector_file=str(path)))
        if j%10==0:
            tmp=out.with_suffix('.tmp');tmp.write_text(json.dumps(state,indent=2)+'\n');tmp.replace(out)
    step=.03;accumulated=0.;half=0
    if a.resume:
        state=json.loads(out.read_text());state.update(status='running',finished=False);last=state['records'][-1]
        z=np.array(last['z']);q=np.load(last['vector_file']);r=last['point'];t=np.array(last['tangent']);step=.01
        half=int('reverse_half_start_index' in state)
        accumulated=sum(np.linalg.norm(np.array(b['z'])-np.array(c['z'])) for b,c in zip(state['records'][1:],state['records'][:-1]))
    else:
        assert not out.exists();store(z,r,q,t,step,0)
    for loop in range(6000):
        if abs(z[0]*SCALE[0]-round(z[0]*SCALE[0]))<1e-9:
            r,q=family.evaluate(z,q,r['beta_r'],side=np.sign(t[0]));gg=gradient(r);tt=np.array([gg[1],-gg[0]]);tt/=np.linalg.norm(tt)
            if tt[0]*t[0]<0:tt=-tt
            t=tt
        oldz=z.copy();oldq=q.copy();oldsection=float(np.dot(z-startz,startt));accepted=False;edge=False
        for retry in range(14):
            try:
                pred=z+step*t;xi=z[0]*SCALE[0];px=pred[0]*SCALE[0]
                boundary=int(np.floor(xi+1e-8))+1 if t[0]>0 else int(np.ceil(xi-1e-8))-1
                if not 1<=boundary<=family.count:
                    edge=True;break
                bridge=(px-boundary)*t[0]>=0
                if bridge:
                    plane=np.array([1.,0.]);pred=z+(boundary/SCALE[0]-z[0])/t[0]*t
                    zz,rr,v,trace=correct(pred,plane,q,r['beta_r'],30)
                    assert rr['beta_alpha_i']*r['beta_alpha_i']>0,'Native landing skips fold'
                else:zz,rr,v,trace=correct(pred,t,q,r['beta_r'])
                overlap=mac(oldq,v)
                assert overlap>=.95 and np.linalg.norm(zz-oldz)<=2*step,('Arc continuity',overlap)
                gg=gradient(rr);tt=np.array([gg[1],-gg[0]]);tt/=np.linalg.norm(tt)
                if bridge:
                    if tt[0]*t[0]<0:tt=-tt
                    rr['native_grid_corner_bridge']=True
                else:
                    if np.dot(tt,t)<0:tt=-tt
                    assert np.dot(tt,t)>.35,'Sharp turn'
                z,r,q,t=zz,rr,v,tt;accumulated+=np.linalg.norm(z-oldz);r['adjacent_raw_MAC']=overlap
                store(z,r,q,t,step,len(trace));accepted=True;break
            except Exception as e:
                state['failures'].append(dict(index=len(state['records'])-1,step=step,reason=repr(e)));step/=2
                if step<1e-6:break
        if edge:
            state.setdefault('domain_endpoints',[]).append(dict(index=len(state['records'])-1,point=r))
            if half==0:
                half=1;z=startz.copy();q=startq.copy();r=startr.copy();t=-startt.copy();step=.03
                state['reverse_half_start_index']=len(state['records']);store(z,r,q,t,step,0);continue
            state.update(status='open_at_domain_boundaries',finished=True);break
        if not accepted:state.update(status='unresolved_minimum_step',finished=True);break
        section=float(np.dot(z-startz,startt));distance=float(np.linalg.norm(z-startz))
        if half==0 and len(state['records'])>30 and accumulated>.5 and oldsection<0<=section and distance<max(2*step,.02) and np.dot(t,startt)>.7:
            pred=oldz+(-oldsection)/(section-oldsection)*(z-oldz)+1e-3*np.array([-startt[1],startt[0]])
            # Independent return section through start; predictor perturbed along section tangent.
            pred-=np.dot(pred-startz,startt)*startt
            zz,rr,v,trace=correct(pred,startt,q,r['beta_r'])
            closure=dict(scaled_distance=float(np.linalg.norm(zz-startz)),raw_MAC=mac(startq,v),beta_difference=abs(rr['beta_r']-startr['beta_r']),trace=trace)
            closure['passed']=closure['scaled_distance']<1e-7 and closure['raw_MAC']>1-1e-7 and closure['beta_difference']<1e-7
            state['closure']=closure
            if closure['passed']:store(zz,rr,v,startt,0,len(trace));state.update(status='closed_verified',finished=True);break
        if len(state['records'])%25==0:print('ARC',case,a.frequency,a.component,len(state['records']),r['xi'],r['x_mm'],r['beta_r'],r['alpha_r'],step,flush=True)
        if len(trace)<=3 and overlap>.99:step=min(.1,step*1.15)
    if not state['finished']:state.update(status='point_limit',finished=True)
    state['original_inputs']=snapshot();out.write_text(json.dumps(state,indent=2)+'\n');print('FINAL',case,a.frequency,a.component,state['status'],len(state['records']),flush=True)

if __name__=='__main__':main()
