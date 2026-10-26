"""Neutral continuation in real (station,alpha,beta), with no eigenvalue elimination.

A bordered original polynomial supplies a complex scalar zero and its exact
parameter derivatives. This remains regular at a fold of an alpha or beta
eigenvalue projection whenever the full two-by-three neutral Jacobian has rank2.
"""
import argparse,json,time
import numpy as np
from scipy.sparse import bmat,csc_matrix
from scipy.sparse.linalg import splu
from common import ROOT,CASES,snapshot
from solver import Pencil,qep,mac
from contour import Family,seed_root
SCALE=np.array([50.,20.,10.])

def bordered(C,dC,z,qref,pref):
    xi,alpha,beta=z*SCALE
    Q=sum(alpha**k*beta**j*C[k][j] for k in range(3) for j in range(3))
    derivatives=[sum(alpha**k*beta**j*dC[k][j] for k in range(3) for j in range(3)),
        sum(beta**j*(C[1][j]+2*alpha*C[2][j]) for j in range(3)),
        sum(alpha**k*(C[k][1]+2*beta*C[k][2]) for k in range(3))]
    n=len(qref);M=bmat([[Q,csc_matrix(pref[:,None])],[csc_matrix(qref.conj()[None,:]),csc_matrix((1,1),dtype=complex)]],format='csc')
    lu=splu(M);rhs=np.r_[np.zeros(n,dtype=complex),1.+0j];solution=lu.solve(rhs);adj=lu.solve(rhs,trans='H');q=solution[:-1];p=adj[:-1];error=solution[-1]
    gradient=np.array([-np.vdot(p,A@q)*scale for A,scale in zip(derivatives,SCALE)])
    J=np.array([gradient.real,gradient.imag]);_,singular,V=np.linalg.svd(J,full_matrices=True)
    assert singular[-1]>1e-10*singular[0],('Neutral Jacobian rank loss',singular,z)
    tangent=V[-1];q/=np.linalg.norm(q);p/=np.linalg.norm(p)
    return error,J,tangent,q,p,float(np.linalg.norm(Q@q)),float(np.linalg.norm(Q.conj().T@p))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('case');parser.add_argument('frequency',type=int,choices=[50,100]);parser.add_argument('component',choices=['A','B']);args=parser.parse_args()
    case,f,component=args.case,args.frequency,args.component;tag=f'{case}_f{f}kHz_{component}';omega=next(r['omega_nd'] for r in json.loads((ROOT/'outputs/frequency_conversion.json').read_text())['frequencies'] if r['f_kHz']==f)
    family=Family(case,omega);seed,q=seed_root(case,f,component);s=seed['point'];pencil=Pencil(case,seed['station'],omega)
    _,q,_,p=qep(pencil.spatial(s['beta_r']),complex(s['alpha_r'],s['alpha_i']),q,True)
    z=np.array([seed['station'],s['alpha_r'],s['beta_r']])/SCALE;out=ROOT/'outputs'/f'{tag}_v3_guide.json';wave=ROOT/'work'/f'{tag}_v3_guide';wave.mkdir(exist_ok=True);assert not out.exists()
    def evaluate(z,q,p,side=1):
        xi=z[0]*SCALE[0];C,dC,i,theta=family.pencil(float(xi),side);e,J,t,v,w,res,lr=bordered(C,dC,z,q,p)
        r=dict(xi=float(xi),x_mm=float((1-theta)*family.x[i-1]+theta*family.x[i]),alpha_r=float(z[1]*SCALE[1]),alpha_i=0.,beta_r=float(z[2]*SCALE[2]),beta_i=0.,omega_r=omega,omega_i=0.,growth=0.,
            right_residual=res,left_residual=lr,raw_MAC=mac(q,v),border_error_r=float(e.real),border_error_i=float(e.imag),native_bracket=[i,i+1],theta=float(theta))
        return e,J,t,v,w,r
    def correct(pred,plane,q,p,maxit=20):
        zz=pred.copy();tr=[];constrained=np.linalg.norm(plane[1:])>1e-12;px=pred[0]*SCALE[0];cell=min(int(np.floor(px)),family.count-1)
        if abs(px-round(px))<1e-9 and plane[0]<0:cell=max(1,int(round(px))-1)
        lo,hi=cell/SCALE[0],(cell+1)/SCALE[0]
        for it in range(maxit):
            side=-1 if constrained and abs(zz[0]*SCALE[0]-(cell+1))<1e-8 else 1
            e,J,t,v,w,r=evaluate(zz,q,p,side);F=np.array([e.real,e.imag,np.dot(zz-pred,plane)]);tr.append(dict(point=r,F=F.tolist()))
            if abs(e)<1e-12 and abs(F[2])<1e-9 and r['right_residual']<1e-8:return zz,r,v,w,t,tr
            A=np.vstack([J,plane]);weights=1/np.maximum(np.linalg.norm(A,axis=1),1e-30);dz=np.linalg.solve(A*weights[:,None],-F*weights)
            if np.linalg.norm(dz)>.2:dz*=.2/np.linalg.norm(dz)
            factor=1.
            if constrained and dz[0]>0 and zz[0]+dz[0]>hi:factor=.8*(hi-zz[0])/dz[0]
            if constrained and dz[0]<0 and zz[0]+dz[0]<lo:factor=.8*(lo-zz[0])/dz[0]
            assert factor>1e-9,'3D corrector exits native cell';zz+=factor*dz
        raise RuntimeError('3D neutral corrector iteration limit')
    e,J,t,q,p,r=evaluate(z,q,p);assert abs(e)<1e-10
    if t[0]<0:t=-t
    z0,q0,p0,t0,r0=z.copy(),q.copy(),p.copy(),t.copy(),dict(r);step=.025
    state=dict(case=case,f_kHz=f,component=component,omega_r=omega,scale=SCALE.tolist(),seed=seed,records=[],failures=[],status='running',finished=False,model='Guide only: bordered original real-alpha/real-beta polynomial; piecewise native-operator interpolation')
    def store(z,r,q,p,t,ds,it):
        j=len(state['records']);qfile=wave/f'q{j:04d}.npy';pfile=wave/f'p{j:04d}.npy';np.save(qfile,q);np.save(pfile,p)
        state['records'].append(dict(index=j,z=z.tolist(),point=r,tangent=t.tolist(),step=ds,iterations=it,vector_file=str(qfile),left_vector_file=str(pfile)))
        if j%10==0:
            tmp=out.with_suffix('.tmp');tmp.write_text(json.dumps(state,indent=2)+'\n');tmp.replace(out)
    store(z,r,q,p,t,step,0);arc=0.
    for loop in range(4000):
        if abs(z[0]*SCALE[0]-round(z[0]*SCALE[0]))<1e-9:
            e,J,tt,q,p,r=evaluate(z,q,p,side=np.sign(t[0]));
            if tt[0]*t[0]<0:tt=-tt
            t=tt
        oldz=z.copy();oldq=q.copy();oldsection=float(np.dot(z-z0,t0));accepted=False
        for retry in range(16):
            try:
                pred=z+step*t;xi=z[0]*SCALE[0];px=pred[0]*SCALE[0];boundary=int(np.floor(xi+1e-8))+1 if t[0]>0 else int(np.ceil(xi-1e-8))-1
                assert 1<=boundary<=family.count,'Neutral curve reaches domain edge'
                bridge=(px-boundary)*t[0]>=0
                if bridge:
                    plane=np.array([1.,0.,0.]);pred=z+(boundary/SCALE[0]-z[0])/t[0]*t
                    zz,rr,v,w,tt,tr=correct(pred,plane,q,p,30)
                    if tt[0]*t[0]<0:tt=-tt
                else:
                    zz,rr,v,w,tt,tr=correct(pred,t,q,p)
                    if np.dot(tt,t)<0:tt=-tt
                    assert np.dot(tt,t)>.35,'3D sharp turn'
                overlap=mac(oldq,v);assert overlap>.95 and np.linalg.norm(zz-oldz)<2*step,('3D continuity gate',overlap)
                z,r,q,p,t=zz,rr,v,w,tt;arc+=np.linalg.norm(z-oldz);store(z,r,q,p,t,step,len(tr));accepted=True;break
            except Exception as error:
                state['failures'].append(dict(index=len(state['records'])-1,step=step,reason=repr(error)));step/=2
                if step<1e-6:break
        if not accepted:state.update(status='unresolved_minimum_step',finished=True);break
        section=float(np.dot(z-z0,t0));distance=float(np.linalg.norm(z-z0))
        if len(state['records'])>30 and arc>.5 and oldsection<0<=section and distance<max(2*step,.03) and np.dot(t,t0)>.7:
            normal=np.array([1.,0.,0.]);normal-=np.dot(normal,t0)*t0;normal/=np.linalg.norm(normal)
            pred=oldz+(-oldsection)/(section-oldsection)*(z-oldz)+1e-3*normal;pred-=np.dot(pred-z0,t0)*t0
            zz,rr,v,w,tt,tr=correct(pred,t0,q,p)
            closure=dict(scaled_distance=float(np.linalg.norm(zz-z0)),raw_MAC=mac(q0,v),trace=tr,initial_distance=float(np.linalg.norm(pred-z0)))
            closure['passed']=closure['scaled_distance']<1e-7 and closure['raw_MAC']>1-1e-7;state['closure']=closure
            if closure['passed']:store(zz,rr,v,w,t0,0,len(tr));state.update(status='closed_verified',finished=True);break
        if len(state['records'])%50==0:print('ARC3',tag,len(state['records']),r['x_mm'],r['alpha_r'],r['beta_r'],flush=True)
        if len(tr)<=4 and overlap>.99:step=min(.08,step*1.15)
    if not state['finished']:state.update(status='point_limit',finished=True)
    state['original_inputs']=snapshot();out.write_text(json.dumps(state,indent=2)+'\n');print('FINAL3',tag,state['status'],len(state['records']),flush=True)
if __name__=='__main__':main()
