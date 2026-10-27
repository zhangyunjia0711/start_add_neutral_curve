"""Resolve a physical x fold without eliminating the alpha eigenvalue.

At a real neutral point the two real bordered equations have a 2x3 Jacobian
in (x,alpha,beta). An x tangent vanishes when its (alpha,beta) minor vanishes.
This avoids division by p^H P_alpha q near an alpha projection fold.
"""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.sparse import bmat,csc_matrix
from scipy.sparse.linalg import splu,norm
from common import ROOT,CASES
from solver import Pencil,qep,original,mac
from physical import Pair

def solve_tip(case,left,right,index,omega,tag):
    guess=left.get('turning_guess') or {k:(left['point'][k]+right['point'][k])/2 for k in ['xi','alpha_r','beta_r']}
    xi=float(guess['xi']);beta=float(guess['beta_r']);alpha=float(guess['alpha_r']);reference=np.load(left['vector_file']);pairs={};trace=[]
    def pencil(x):
        assert 1<=x<=CASES[case]['stations']
        i=min(int(np.floor(x)),CASES[case]['stations']-1);theta=x-i
        if i not in pairs:pairs[i]=Pair(case,i,omega,tag)
        return pairs[i].pencil(float(theta)),i,theta
    native,i,theta=pencil(xi)
    aa,q0,r,p0=qep(native.spatial(beta),complex(alpha),reference,True)
    assert mac(q0,reference)>.95
    n=len(q0);rhs=np.r_[np.zeros(n,dtype=complex),1.+0j]
    def evaluate(z):
        x,a,b=z*np.array([1.,20.,10.]);p,i,t=pencil(float(x))
        C=p.C;Q=sum(a**k*b**j*C[k][j] for k in range(3) for j in range(3))
        M=bmat([[Q,csc_matrix(p0[:,None])],[csc_matrix(q0.conj()[None,:]),csc_matrix((1,1),dtype=complex)]],format='csc')
        lu=splu(M);v=lu.solve(rhs);w=lu.solve(rhs,trans='H')
        da=sum(b**j*(C[1][j]+2*a*C[2][j]) for j in range(3))
        db=sum(a**k*(C[k][1]+2*b*C[k][2]) for k in range(3))
        ga=-np.vdot(w[:-1],da@v[:-1])*20;gb=-np.vdot(w[:-1],db@v[:-1])*10
        determinant=ga.real*gb.imag-ga.imag*gb.real
        q=v[:-1]/np.linalg.norm(v[:-1]);res=float(np.linalg.norm(Q@q))
        scaled=res/(max(1.,sum(abs(a)**k*norm(K) for k,K in enumerate(p.spatial(b)))))
        return np.array([v[-1].real,v[-1].imag,determinant]),q,i,t,res,scaled,np.array([ga,gb])
    z0=np.array([xi,alpha/20,beta/10]);F,q,i,t,_,_,g=evaluate(z0)
    scales=np.array([max(np.linalg.norm(g.real),1e-8),max(np.linalg.norm(g.imag),1e-8),max(abs(g[0])*abs(g[1]),1e-8)])
    def fun(z):
        F,q,i,t,res,scaled,g=evaluate(z);rr=F/scales
        trace.append(dict(xi=float(z[0]),alpha_r=float(20*z[1]),beta_r=float(10*z[2]),F=F.tolist(),scaled_F=rr.tolist(),original_matrix_residual=res,coefficient_scaled_residual=scaled,raw_MAC=mac(q0,q)))
        assert mac(q0,q)>.90
        return rr
    lower=np.array([max(1.,xi-1.99),alpha/20-2.,beta/10-2.])
    upper=np.array([min(CASES[case]['stations'],xi+1.99),alpha/20+2.,beta/10+2.])
    opt=least_squares(fun,z0,bounds=(lower,upper),diff_step=1e-5,x_scale='jac',ftol=1e-12,xtol=1e-12,gtol=1e-12,max_nfev=100)
    F,q,i,theta,res,scaled,g=evaluate(opt.x);xi=float(opt.x[0]);alpha=float(20*opt.x[1]);beta=float(10*opt.x[2])
    assert np.max(abs(F/scales))<1e-7 and res<1e-7,('Bordered physical fold unresolved',tag,index,opt.message,F/scales,trace[-1])
    pair=pairs[i];directory=pair.build(theta);generation=json.loads((directory/'generation.json').read_text())
    checks=[]
    for t,j in [(0.,i),(1.,i+1)]:
        a,b=pair.pencil(t),Pencil(case,j,omega)
        errors=[float(norm(x-y))/max(1.,float(norm(x))) for x,y in zip(a.spatial(beta),b.spatial(beta))]
        assert max(errors)<1e-9;checks.append(dict(station=j,relative_matrix_errors=errors))
    state=dict(omega_r=omega,beta_r=beta,alpha_r=alpha,alpha_i=0.)
    point,v,receipt=original(directory,5,state,beta,f'{tag}_physical_bordered_tip{index}_spatial')
    assert abs(point['growth'])<1e-8 and mac(v,q)>1-1e-7
    qpath=ROOT/'outputs'/f'{tag}_tip{index}.npy';np.save(qpath,v)
    return dict(case=case,kind=1,xi=xi,x_mm=generation['x_mm'],x_flat_mm=generation['x_mm']-CASES[case]['flat_start_mm'],
        alpha_r=point['alpha_r'],alpha_i=point['alpha_i'],beta_r=beta,beta_i=0.,omega_r=omega,omega_i=0.,growth=point['growth'],
        right_residual=res,left_residual=None,coefficient_scaled_residual=scaled,raw_MAC=mac(v,q),guide_order=(left['order']+right['order'])/2,
        guide_left_index=left['index'],guide_right_index=right['index'],guide_xi=guess['xi'],native_bracket=[i,i+1],theta=float(theta),trace=trace,
        endpoint_checks=checks,original_receipt=receipt,generation=str(directory/'generation.json'),vector_file=str(qpath),
        turning_condition='Re(border_error)=Im(border_error)=0; real(alpha,beta) Jacobian minor=0',bordered_scaled_fold_residual=(F/scales).tolist())
