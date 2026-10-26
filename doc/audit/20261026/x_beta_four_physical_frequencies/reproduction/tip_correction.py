"""Physical x folds from Im(alpha)=0 and Im(d alpha/d beta)=0.

Both beta and omega remain real. Adjacent local profiles are reassembled by
START, rather than interpolating an operator and publishing its guide point.
"""
import json
import numpy as np
from scipy.sparse.linalg import norm
from common import ROOT,CASES
from solver import Pencil,original,mac
from physical import Pair

def solve_tip(case,left,right,index,omega,tag):
    guess=left.get('turning_guess')
    if guess is None:
        guess={key:(left['point'][key]+right['point'][key])/2 for key in ['xi','alpha_r','beta_r']}
    xi=float(guess['xi']);beta=float(guess['beta_r']);target=complex(guess['alpha_r']);q=np.load(left['vector_file']);pairs={};trace=[]
    def evaluate(x,b,a,q):
        assert 1<=x<=CASES[case]['stations']
        i=min(int(np.floor(x)),CASES[case]['stations']-1);theta=x-i
        if i not in pairs:pairs[i]=Pair(case,i,omega,tag)
        p=pairs[i].pencil(float(theta));aa,v,r=p.solve_alpha(float(b),a,q,True)
        return np.array([aa.imag,r['alpha_beta_i']]),aa,v,r,i,theta
    F,a,q,r,i,theta=evaluate(xi,beta,target,q)
    for it in range(25):
        trace.append(dict(xi=xi,beta_r=beta,alpha_r=a.real,F=F.tolist(),point=r))
        if max(abs(F))<2e-10:break
        hx=1e-4;hb=1e-3
        Fx,*_=evaluate(xi+hx,beta,a,q);Fb,*_=evaluate(xi,beta+hb,a,q)
        J=np.column_stack([(Fx-F)/hx,(Fb-F)/hb]);delta=np.linalg.solve(J,-F)
        size=max(abs(delta[0])/.5,abs(delta[1])/1.5,1.);delta/=size
        accepted=False;damping=1.;oldnorm=np.linalg.norm(F*np.array([1.,3.]))
        for attempt in range(14):
            xx=xi+damping*delta[0];bb=beta+damping*delta[1]
            if not 1<=xx<=CASES[case]['stations']:damping/=2;continue
            try:FF,aa,v,rr,ii,tt=evaluate(xx,bb,a,q)
            except (AssertionError,RuntimeError):damping/=2;continue
            if np.linalg.norm(FF*np.array([1.,3.]))<oldnorm:
                xi,beta,F,a,q,r,i,theta=xx,bb,FF,aa,v,rr,ii,tt;accepted=True;break
            damping/=2
        assert accepted,('Physical spatial fold corrector stalled',tag,index,trace[-1])
    assert max(abs(F))<1e-8,(tag,index,F)
    assert abs(xi-guess['xi'])<2.,('Physical fold moved too far from guide',xi,guess['xi'])
    pair=pairs[i];generation=json.loads((pair.build(theta)/'generation.json').read_text());endpoint_checks=[]
    for t,j in [(0.,i),(1.,i+1)]:
        physical=pair.pencil(t);native=Pencil(case,j,omega)
        errors=[float(norm(aa-bb))/max(1.,float(norm(aa))) for aa,bb in zip(native.spatial(beta),physical.spatial(beta))]
        assert max(errors)<1e-9,('Profile endpoint mismatch',tag,errors)
        endpoint_checks.append(dict(station=j,relative_matrix_errors=errors))
    state=dict(omega_r=omega,beta_r=beta,alpha_r=a.real,alpha_i=a.imag)
    point,v,receipt=original(pair.build(theta),5,state,beta,f'{tag}_physical_tip{index}_spatial')
    assert abs(point['growth'])<1e-8 and mac(v,q)>1-1e-7
    # Arc order brackets the true local x extremum, even if its physical-profile
    # correction moves slightly relative to the operator guide.
    order=(left['order']+right['order'])/2
    out=dict(case=case,kind=1,xi=float(xi),x_mm=generation['x_mm'],x_flat_mm=generation['x_mm']-CASES[case]['flat_start_mm'],
        alpha_r=point['alpha_r'],alpha_i=point['alpha_i'],beta_r=float(beta),beta_i=0.,omega_r=omega,omega_i=0.,growth=point['growth'],
        right_residual=r['right_residual'],left_residual=r['left_residual'],raw_MAC=mac(v,q),guide_order=order,
        guide_left_index=left['index'],guide_right_index=right['index'],guide_xi=guess['xi'],native_bracket=[i,i+1],theta=float(theta),trace=trace,
        endpoint_checks=endpoint_checks,original_receipt=receipt,generation=str(pair.build(theta)/'generation.json'),
        turning_condition='Im(alpha)=0 and Im(alpha_beta)=0; real beta and omega')
    path=ROOT/'outputs'/f'{tag}_tip{index}.npy';np.save(path,v);out['vector_file']=str(path);return out
