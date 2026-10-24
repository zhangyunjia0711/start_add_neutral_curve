"""Finite complex-alpha saddle probe on the actual exported spatial PEP.

No core/CFD edits; complex temporal simple branch is tracked by actual raw q.
Finding a saddle and two nearby spatial roots does NOT establish pinching.
"""
import json,sys,time,traceback
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import LinearOperator,eigs,splu
HERE=Path(__file__).resolve().parent;RUN=HERE.parent.parent
sys.path.insert(0,str(RUN/'work/pep_matrix_temporal_probe_20261023'))
import temporal_pencil as tp
LOG=RUN/'outputs/工作日志_2026-10-24_空间根方向与复alpha鞍点核验.md'
OUT=RUN/'outputs/complex_alpha_saddle_native10_20261024.json'

def log(t):
    with LOG.open('a') as f:f.write('\n'+tp.stamp()+' '+t+'\n')
def save(j):OUT.write_text(json.dumps(j,indent=2,allow_nan=False)+'\n')
def z(r,k):return complex(r[k]['real'],r[k]['imag'])
def evaluate(model,alpha,beta,qref,omega_shift):
    t=time.time();K,Ka,Kb=model.matrices(alpha,beta);M=model.M
    shift=complex(omega_shift);lu=splu((K-shift*M).tocsc())
    op=LinearOperator(K.shape,matvec=lambda v:lu.solve(M@v),dtype=np.complex128)
    mu,qq=eigs(op,k=2,which='LM',v0=qref,ncv=12,tol=1e-12,maxiter=1000)
    w=shift+1/mu;idx=max(range(len(mu)),key=lambda k:tp.mac(qref,qq[:,k]));q=qq[:,idx];q/=np.linalg.norm(q)
    omega=w[idx];res=float(np.linalg.norm((K-omega*M)@q));mac=tp.mac(qref,q)
    assert res<1e-7 and mac>=.95,('right eigengates',res,mac)
    oph=LinearOperator(K.shape,matvec=lambda v:M.conj().T@lu.solve(v,trans='H'),dtype=np.complex128)
    lm,lv=eigs(oph,k=2,which='LM',ncv=12,tol=1e-12,maxiter=1000)
    li=int(np.argmin(abs(lm-np.conj(mu[idx]))));p=lu.solve(lv[:,li],trans='H');p/=np.linalg.norm(p)
    lres=float(np.linalg.norm((K-omega*M).conj().T@p));den=np.vdot(p,M@q)
    assert lres<1e-7 and abs(den)>1e-14,('left eigengates',lres,abs(den))
    wa=np.vdot(p,Ka@q)/den
    r=dict(alpha=tp.cn(alpha),beta=beta,omega=tp.cn(omega),omega_alpha=tp.cn(wa),
        raw_MAC_to_reference=mac,right_absolute_pencil_residual=res,left_absolute_pencil_residual=lres,
        pH_M_q=tp.cn(den),temporal_condition_indicator=1/abs(den),seconds=time.time()-t,
        shift=tp.cn(shift),other_temporal_values=[tp.cn(v) for k,v in enumerate(w) if k!=idx])
    return r,q,p

def saddle(model,alpha,beta,q,omega,maxiter=15,backtracks=4):
    trace=[];currentq=q;shift=omega+1e-3
    for it in range(maxiter):
        r,v,_=evaluate(model,alpha,beta,currentq,shift);wa=z(r,'omega_alpha');h=1e-3
        samples=[]
        for da in [-h,h]:
            rr,_,_=evaluate(model,alpha+da,beta,v,z(r,'omega')+1e-3);samples.append(rr)
        waa=(z(samples[1],'omega_alpha')-z(samples[0],'omega_alpha'))/(2*h)
        assert abs(waa)>1e-10,'unresolved saddle second derivative'
        delta=-wa/waa;record=dict(iteration=it,center=r,omega_alphalpha_FD=tp.cn(waa),h=h,samples=samples,
            Newton_delta=tp.cn(delta),trials=[]);trace.append(record)
        if abs(wa)<1e-9 and abs(delta)<1e-6:
            return dict(status='complex_alpha_saddle_candidate',point=r,trace=trace,iterations=it+1,
                absolute_omega_alpha_tolerance=1e-9,alpha_correction_tolerance=1e-6),v
        delta*=min(1.,.5/max(abs(delta),1e-300));accepted=False
        for bt in range(backtracks):
            factor=.5**bt
            try:
                rr,vv,_=evaluate(model,alpha+factor*delta,beta,v,z(r,'omega')+1e-3)
                record['trials'].append(dict(factor=factor,point=rr))
                if abs(z(rr,'omega_alpha'))<abs(wa):
                    alpha+=factor*delta;currentq=vv;shift=z(rr,'omega')+1e-3;accepted=True;break
            except Exception as exc:record['trials'].append(dict(factor=factor,error=str(exc)))
        if not accepted:return dict(status='unknown_saddle_no_descent',trace=trace,iterations=it+1),v
    return dict(status='unknown_saddle_max_iterations',trace=trace,iterations=maxiter),currentq

def main():
    assert not OUT.exists(),'inspect saved bounded run before any retry'
    source_path=RUN/'outputs/matrix_temporal_M1c_native10_secant_20261023.json';source=json.loads(source_path.read_text());native=source['native_fresh'];p=source['result']['point']
    oldpath=RUN/'outputs/matrix_temporal_M1c_native_20261023.json';old=json.loads(oldpath.read_text());bp=old['source_temporal_refinement']['point']
    model=tp.Pencil(RUN/'work/pep_matrix_temporal_probe_20261023/M1c_native/i10_matrices',old['source']['actual_native'],10,bp['beta'],bp['omega']['real'])
    checks=model.build();q=tp.vf.vectors(Path(native['native_vector_file']),400)[native['raw_index']].ravel();q/=np.linalg.norm(q)
    a=complex(native['alpha_r'],native['alpha_i']);b=native['fixed_parameters']['beta_r'];w=native['fixed_parameters']['omega_r']
    report=dict(created=tp.stamp(),status='running',station=10,fixed_beta=b,
        source_path=str(source_path),source_sha256=tp.driver.sha(source_path),source_actual=native,
        matrix_provenance_path=str(oldpath),matrix_provenance_sha256=tp.driver.sha(oldpath),matrix_checks=checks,
        solver_sha256=tp.SHA,limits=dict(max_iterations=15,max_backtracks=4,minimum_MAC=.95),
        physical_absolute_instability_claim=False,pinch_origin_verified=False,
        qualification='Finite fixed-real-beta complex-alpha saddle search on one temporal simple branch. Neither Re omega_alpha nor saddle sign alone determines causal direction or absolute/convective class.')
    save(report);tp.driver.snapshot();assert tp.driver.sha(tp.BIN)==tp.SHA
    try:
        initial,qi,_=evaluate(model,a,b,q,w+1e-3)
        report['source_complex_alpha_parity']=dict(point=initial,frequency_error=abs(z(initial,'omega')-w),raw_MAC=tp.mac(q,qi))
        assert report['source_complex_alpha_parity']['frequency_error']<1e-8 and report['source_complex_alpha_parity']['raw_MAC']>1-1e-10
        result,qs=saddle(model,a,b,qi,w);report['result']=result;report['status']=result['status'];save(report)
        if result['status']=='complex_alpha_saddle_candidate':
            r=result['point'];aa=z(r,'alpha');ww=z(r,'omega');verification=[]
            for h in [1e-3,5e-4,1e-4]:
                vals=[]
                for da in [-h,h,-1j*h,1j*h]:
                    v,_,_=evaluate(model,aa+da,b,qs,ww+1e-3);vals.append(v)
                waa_r=(z(vals[1],'omega_alpha')-z(vals[0],'omega_alpha'))/(2*h)
                waa_i=(z(vals[3],'omega_alpha')-z(vals[2],'omega_alpha'))/(2j*h)
                wa_r=(z(vals[1],'omega')-z(vals[0],'omega'))/(2*h)
                wa_i=(z(vals[3],'omega')-z(vals[2],'omega'))/(2j*h)
                verification.append(dict(h=h,omega_alphalpha_real_direction=tp.cn(waa_r),omega_alphalpha_imag_direction=tp.cn(waa_i),
                    second_derivative_direction_relative_difference=abs(waa_r-waa_i)/max(abs(waa_r),1e-30),
                    omega_alpha_FD_real=tp.cn(wa_r),omega_alpha_FD_imag=tp.cn(wa_i),samples=vals))
            report['saddle_h_h2_verification']=verification;report['saddle_vector_file']=str(HERE/'native10_saddle_q.npy');np.save(Path(report['saddle_vector_file']),qs)
            report['saddle_vector_sha256']=tp.driver.sha(Path(report['saddle_vector_file']))
            report['local_second_derivative_verification_passed']=max(v['second_derivative_direction_relative_difference'] for v in verification)<1e-3
            save(report);log('native10实际鞍点 '+json.dumps(r)+'；h/h2 Cauchy方向检查='+str(report['local_second_derivative_verification_passed'])+'。尚无pinch两半平面起源证据，不能称absolute。')
    except Exception:
        report['status']='bounded_saddle_stage_failed';report['exception']=traceback.format_exc();log(report['exception']);raise
    finally:
        tp.driver.snapshot();assert tp.driver.sha(tp.BIN)==tp.SHA and tp.driver.sha(source_path)==report['source_sha256'] and tp.driver.sha(oldpath)==report['matrix_provenance_sha256']
        report.update(completed=tp.stamp(),original_files_unchanged=True);save(report)
        print('SADDLE '+json.dumps(dict(status=report['status'],point=report.get('result',{}).get('point'),iterations=report.get('result',{}).get('iterations'))),flush=True)

if __name__=='__main__':main()
