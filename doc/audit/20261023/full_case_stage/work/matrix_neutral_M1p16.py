#!/usr/bin/env python3
"""Bounded no-core native17 validation and native16 temporal-pencil candidate."""
import json,sys,time
from pathlib import Path
import numpy as np
import full_driver as driver
from local_peak_x_continue import stamp

ROOT=driver.ROOT
sys.path.insert(0,str(ROOT/'work/pep_matrix_temporal_probe_20261023'))
import temporal_pencil as tp

def z(d):return complex(d['real'],d['imag'])

def derivatives(pencil,point,q,label):
    alpha,beta=point['alpha'],point['beta'];omega=z(point['omega'])
    hs=[max(1e-4,abs(alpha)*1e-5),max(1e-5,abs(beta)*1e-5)]
    output=[]
    for factor in [1.,.5]:
        J=np.empty((2,2));pairs=[];fd=[]
        for axis,h in enumerate(hs):
            h*=factor;records=[]
            for sign in [-1,1]:
                args=[alpha,beta];args[axis]+=sign*h
                r,_,_=pencil.evaluate(*args,q,omega+1e-3,left=True);records.append(r)
            J[:,axis]=(tp.objective(records[1])-tp.objective(records[0]))/(2*h)
            fd.append((z(records[1]['omega'])-z(records[0]['omega']))/(2*h));pairs.append(records)
        F=tp.objective(point);delta=np.linalg.solve(J,-F)
        ab=z(point['implicit_alpha_beta']).real
        curvature=float(J[1,1]+ab*J[1,0])
        wa,wb=z(point['omega_alpha']),z(point['omega_beta'])
        output.append(dict(label=label,factor=factor,h=[h*factor for h in hs],J=J.tolist(),F=F.tolist(),
            condition=float(np.linalg.cond(J)),alpha_beta_real=ab,spatial_beta_curvature=curvature,
            proposed_alpha_beta_correction=delta.tolist(),
            proposed_omega_correction=tp.cn(wa*delta[0]+wb*delta[1]),
            analytic_omega_alpha=tp.cn(wa),analytic_omega_beta=tp.cn(wb),
            finite_difference_omega_alpha=tp.cn(fd[0]),finite_difference_omega_beta=tp.cn(fd[1]),
            alpha_derivative_relative_error=abs(fd[0]-wa)/max(abs(wa),1e-30),
            beta_derivative_relative_error=abs(fd[1]-wb)/max(abs(wb),1e-30),stencils=pairs))
    return output

def main():
    out=ROOT/'outputs/matrix_neutral_M1p16_20261023.json'
    if out.exists():raise RuntimeError('Existing matrix phase preserved')
    source_path=ROOT/'outputs/nose_neutral_native_vectors_20261023.json'
    source=next(v for v in json.loads(source_path.read_text())['cases'] if v['family']=='M1p')
    native=source['actual_native'];point=source['source_point']
    init_path=ROOT/'outputs/matrix_M1p16_initialization_20261023.json'
    initial=json.loads(init_path.read_text())['predicted_native16']
    q=tp.vf.vectors(Path(native['native_vector_file']),400)[native['raw_index']].ravel()
    here=ROOT/'work/matrix_neutral_M1p16_20261023';here.mkdir(parents=True,exist_ok=True)
    start=time.monotonic();report=dict(created=stamp(),status='running',source_path=str(source_path),
        source_sha256=driver.sha(source_path),initialization_path=str(init_path),initialization_sha256=driver.sha(init_path),
        solver=str(tp.BIN),solver_sha256=tp.SHA,helper_path=str(Path(tp.__file__)),helper_sha256=driver.sha(Path(tp.__file__)),
        wrapper_sha256=driver.sha(Path(__file__)),qualification='Finite native spatial-neutral beta-stationary candidate, not a complete spectrum or a global beta envelope')
    def save():
        report.update(updated=stamp(),seconds=time.monotonic()-start)
        temp=out.with_suffix('.tmp');temp.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');temp.replace(out)
    driver.snapshot();save()
    try:
        p17=tp.Pencil(here/'native17_matrices',native,17,point['beta_r'],point['omega_r'])
        report['native17_build']=p17.build();save()
        r17,q17,_=p17.evaluate(native['alpha_r'],point['beta_r'],q,point['omega_r']+1e-3,left=True)
        report['native17_reference']=r17;save()
        diagnostics=derivatives(p17,r17,q17,'native17_reference');report['native17_h_h2']=diagnostics
        omega_error=abs(z(r17['omega'])-point['omega_r'])
        conditions=bool(omega_error<1e-7 and all(d['spatial_beta_curvature']<0 and
            abs(d['proposed_alpha_beta_correction'][1])<1e-5 and abs(z(d['proposed_omega_correction']))<1e-4 and
            d['alpha_derivative_relative_error']<.01 and d['beta_derivative_relative_error']<.01 for d in diagnostics))
        report['native17_validation']=dict(omega_difference=omega_error,passed=conditions,calls=p17.calls);save()
        if not conditions:raise RuntimeError('Known native17 temporal/derivative validation did not pass')
        p16=tp.Pencil(here/'native16_matrices',native,16,initial['beta_r'],initial['omega_r'])
        report['native16_build']=p16.build();save()
        r16,q16,_=p16.evaluate(initial['alpha_r'],initial['beta_r'],q17,initial['omega_r']+1e-3,left=True)
        report['native16_initial_actual_temporal']=r16;save()
        result,finalq=tp.correct(p16,initial['alpha_r'],initial['beta_r'],q16,
            z(r16['omega'])+1e-3,maxiter=20,backtracks=4)
        report['native16_correction']=result;report['native16_temporal_calls']=p16.calls;save()
        if result['status']!='temporal_candidate':
            report['status']='native16_unknown_bounded_matrix_corrector_stopped';save();return
        candidate=result['point'];diagnostics=derivatives(p16,candidate,finalq,'native16_candidate')
        report['native16_h_h2']=diagnostics;save()
        if not all(d['spatial_beta_curvature']<0 and abs(d['proposed_alpha_beta_correction'][1])<1e-5 and
            abs(z(d['proposed_omega_correction']))<1e-4 for d in diagnostics):
            report['status']='native16_stationary_candidate_not_verified_beta_maximum';save();return
        w=z(candidate['omega']);expected=dict(omega_r=w.real,beta_r=candidate['beta'],alpha_r=candidate['alpha'],alpha_i=0.)
        fresh=tp.vf.solve('native16_spatial_candidate_fresh',driver.CASE/'M1/in',400,16,expected,
            here=here,binary=tp.BIN,sha=tp.SHA)
        report['native16_spatial_fresh']=fresh
        if fresh.get('native_vector_file'):
            sv=tp.vf.vectors(Path(fresh['native_vector_file']),400)[fresh['raw_index']].ravel()
            rawmac=tp.mac(finalq,sv)
        else:rawmac=None
        delta=abs(complex(fresh.get('alpha_r',1e100)-candidate['alpha'],fresh.get('alpha_i',1e100)))
        passed=bool(fresh['passed'] and abs(fresh['growth'])<1e-7 and delta<1e-7 and rawmac>=1-1e-8)
        report['native16_independent_verification']=dict(alpha_difference=delta,raw_primitive_MAC=rawmac,passed=passed)
        if passed:
            point=dict(omega_r=w.real,beta_r=candidate['beta'],alpha_r=fresh['alpha_r'],alpha_i=fresh['alpha_i'],
                growth=fresh['growth'],iloc=16,x=float(fresh['metadata']['x_global_N']))
            report.update(status='native16_beta_stationary_neutral_candidate_verified',point=point)
        else:report['status']='native16_temporal_candidate_failed_spatial_verification'
        save()
    except Exception as error:report.update(status='stopped_unresolved',reason=repr(error));save()
    finally:driver.snapshot();report['finished']=stamp();save()

if __name__=='__main__':main()
