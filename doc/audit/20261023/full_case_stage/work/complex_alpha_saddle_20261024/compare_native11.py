"""Native i11 counterpart, with the same bounded complex-alpha algorithm."""
import json
from pathlib import Path
import numpy as np
import saddle_probe as sa
import verify_neighbors as vn
tp=sa.tp;RUN=sa.RUN;HERE=sa.HERE
OUT=RUN/'outputs/complex_alpha_saddle_native11_20261024.json'

def main():
    assert not OUT.exists()
    path=RUN/'outputs/nose_neutral_native_vectors_20261023.json';src=next(x for x in json.loads(path.read_text())['cases'] if x['family']=='M1c');native=src['actual_native'];p=src['actual_point']
    model=tp.Pencil(RUN/'work/pep_matrix_temporal_probe_20261023/M1c_native/i11_matrices',native,11,p['beta_r'],p['omega_r']);checks=model.build()
    q=tp.vf.vectors(Path(native['native_vector_file']),400)[native['raw_index']].ravel();q/=np.linalg.norm(q)
    a=complex(native['alpha_r'],native['alpha_i']);b=p['beta_r'];w=p['omega_r']
    report=dict(created=tp.stamp(),status='running',station=11,fixed_beta=b,source_path=str(path),source_sha256=tp.driver.sha(path),source_actual=native,
        matrix_checks=checks,solver_sha256=tp.SHA,limits=dict(max_iterations=15,max_backtracks=4,minimum_MAC=.95),
        pinch_origin_verified=False,physical_absolute_instability_claim=False,
        qualification='Native11 fixed-beta finite saddle counterpart only; no causal half-plane-origin proof')
    OUT.write_text(json.dumps(report,indent=2)+'\n');tp.driver.snapshot()
    parity,qq,_=sa.evaluate(model,a,b,q,w+1e-3);report['source_complex_alpha_parity']=dict(point=parity,frequency_error=abs(sa.z(parity,'omega')-w),raw_MAC=tp.mac(q,qq))
    assert report['source_complex_alpha_parity']['frequency_error']<1e-8
    ans,qs=sa.saddle(model,a,b,qq,w,maxiter=15,backtracks=4);report['result']=ans;report['status']=ans['status'];OUT.write_text(json.dumps(report,indent=2)+'\n')
    if ans['status']=='complex_alpha_saddle_candidate':
        r=ans['point'];a0=sa.z(r,'alpha');w0=sa.z(r,'omega');report['h_h2_verification']=[]
        for h in [1e-3,5e-4]:
            v=[]
            for da in [-h,h,-1j*h,1j*h]:
                rr,_,_=sa.evaluate(model,a0+da,b,qs,w0+1e-3);v.append(rr)
            wr=(sa.z(v[1],'omega_alpha')-sa.z(v[0],'omega_alpha'))/(2*h);wi=(sa.z(v[3],'omega_alpha')-sa.z(v[2],'omega_alpha'))/(2j*h)
            report['h_h2_verification'].append(dict(h=h,omega_alphalpha=tp.cn(wr),imag_direction_omega_alphalpha=tp.cn(wi),direction_relative_difference=abs(wr-wi)/abs(wr),samples=v))
        np.save(HERE/'native11_saddle_q.npy',qs);report['saddle_vector_file']=str(HERE/'native11_saddle_q.npy');report['saddle_vector_sha256']=tp.driver.sha(HERE/'native11_saddle_q.npy')
        waa=sa.z(report['h_h2_verification'][-1],'omega_alphalpha');report['native_neighbor_pairs']=[]
        # One frequency supplies two independent original-PEP residual-qualified
        # branches; native10 already has the three-scale convergence series.
        eps=1e-5;omega=w0+1j*eps;da=np.sqrt(2j*eps/waa)
        for sign in [-1,1]:
            pair,_=vn.native_pair(HERE/'runs'/('native11_eps1e-5_s%+d'%sign),native,a0+sign*da,omega,b,qs,model);report['native_neighbor_pairs'].append(pair)
        roots=[sa.z(v,'actual_alpha') for v in report['native_neighbor_pairs']];assert abs(roots[0]-roots[1])>1e-5
        report['native_neighbor_diagnostic']=dict(epsilon=eps,root_separation=abs(roots[0]-roots[1]),midpoint_difference_to_saddle=abs((roots[0]+roots[1])/2-a0),
            square_law_relative_errors=[abs(.5*waa*(v-a0)**2-1j*eps)/eps for v in roots])
    tp.driver.snapshot();assert tp.driver.sha(tp.BIN)==tp.SHA and tp.driver.sha(path)==report['source_sha256']
    report.update(completed=tp.stamp(),original_files_unchanged=True);OUT.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    sa.log('native11有界对照 '+json.dumps(dict(status=report['status'],point=report.get('result',{}).get('point')))+'。同样未证pinch起源，不作absolute分类。')
    print('SADDLE11 '+json.dumps(dict(status=report['status'],point=report.get('result',{}).get('point'),iterations=report.get('result',{}).get('iterations'))),flush=True)

if __name__=='__main__':main()
