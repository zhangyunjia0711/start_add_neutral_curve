"""One separately authorized native10 restart with actual11/12 x-secant."""
import json,traceback
from pathlib import Path
import h5py,numpy as np
import temporal_pencil as tp
RUN=tp.RUN;HERE=Path(__file__).resolve().parent/'M1c_native'
OUT=RUN/'outputs/matrix_temporal_M1c_native10_secant_20261023.json'
LOG=RUN/'outputs/工作日志_2026-10-23_矩阵频率换元M1c原生鼻部延拓.md'
def save(j):OUT.write_text(json.dumps(j,indent=2,allow_nan=False)+'\n')
def log(t):
    with LOG.open('a') as f:f.write('\n'+tp.stamp()+' '+t+'\n')

def main():
    assert not OUT.exists()
    src=RUN/'outputs/nose_neutral_native_vectors_20261023.json';old=next(v for v in json.loads(src.read_text())['cases'] if v['family']=='M1c')
    past=RUN/'outputs/local_peak_nose_M1c_ridge_v10_1616_progress.json';p12=next(v['result']['point'] for v in json.loads(past.read_text())['records'] if v.get('station')==12)
    p11=old['actual_point'];native=old['actual_native']
    with h5py.File(Path(native['input_dir'])/'grid_vertical2d_N.h5','r') as f:x10=float(f['grid'][0,9,0,0])
    frac=(x10-p11['x'])/(p11['x']-p12['x'])
    predict={k:p11[k]+frac*(p11[k]-p12[k]) for k in ['alpha_r','beta_r','omega_r']}
    q=tp.vf.vectors(Path(native['native_vector_file']),400)[native['raw_index']].ravel();q/=np.linalg.norm(q)
    first=json.loads((RUN/'outputs/matrix_temporal_M1c_native_20261023.json').read_text());b0=first['source_temporal_refinement']['point']['beta'];w0=first['source_temporal_refinement']['point']['omega']['real']
    model=tp.Pencil(HERE/'i10_matrices',native,10,b0,w0);checks=model.build()
    report=dict(created=tp.stamp(),status='running',station=10,phase='separate_authorized_actual_x_secant_restart',
        provenance=[dict(path=str(p),sha256=tp.driver.sha(p)) for p in [src,past]],source11=p11,source12=p12,
        target_x=x10,physical_x_secant_fraction=frac,prediction=predict,matrix_checks=checks,
        limits=dict(max_iterations=20,max_backtracks=4,minimum_raw_MAC=.95,minimum_abs_omega_alpha=1e-9),
        physical_candidate=False)
    save(report);log('独立重策：前一same11初猜失败保留；root明确授权一次实际11/12 x-secant、20iter/4回溯。初值 '+json.dumps(predict))
    try:
        ans,qq=tp.correct(model,predict['alpha_r'],predict['beta_r'],q,predict['omega_r']+1e-3,maxiter=20,backtracks=4)
        report['result']=ans;report['status']='native10_'+ans['status'];save(report)
        if ans['status']=='temporal_candidate':
            t=ans['point'];point=dict(alpha_r=t['alpha'],alpha_i=0.,omega_r=t['omega']['real'],beta_r=t['beta'])
            actual=tp.vf.solve('M1c_i10_secant_matrix_candidate_native',Path(native['input_dir']),400,10,point,here=HERE,binary=tp.BIN,sha=tp.SHA)
            v=tp.vf.vectors(Path(actual['native_vector_file']),400)[actual['raw_index']].ravel()
            report.update(native_fresh=actual,native_matrix_alpha_difference=abs(complex(actual['alpha_r'],actual['alpha_i'])-point['alpha_r']),native_matrix_raw_MAC=tp.mac(qq,v))
            report['native_fresh_passed']=actual['passed'] and report['native_matrix_alpha_difference']<1e-8 and report['native_matrix_raw_MAC']>1-1e-10 and abs(actual['growth'])<1e-7
            np.save(HERE/'i10_secant_temporal_candidate_vector.npy',qq)
            report['status']='native10_fresh_spatial_passed_curvature_pending' if report['native_fresh_passed'] else 'native10_fresh_spatial_failed'
            save(report)
    except Exception:
        report['status']='bounded_secant_probe_failed';report['exception']=traceback.format_exc();log(report['exception']);raise
    finally:
        tp.driver.snapshot();assert tp.driver.sha(tp.BIN)==tp.SHA and all(tp.driver.sha(Path(p['path']))==p['sha256'] for p in report['provenance'])
        report.update(completed=tp.stamp(),original_files_unchanged=True,temporal_calls=model.calls)
        save(report);log('secant阶段结束 '+report['status']+'；空间负曲率未核验前不列candidate_converged，无追加无依据迭代。')
        print('SECANT RESULT '+json.dumps(dict(status=report['status'],prediction=predict,iterations=report.get('result',{}).get('iterations'))),flush=True)

if __name__=='__main__':main()
