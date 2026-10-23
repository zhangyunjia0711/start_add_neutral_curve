"""Bounded M1c matrix-coordinate validation and native i10 neutral attempt."""
import datetime,json,sys,time,traceback
from pathlib import Path
import numpy as np
import temporal_pencil as tp
RUN=tp.RUN;HERE=Path(__file__).resolve().parent/'M1c_native';HERE.mkdir(exist_ok=True)
OUT=RUN/'outputs/matrix_temporal_M1c_native_20261023.json'
LOG=RUN/'outputs/工作日志_2026-10-23_矩阵频率换元M1c原生鼻部延拓.md'
def log(t):
    with LOG.open('a') as f:f.write('\n'+tp.stamp()+' '+t+'\n')
def save(r):OUT.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
def z(r,k):return complex(r[k]['real'],r[k]['imag'])

def main():
    assert not OUT.exists(),'inspect saved state; no blind repeat'
    LOG.write_text('# 矩阵频率换元M1c原生鼻部延拓\n\n目标：在不改core/原H5的约束下，以实际空间PEP系数作数值换元，先独立核验beta二次重建/G独立性及左右向量解析导数，再最多20迭代、4回溯实际尝试native10。成功必须原H5 fresh spatial α/向量/残差和实际beta负曲率复核；辅助/矩阵点未通过前不列物理中性结果。单线程；所有失败和边界保留。原i1..10未知不会由换元假定补齐。理论公式为我们对简单解析分支的隐函数代数推断，不用Gaster近似；接口依据见前PEP矩阵小试及本任务实际矩阵验证。\n')
    tp.driver.snapshot();assert tp.driver.sha(tp.BIN)==tp.SHA
    srcfile=RUN/'outputs/nose_neutral_native_vectors_20261023.json'
    source=next(v for v in json.loads(srcfile.read_text())['cases'] if v['family']=='M1c')
    native=source['actual_native'];p=source['actual_point'];q=tp.vf.vectors(Path(native['native_vector_file']),400)[native['raw_index']].ravel();q/=np.linalg.norm(q)
    report=dict(created=tp.stamp(),status='running',source_path=str(srcfile),source_sha256=tp.driver.sha(srcfile),source=source,
        solver_sha256=tp.SHA,source_checks={},records=[],qualification='Finite local starting family; not global/all-spectrum envelope')
    save(report)
    try:
        model=tp.Pencil(HERE/'i11_matrices',native,11,p['beta_r'],p['omega_r']);report['source_checks']['polynomial']=model.build();save(report)
        r,q11,_=model.evaluate(p['alpha_r'],p['beta_r'],q,p['omega_r']+1e-3)
        deriv=[]
        for h in [1e-3,5e-4,1e-4]:
            rr=[]
            for da,db in [(-h,0),(h,0),(0,-h),(0,h)]:
                t,_,_=model.evaluate(p['alpha_r']+da,p['beta_r']+db,q11,p['omega_r']+1e-3,left=False);rr.append(t)
            wa=(z(rr[1],'omega')-z(rr[0],'omega'))/(2*h);wb=(z(rr[3],'omega')-z(rr[2],'omega'))/(2*h)
            deriv.append(dict(h=h,omega_alpha_FD=tp.cn(wa),omega_beta_FD=tp.cn(wb),
                alpha_relative_error=abs(wa-z(r,'omega_alpha'))/abs(z(r,'omega_alpha')),
                beta_relative_error=abs(wb-z(r,'omega_beta'))/max(abs(z(r,'omega_beta')),1e-30),samples=rr))
        report['source_checks'].update(point=r,temporal_FD=deriv,original_spatial_verification=source['source_neutral_record']['verification'])
        report['source_checks']['passed']=max(deriv[-1][k] for k in ['alpha_relative_error','beta_relative_error'])<1e-4 and abs(z(r,'omega')-p['omega_r'])<1e-8
        save(report);log('i11重建/导数 '+json.dumps({k:v for k,v in report['source_checks'].items() if k not in ['temporal_FD','original_spatial_verification']}))
        assert report['source_checks']['passed'],'derivative or source frequency mismatch'
        # Retain source finite-difference condition without pretending its coarse
        # h gives the exact analytic beta stationary point.
        source_refine,qref=tp.correct(model,p['alpha_r'],p['beta_r'],q11,p['omega_r']+1e-3,maxiter=20,backtracks=4)
        report['source_temporal_refinement']=source_refine;save(report)
        assert source_refine['status']=='temporal_candidate','known native point refinement failed'
        # Start target from actual old source; this is an initialization only.
        a0=source_refine['point']['alpha'];b0=source_refine['point']['beta'];w0=z(source_refine['point'],'omega').real
        target=tp.Pencil(HERE/'i10_matrices',native,10,b0,w0);checks=target.build()
        ans,q10=tp.correct(target,a0,b0,qref,w0+1e-3,maxiter=20,backtracks=4)
        report['records'].append(dict(station=10,matrix_checks=checks,result=ans,physical_candidate=False))
        save(report)
        if ans['status']=='temporal_candidate':
            t=ans['point'];point=dict(alpha_r=t['alpha'],alpha_i=0.,omega_r=t['omega']['real'],beta_r=t['beta'])
            actual=tp.vf.solve('M1c_i10_matrix_candidate_native',Path(native['input_dir']),400,10,point,here=HERE,binary=tp.BIN,sha=tp.SHA)
            qq=tp.vf.vectors(Path(actual['native_vector_file']),400)[actual['raw_index']].ravel()
            delta=abs(complex(actual['alpha_r'],actual['alpha_i'])-point['alpha_r']);overlap=tp.mac(q10,qq)
            report['records'][-1].update(native_fresh=actual,native_matrix_alpha_difference=delta,native_matrix_raw_MAC=overlap)
            report['records'][-1]['native_fresh_passed']=actual['passed'] and delta<1e-8 and overlap>1-1e-10 and abs(actual['growth'])<1e-7
            np.save(HERE/'i10_temporal_candidate_vector.npy',q10)
            save(report)
            assert report['records'][-1]['native_fresh_passed'],'fresh native branch/neutral check failed'
            report['status']='native_i10_neutral_and_MAC_passed_beta_curvature_pending'
        else:report['status']='native_i10_unresolved_bounded_temporal_corrector'
    except Exception:
        report['status']='bounded_stage_failed';report['exception']=traceback.format_exc();log(report['exception']);raise
    finally:
        tp.driver.snapshot();assert tp.driver.sha(tp.BIN)==tp.SHA and tp.driver.sha(srcfile)==report['source_sha256']
        report.update(completed=tp.stamp(),original_files_unchanged=True)
        save(report);log('阶段审查 '+report['status']+'；未通过独立空间负曲率核验前不能写candidate_converged。必要矩阵/配置/脚本/向量/失败均保留。')

if __name__=='__main__':main()
