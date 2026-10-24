"""Two bounded auxiliary-profile checks of stationary-root type change.

These columns interpolate existing projected profiles, not physical CFD stations.
"""
import configparser,json,traceback
from pathlib import Path
import numpy as np
import temporal_pencil as tp
RUN=tp.RUN;HERE=Path(__file__).resolve().parent/'M1c_profile_flip';HERE.mkdir(exist_ok=True)
OUT=RUN/'outputs/matrix_temporal_M1c_profile_type_change_20261023.json'
LOG=RUN/'outputs/工作日志_2026-10-23_矩阵频率换元M1c原生鼻部延拓.md'
def save(r):OUT.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
def point(record):
    p=record['result']['point'];return p
def main():
    assert not OUT.exists()
    manifestpath=RUN/'outputs/nose_M1c_profile_homotopy_inputs_11_10_steps64_20261023.json';manifest=json.loads(manifestpath.read_text())
    oldpath=RUN/'outputs/aux_neutral_M1c_11_10_s64_20261023.json';old=json.loads(oldpath.read_text())
    endpointpath=RUN/'outputs/matrix_temporal_M1c_native10_secant_20261023.json';endpoint=json.loads(endpointpath.read_text())['result']['point']
    # JSON schema is inspected rather than treating a failed iteration as seed.
    accepted=[x for x in old['records'] if x['result']['status']=='candidate_converged']
    p40=accepted[-1]['result']['point'];p39=accepted[-2]['result']['point']
    native=manifest['source_anchor_actual'];source=dict(native)
    cfg=configparser.ConfigParser();cfg.read(native['config_path'])
    cfg['filenames']['dir']=manifest['input_dir']+'/'
    cfg['domain']['in']='65';cfg['domain']['istart']='1';cfg['domain']['iend']='65'
    template=HERE/'aux_template.ini'
    with template.open('w') as f:cfg.write(f)
    source['config_path']=str(template)
    q=tp.vf.vectors(Path(native['native_vector_file']),400)[native['raw_index']].ravel();q/=np.linalg.norm(q)
    report=dict(created=tp.stamp(),status='running',physical_candidate=False,
        provenance=[dict(path=str(p),sha256=tp.driver.sha(p)) for p in [manifestpath,oldpath,endpointpath]],
        auxiliary_input_manifest=manifest,records=[],qualification='Two interpolated local-profile columns, not physical CFD x stations and not exact degeneracy localization')
    save(report)
    for column,theta in [(41,.625),(49,.75)]:
        if column==41:
            pred={k:p40[k]+(p40[k]-p39[k]) for k in ['alpha_r','beta_r','omega_r']}
        else:
            f=(theta-.609375)/(1-.609375)
            pred=dict(alpha_r=p40['alpha_r']+f*(endpoint['alpha']-p40['alpha_r']),
                beta_r=p40['beta_r']+f*(endpoint['beta']-p40['beta_r']),omega_r=p40['omega_r']+f*(endpoint['omega']['real']-p40['omega_r']))
        model=tp.Pencil(HERE/('aux%d_matrices'%column),source,column,pred['beta_r'],pred['omega_r'])
        rec=dict(auxiliary_column=column,theta=theta,physical_station=False,prediction=pred)
        report['records'].append(rec);save(report)
        try:
            rec['matrix_checks']=model.build()
            ans,qq=tp.correct(model,pred['alpha_r'],pred['beta_r'],q,pred['omega_r']+1e-3,maxiter=20,backtracks=4)
            rec['result']=ans
            if ans['status']=='temporal_candidate':
                p=ans['point'];wa=complex(p['omega_alpha']['real'],p['omega_alpha']['imag']);ab=p['implicit_alpha_beta']['real']
                rec['g_omega']=-float((1/wa).imag);rec['curvature_checks']=[]
                for h in [1e-4,5e-5]:
                    vv=[]
                    for da,db in [(-h,0),(h,0),(0,-h),(0,h)]:
                        v,_,_=model.evaluate(p['alpha']+da,p['beta']+db,qq,p['omega']['real']+1e-3);vv.append(v)
                    f=[v['cross_divided_abs_omega_alpha_squared'] for v in vv];ca=(f[1]-f[0])/(2*h);cb=(f[3]-f[2])/(2*h)
                    rec['curvature_checks'].append(dict(h=h,g_betabeta=cb+ab*ca,partial_alpha=ca,partial_beta=cb,samples=vv))
                np.save(HERE/('aux%d_stationary_vector.npy'%column),qq)
                signs=[v['g_betabeta'] for v in rec['curvature_checks']]
                rec['classification']='beta_maximum' if max(signs)<0 else ('beta_minimum' if min(signs)>0 else 'unresolved_curvature')
            else:rec['classification']='unknown_bounded_corrector'
        except Exception:rec['exception']=traceback.format_exc();rec['classification']='unknown_actual_exception'
        rec['temporal_calls']=model.calls;save(report)
        print('AUX TYPE '+json.dumps(dict(column=column,theta=theta,classification=rec['classification'],g_omega=rec.get('g_omega'),point=rec.get('result',{}).get('point'))),flush=True)
    tp.driver.snapshot();assert tp.driver.sha(tp.BIN)==tp.SHA and all(tp.driver.sha(Path(p['path']))==p['sha256'] for p in report['provenance'])
    report.update(status='two_auxiliary_checks_complete',completed=tp.stamp(),original_files_unchanged=True);save(report)
    with LOG.open('a') as f:f.write('\n'+tp.stamp()+' 有限辅助profile类型审查 '+json.dumps([(v['theta'],v['classification'],v.get('g_omega')) for v in report['records']])+'。辅助列不作物理x站；没有精确定位omega_alpha零点或证明谱重根。必要矩阵/脚本/记录保留。\n')

if __name__=='__main__':main()
