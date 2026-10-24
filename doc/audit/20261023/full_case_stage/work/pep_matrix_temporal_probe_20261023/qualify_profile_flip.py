"""Curvature at saved near-stationary auxiliary iterates, without more Newton."""
import json
from pathlib import Path
import numpy as np
import temporal_pencil as tp
RUN=tp.RUN;HERE=Path(__file__).resolve().parent/'M1c_profile_flip'
OUT=RUN/'outputs/matrix_temporal_M1c_profile_type_diagnostic_20261023.json'
LOG=RUN/'outputs/工作日志_2026-10-23_矩阵频率换元M1c原生鼻部延拓.md'
def main():
    assert not OUT.exists()
    srcpath=RUN/'outputs/matrix_temporal_M1c_profile_type_change_20261023.json';src=json.loads(srcpath.read_text())
    manifest=src['auxiliary_input_manifest'];native=manifest['source_anchor_actual'];q=tp.vf.vectors(Path(native['native_vector_file']),400)[native['raw_index']].ravel();q/=np.linalg.norm(q)
    result=dict(created=tp.stamp(),source=str(srcpath),source_sha256=tp.driver.sha(srcpath),physical_candidate=False,formal_auxiliary_root_converged=False,records=[])
    for old in src['records']:
        best=min(old['result']['trace'],key=lambda v:np.linalg.norm(v['F']));point=best['center'];wa=complex(point['omega_alpha']['real'],point['omega_alpha']['imag']);ab=point['implicit_alpha_beta']['real']
        source=dict(native);source['config_path']=str(HERE/'aux_template.ini')
        pred=old['prediction'];model=tp.Pencil(HERE/('aux%d_matrices'%old['auxiliary_column']),source,old['auxiliary_column'],pred['beta_r'],pred['omega_r']);model.build()
        rec=dict(auxiliary_column=old['auxiliary_column'],theta=old['theta'],selected_saved_iteration=best['iteration'],actual_iterate=point,
            original_corrector_status=old['result']['status'],F=best['F'],g_omega=-float((1/wa).imag),curvature_checks=[])
        for h in [1e-4,5e-5]:
            values=[]
            for da,db in [(-h,0),(h,0),(0,-h),(0,h)]:
                rr,_,_=model.evaluate(point['alpha']+da,point['beta']+db,q,point['omega']['real']+1e-3);values.append(rr)
            f=[v['cross_divided_abs_omega_alpha_squared'] for v in values];Ja=(f[1]-f[0])/(2*h);Jb=(f[3]-f[2])/(2*h);curv=Jb+ab*Ja
            rec['curvature_checks'].append(dict(h=h,approximate_g_betabeta_at_saved_near_stationary_iterate=curv,
                implied_beta_location_correction=-best['F'][1]/curv,samples=values))
        result['records'].append(rec)
    tp.driver.snapshot();assert tp.driver.sha(tp.BIN)==tp.SHA and tp.driver.sha(srcpath)==result['source_sha256']
    result.update(status='finite_near_stationary_type_diagnostic_not_formal_roots',completed=tp.stamp(),original_files_unchanged=True,
        qualification='The 1e-7 exact implicit-beta-gradient threshold was not reached reliably: ~1e-6 derivative noise at tiny omega_alpha. No threshold relaxed, no more Newton, no physical station accepted. Curvature signs only characterize saved near-stationary iterates and do not prove exact degeneracy location or root multiplicity.')
    OUT.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    with LOG.open('a') as f:f.write('\n'+tp.stamp()+' 两辅助列20iter/4回溯内均no-descent；F2约1e−6而严格门限1e−7，保留unknown。未追加Newton/放宽门限，仅原已算近驻点的h/h2曲率诊断 '+str([(v['theta'],v['g_omega'],[(x['h'],x['approximate_g_betabeta_at_saved_near_stationary_iterate']) for x in v['curvature_checks']]) for v in result['records']])+'；不是正式root或物理x点。\n')
    print(json.dumps({k:v for k,v in result.items() if k!='records'}));print([(v['theta'],v['g_omega'],[(x['h'],x['approximate_g_betabeta_at_saved_near_stationary_iterate']) for x in v['curvature_checks']]) for v in result['records']])
if __name__=='__main__':main()
