"""Export the verified native stationary minimum separately from beta maxima."""
import csv,json
from pathlib import Path
import temporal_pencil as tp
RUN=tp.RUN;HERE=Path(__file__).resolve().parent
LOG=RUN/'outputs/工作日志_2026-10-23_矩阵频率换元M1c原生鼻部延拓.md'
def main():
    files=[RUN/'outputs'/name for name in ['matrix_temporal_M1c_native_20261023.json','matrix_temporal_M1c_native10_secant_20261023.json','matrix_temporal_M1c_native10_curvature_20261023.json','matrix_temporal_M1c_profile_type_change_20261023.json','matrix_temporal_M1c_profile_type_diagnostic_20261023.json']]
    initial,secant,curv,flip,diag=[json.loads(p.read_text()) for p in files]
    actual=secant['native_fresh'];t=secant['result']['point'];a=complex(t['omega_alpha']['real'],t['omega_alpha']['imag']);gomega=-float((1/a).imag)
    assert secant['native_fresh_passed'] and curv['status']=='native_beta_local_minimum_verified_not_envelope_peak'
    p=dict(outer=t['omega']['real'],inner=t['beta'],iloc=10.,x=float(actual['metadata']['x_global_N']),Re=79100.,
        alpha_r=actual['alpha_r'],alpha_i=actual['alpha_i'],beta_r=t['beta'],beta_i=0.,omega_r=t['omega']['real'],omega_i=0.,growth=actual['growth'],branch=0.,status=0.)
    sourcehashfile=RUN/'outputs/temporal_frozen_global_source_hashes_20261023.json';sourcehashes=json.loads(sourcehashfile.read_text())
    assert all(tp.driver.sha(Path(x['path']))==x['sha256'] for x in sourcehashes['files'])
    tp.driver.snapshot();assert tp.driver.sha(tp.BIN)==tp.SHA
    report=dict(created=tp.stamp(),status='stationary_neutral_beta_minimum_verified',station=10,source_family='M1c',
        coordinates='wall_local',solver=str(tp.BIN),solver_sha256=tp.SHA,physical_grid_station=True,
        point=p,actual_native=actual,temporal_point=t,raw_temporal_native_MAC=secant['native_matrix_raw_MAC'],
        alpha_temporal_native_difference=secant['native_matrix_alpha_difference'],g_omega=gomega,
        g_beta_matrix=t['cross_divided_abs_omega_alpha_squared'],spatial_beta_curvature=curv,
        neutral_stationary_candidate=True,beta_peak_candidate=False,upper_envelope_candidate=False,
        mathematical_role='Local beta minimum on a finite source family, with positive fixed-beta g_omega; not a max-beta envelope point',
        provenance=[dict(path=str(p),sha256=tp.driver.sha(p)) for p in files],
        core_hash_check=dict(source_manifest=str(sourcehashfile),source_manifest_sha256=tp.driver.sha(sourcehashfile),files_unchanged=True),
        helper_files=[dict(path=str(p),sha256=tp.driver.sha(p)) for p in HERE.glob('*.py')],
        unresolved=dict(native_i1_to_i9='not extended from beta minimum',native_i10_max_beta_envelope='unknown; no beta maximum verified by this continuation',
            exact_profile_degeneracy_location='not localized; auxiliary near-stationary diagnostics are not formal converged roots'),
        references=['https://www.cs.cornell.edu/courses/cs6210/2025fa/lec/2025-10-22.html','https://slepc.upv.es/release/src/pep/interface/pepsolve.c.html','https://docs.scipy.org/doc/scipy-1.10.1/reference/generated/scipy.sparse.linalg.eigs.html'],
        qualification='Actual finite-family stationary neutral candidate under the existing local parallel wall-coordinate spatial operator, not full-spectrum/stability proof or a claim of physical absolute instability')
    out=RUN/'outputs/native10_stationary_beta_minimum_20261023.json';assert not out.exists();out.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    table=out.with_suffix('.csv')
    with table.open('w') as f:
        row=dict(station=10,source='M1c',classification='beta_local_minimum',x_N=p['x'],omega=p['omega_r'],beta=p['beta_r'],alpha_r=p['alpha_r'],alpha_i=p['alpha_i'],growth=p['growth'],g_omega=gomega,g_betabeta=curv['chain_checks'][-1]['g_betabeta'],upper_envelope_candidate=False)
        w=csv.DictWriter(f,fieldnames=list(row));w.writeheader();w.writerow(row)
    text=['# M1c 矩阵频率换元与原生 i10 驻点审查',tp.stamp(),
        '本阶段在原 v10 wall_local 空间 PEP 导出的同一系数矩阵上换元，没有修改 core、原 H5 或已有求解门限。独立 beta=beta0+0.37 的三系数矩阵核验相对误差最高7.75e−17，G 的 beta 变化范数1.82e−14。M1c i11 左右向量回代残差约3e−13；解析 omega_alpha 与三个有限步长 temporal FD 的相对误差低于9e−7，omega_beta低于2e−9。',
        '从原生 i11 与 i12 的实际物理x坐标构造割线初值 alpha=−22.217562401419855、beta=12.16214490428434、omega=8.920739732795656，独立重策11迭代收敛。此前直接沿用i11参数的18迭代no-descent失败完整保留；没有把失败解释为无根。',
        '|原生站|omega|beta|fresh alpha|g|g_omega|空间beta曲率|',
        '|---|---:|---:|---:|---:|---:|---:|',
        '|10|%.15g|%.15g|%.15g %+.4gi|%.5g|%.9g|+1.70155e7|'%(p['omega_r'],p['beta_r'],p['alpha_r'],p['alpha_i'],p['growth'],gomega),
        'fresh原H5求解与矩阵alpha之差1.08e−9，原始向量MAC约1；PEP实际残差/向量文件/索引/hash见独立JSON。三个原生beta三点采样步长1e−6、5e−7、2.5e−7的曲率为+1.70294e7、+1.70172e7、+1.70373e7，与矩阵链式h/h2曲率一致。因此这是实际beta局部最小值，不能加入对beta最大增长包络。其g_omega为正，不能沿用i11 upper的频率边界标签。',
        '简单解析分支下，我们由Kq=omega Mq及其左向量推得omega_t=p^H K_t q/(p^H Mq)（本次M独立参数），再由隐函数推得alpha_beta=−omega_beta/omega_alpha。因此实omega的空间驻点条件是Im(omega_beta*conj(omega_alpha))=0，但必须排除omega_alpha不可解析/接近零的退化；本次收敛标准另检查cross/|omega_alpha|²与Imomega/|omega_alpha|，避免未缩放分子很小导致假通过。一般非厄米简单根左右向量一阶公式已实际查阅[Bindel/Cornell 2025课堂原文](https://www.cs.cornell.edu/courses/cs6210/2025fa/lec/2025-10-22.html)；广义pencil、隐函数及曲率链式公式为我们针对此算子的代数推导，并以实际矩阵和原生差分复核。这里没有使用Gaster小增长近似。',
        '为检查类型翻转，额外只取既有投影profile辅助列theta=.625、.75。两列的严格1e−7驻点门限均未可靠通过，F2约1e−6噪声导致no-descent；原unknown和失败保留，未放宽门限或继续Newton。仅对已保存近驻点做h/h2诊断：theta=.625的g_omega≈−644.66、曲率≈−6.452e9；theta=.75的g_omega≈+365.54、曲率≈+1.206e9。它们支持数值上发生类型转换的待核验迹象，但不是正式中性根、物理CFD站或精确omega_alpha零点定位，也未证明谱重根。',
        '交付native10_stationary_beta_minimum_20261023.json/csv单独记录该真实驻点，未写入upper或beta-max表；原生i1…9未据此最小值继续，i10最大beta中性包络仍未知。所有矩阵只覆盖当前来源族，不是全谱搜索或完整联合稳定性判据。鼻部local parallel/zero-curvature模型的物理适用性限制继续有效。',
        '原始输入/bin和只读评估的9个core源文件hash均保持。所有脚本、矩阵、向量、配置、失败trace和日志是必要复现证据，保留；本阶段没有确认不再需要的尝试文件可清理。']
    # Join contiguous table rows without inserting blank rows.
    blocks=[];rows=[]
    for line in text:
        if line.startswith('|'):rows.append(line);continue
        if rows:blocks.append('\n'.join(rows));rows=[]
        blocks.append(line)
    if rows:blocks.append('\n'.join(rows))
    out.with_name('matrix_temporal_M1c_native_20261023_说明.md').write_text('\n\n'.join(blocks)+'\n')
    with LOG.open('a') as f:f.write('\n'+tp.stamp()+' 结束审查：i11矩阵重建/左右导数完成，native10一次x-secant重策得到真中性beta最小值，gomega='+str(gomega)+'；独立JSON/CSV不并upper。有限辅助profile未达严格驻点门限，保留unknown。原H5/bin/core hash保持；必要复现文件保留，无额外清理。i10 beta-max envelope及i1..9仍未解决，后续有限beta cut由另一agent负责。\n')
    print(str(out))
if __name__=='__main__':main()
