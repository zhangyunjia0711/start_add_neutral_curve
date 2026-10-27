"""Create a readable report and archive only this task's verified artifacts."""
import csv,datetime,hashlib,json,shutil,subprocess,sys
from pathlib import Path
import numpy as np
WIN=Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate/75°')
BASE=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°')
DEST=WIN/'outputs/M3_neutral_eigenfunctions_20261027'
LOG=WIN/'outputs/工作日志_2026-10-27_中性边界特征函数频率对比.md'
REPO=BASE/'bluntplate/neutral_beta_full_20261022/work/solver_source'
TARGET=REPO/'doc/audit/20261027/M3_neutral_eigenfunctions'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def js(p,obj):p.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')
def main():
    evidence=json.loads((DEST/'eigenfunction_verification.json').read_text())
    native=json.loads((DEST/'tecplot_validation.json').read_text())
    for r in native['files']:
        assert sha(DEST/r['file'])==r['sha256'] and sha((DEST/r['file']).with_suffix('.dat'))==r['source_sha256']
    # Local code is the primary source for this implementation's actual conventions.
    refs=evidence['source_code']
    for rel in ['src/public/wall_local_coordinates.f90','src/lst/cfgio_adapter.f90','src/lst/global_parameters.f90']:
        path=REPO/rel;shutil.copy2(path,DEST/'source_evidence'/path.name)
        if not any(r['path']==str(path) for r in refs):refs.append(dict(path=str(path),sha256=sha(path)))
    for r in refs:
        path=Path(r['path']);assert sha(path)==r['sha256']
        target=DEST/'source_evidence'/path.name;shutil.copy2(path,target);assert sha(target)==r['sha256']
    js(DEST/'eigenfunction_verification.json',evidence)
    rows=list(csv.DictReader((DEST/'selected_mode_metrics.csv').open()))
    continuous=list(csv.DictReader((DEST/'frequency_evolution_metrics.csv').open()))
    maxres=max(r['matrix_scaled_residual'] for r in evidence['records'])
    maxai=max(abs(r['alpha_i']) for r in evidence['records'])
    text='''# M3 中性边界特征函数：五频率与连续频率变化

## 比较范围

沿用上一任务 x–β 叠加图中的同一个代表家族，f=−100、−50、0、+50、+100 kHz；每个频率取 β 较小和较大两侧中性根，共 20 个特征函数。不额外加入另一组共轭代表。

- bluntplate：原生截面 i=119，x=−25.0705680199 mm。
- rough：原生截面 i=24，x=−25.0270107051 mm。

每个工况内保持 x 相同；两工况 x 相差 0.0435573 mm，直接比较时应保留这个差别。本任务展示这些共同截面的频率变化，没有新增沿整条中性曲线的所有 x 截面特征函数扫描。

## 先看这些图

1. [w 特征函数两工况总对比](M3_eigenfunction_w_comparison.png)：近壁 0–5 mm，各工况、两侧边界分别比较五频率。
2. [连续变化热图](M3_eigenfunction_frequency_heatmap.png)：每工况每侧 139 个实际保存向量，共 556 点，展示 −100～100 kHz 的 |w|/S；没有用五个频率剖面插值代替中间的解。
3. [bluntplate 五分量幅值](bluntplate_eigenfunction_magnitudes.png) 和 [rough 五分量幅值](rough_eigenfunction_magnitudes.png)：ρ、u_t、v_n、w、T。
4. [形状变化指标](M3_eigenfunction_frequency_metrics.png)：w 原生网格峰位置、速度平方加权中心、相对零频的速度形状 MAC。峰位置图的台阶来自原生法向网格取点，不能解释为物理突跳。
5. [中性曲线采样点](M3_eigenfunction_sampling_points.png)：圆点为 β 较小侧，方点为 β 较大侧。

各图同时提供 SVG；`*_full_tail.*` 展示至 n≈48.81 mm 的外侧衰减尾部。`{case}_{lower|upper}_complex_eigenfunctions.*` 展示五个分量的实部和虚部，lower/upper 仅指 β 大小。

## 幅值、相位和坐标规范

实际 START 变量顺序是 q=[ρ,u_t,v_n,w,T]，每法向节点 5 个复数，400 个节点。依据 `source_evidence/module_lst.f90` 的输出标签与索引、`lst_solver_abstract.f90` 的重排方式以及 `lst_eqn.f90` 的系数定义核对；没有套用其他文献的变量顺序。u_t 沿局部壁面切向，v_n 沿局部壁面法向，w 沿展向。n 由原网格投影到冻结壁法向，长度基准 35 mm；坐标变换依据 `wall_local_coordinates.f90` 和各原运行的壁面坐标元数据。

整向量采用同一个实数 S=max_n sqrt(|u_t|²+|v_n|²+|w|²) 归一化，避免把五个分量分别缩放后抹去相对幅值变化。ρ、T 保留原求解器无量纲变量的基准，因此其纵坐标值不能直接与速度值比较物理能量，也不代表给定的真实扰动幅值。

全部 20 个向量再整体乘以一个单位复数，使同一物理位置 n=1.2071191626 mm 的 w 插值值为正实数。该位置由 20 个向量的最小 |w|/S 最大化选取；最弱点仍有 |w|/S≈0.72577。单一整体相位旋转保留分量间相位差，幅值图完全不受该相位选择影响。

上述处理依据特征方程的齐次性：若 P(α)q=0，则任意非零复常数 c 均有 P(α)(cq)=0。这是由方程直接得到的尺度自由度，而不是物理幅值结论。多项式特征问题形式见 [SLEPc 官方 PEP 手册，式(1)、(9)–(10)](https://slepc.upv.es/release/documentation/manual/pep.html)（2026-10-27 实际查阅）。相位和幅值方案为本任务选取的比较规范，未声称来自某一特定论文。

速度形状 MAC 的计算是 |∫v₀ᴴv_f dn|² / [(∫|v₀|² dn)(∫|v_f|² dn)]，其中 v=[u_t,v_n,w]。它不依赖整体复数尺度；不使用密度权重，不将其称作扰动能量。速度平方加权中心为 ∫n|v|² dn / ∫|v|² dn。

## 从本次数值中能直接看到的变化

- 两工况所取截面的形状和随频率变化的趋势接近；这里是局部观察，不能扩展成整个流场内粗糙效应很小。
- β 较小侧的 w 主峰主要在 n≈1.2–1.33 mm，形状变化较缓；u_t 的相对峰值随频率从 −100 到 +100 kHz 明显降低，五分量图可直接比较。
- β 较大侧的近壁形状、分量比例和外侧尾部变化更明显。−50 kHz 的外侧尾部比零频与正频率更长，速度平方加权中心也更外移；+50、+100 kHz 的 w 峰较零频更靠外。
- 连续热图和原追踪记录表明这些剖面沿已跟踪路径连续变化。相对零频的 MAC 较低不能单独证明模态切换；也不能由幅值或实/虚部曲线单独判定失稳机制。

这些观察同时伴随中性 β 和 α 的变化，不能把差别完全归因于频率。此处比较的是中性边界特征函数，不是固定 β 扫频，也不是失稳区内固定参数的增长模态。

## 五频率检查点

|工况|β 边界|f (kHz)|β|Re α|w 峰 n (mm)|速度平方中心 n (mm)|速度 MAC 相对 f=0|
|---|---|---:|---:|---:|---:|---:|---:|
'''
    for r in rows:
        text+='|'+r['case']+'|'+r['side']+'|'+r['f_kHz']+'|'+f"{float(r['beta']):.6f}|{float(r['alpha_r']):.6f}|{float(r['w_peak_n_mm']):.6f}|{float(r['velocity_squared_centroid_n_mm']):.6f}|{float(r['velocity_weighted_MAC_to_f0']):.6f}|\n"
    text+='\nβ 较大侧连续频率路径中，速度平方中心最大值（实际保存点，不作连续极值优化）：\n\n'
    for case in ['bluntplate','rough']:
        r=max([r for r in continuous if r['case']==case and r['side']=='upper'],key=lambda r:float(r['velocity_squared_centroid_n_mm']))
        text+=f"- {case}: f={float(r['f_kHz']):.4f} kHz，中心 n={float(r['velocity_squared_centroid_n_mm']):.6f} mm。\n"
    text+=f'''
## 数值核验与证据局限

- 20 个原始向量逐一核对 SHA256，与上一任务中性曲线的 β 和 α 相差均小于 1e−7。
- 原空间矩阵 P(α)=K₀+αK₁+α²K₂ 的本次复核残差最大 {maxres:.4e}。定义为 ||P(α)q||₂ / [(||K₀||F+|α|||K₁||F+|α|²||K₂||F)||q||₂]；这是系数矩阵 Frobenius 尺度化残差，不能与其他残差定义直接混比。
- 最大 |Im α|={maxai:.4e}，符合所取中性根的误差水平。
- 复用先前同 20 点的冻结 START 独立求解收据，原向量 MAC 最小 0.9999999999999987；本任务重新检查矩阵和向量，没有再进行这 20 次原求解。连续图 556 个原向量全部核对哈希，本次新增矩阵残差核验仅覆盖五频率的 20 个点。
- 壁面 u_t,v_n,w,T 最大归一化值 5.83e−20，远场相应值 1.13e−37。原配置为无滑移、默认等温；默认设置核对 `cfgio_adapter.f90`、`global_parameters.f90`。
- 检查任意复数缩放后重新执行同样规范，归一化复特征函数差异小于 1e−11；幅值不受整体相位旋转影响。
- 三个二进制 .plt 经 Preplot 生成，再由 Tecplot 360 EX 2024 R1 读入和输出；逐个双精度值完全一致。合并文件 20 zones、8000 点、184000 个数值，单工况各 10 zones、4000 点、92000 个数值。
- 模型沿用已核验的局部壁面坐标、冻结切平面、零显式曲率的局部 LST，不是新的全局/非平行鼻部计算。原始输入保护快照保持一致。
- 本次重取 Nomura (2006) JStage PDF 超时，没有借该失败检索支撑新增理论结论；实际变量与处理规范均依据本地原代码和成功读取的 SLEPc 一手文档。

## 数据文件与 Tecplot 使用

- `M3_neutral_eigenfunctions.plt`：20 个 zone 合并；`bluntplate_M3_neutral_eigenfunctions.plt`、`rough_M3_neutral_eigenfunctions.plt` 各含 10 个 zone。
- X 选 `n_mm`；Y 可选 `rho_abs,u_abs,v_abs,w_abs,T_abs` 或各自 `real,imag`。`u` 表示 u_t，`v` 表示 v_n。case_id=1/2 对应 bluntplate/rough；boundary_id=1/2 对应 lower/upper。f_kHz、β、α、x 都保存在每个点中。
- `.lay` 默认展示 w 幅值，频率用颜色区分；实线/虚线分别是 lower/upper β 中性根，与额外共轭分支无关。
- `*_f±*kHz.csv` 保存每个剖面的 23 个双精度变量，`.dat` 是 `.plt` 对应的完整 ASCII 输入。
- `eigenfunctions_original_and_normalized.npz` 保留 20 个原始和规范化的 400×5 复向量、两工况法向网格和基流；原始向量保持原尺度。
- `continuous_w_amplitude_maps.npz` 保留连续热图在原法向网格上的数值；`continuous_vector_sources.json` 列出 556 个原始向量路径和哈希。
- `selected_mode_metrics.csv` 与 `frequency_evolution_metrics.csv` 分别保存 20 点和 556 点的形状指标。
- `eigenfunction_verification.json`、`tecplot_validation.json`、`source_evidence/` 保存检查证据；工作日志记录失败尝试、绘图修订及清理。

## 原环境复现

在 WSL Ubuntu 中运行 `eigenfunctions.py`（依赖原冻结求解器、矩阵、追踪向量和前一任务叠加图的原路径，以及 numpy/scipy/h5py/matplotlib），可以重新读取原数据、复核 20 个向量并绘图。该脚本只作后处理，不启动新的特征值计算。Windows 执行 `native_tecplot.ps1`，再在 WSL 运行 `eigenfunctions.py --verify-native` 检查 Tecplot 回读；读回验证需先重新生成回读 .dat。`finalize.py` 创建说明与版本归档，归档脚本会拒绝混入其他任务修改。

归档中已保留图、数值和原始向量，无需重新求解即可查看。本轮没有新改动 Fortran 求解模块。
'''
    (DEST/'README_特征函数对比.md').write_text(text,encoding='utf-8')
    shutil.copy2(__file__,DEST/'finalize.py')
    # Archive a traceable copy of the source continuation report without altering it.
    source=Path(evidence['source_report']);assert sha(source)==evidence['source_report_sha256']
    shutil.copy2(source,DEST/'source_connection_report.json')
    shutil.copy2(LOG,DEST/LOG.name)
    # Keep source DAT and reproducibility artifacts; only remove verified roundtrip copies.
    cleanup=[]
    for r in native['files']:
        p=DEST/(Path(r['file']).stem+'_roundtrip.dat')
        if p.exists():
            assert p.resolve().parent==DEST.resolve() and sha(p)==r['roundtrip_sha256']
            cleanup.append(dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p),reason='Verified exact native Tecplot read-back; input DAT and validation receipt retained'))
            p.unlink()
    old=DEST/'cleanup_receipt.json'
    if old.exists():cleanup=json.loads(old.read_text())['removed']+cleanup
    js(old,dict(removed=cleanup))
    sys.path.insert(0,str(BASE/'neutral_xbeta_four_omega_20261026/work'))
    from common import snapshot
    assert snapshot()==evidence['protected_inputs_before']
    prefix=str(TARGET.relative_to(REPO))+'/'
    for command in [['diff','--cached','--name-only','-z'],['diff','--name-only','-z'],['ls-files','--others','--exclude-standard','-z']]:
        names=subprocess.check_output(['git','-C',str(REPO)]+command).decode().split('\0')
        assert all(p.startswith(prefix) for p in names if p),'Do not mix another task'
    for p in DEST.rglob('*'):
        if p.is_file() and p.suffix in ['.svg','.csv','.json','.lay','.md','.py','.ps1','.mcr']:
            s=p.read_text(encoding='utf-8-sig');p.write_text('\n'.join(v.rstrip() for v in s.splitlines()).rstrip()+'\n',encoding='utf-8')
    (DEST/'.gitattributes').write_text('*.plt binary\n*.log binary\nsource_evidence/*.f90 binary\n',encoding='ascii')
    files=[dict(path=str(p.relative_to(DEST)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(DEST.rglob('*')) if p.is_file() and p.name not in ['delivery_manifest.json','git_version.json'] and not p.name.startswith('工作日志_')]
    js(DEST/'delivery_manifest.json',dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),files=files,protected_inputs=snapshot(),execution_metadata_excluded=['工作日志_*.md','git_version.json','delivery_manifest.json']))
    shutil.copytree(DEST,TARGET,dirs_exist_ok=True)
    for p in TARGET.rglob('*'):
        if p.is_file():p.chmod(0o644)
    assert all(sha(TARGET/r['path'])==r['sha256'] for r in files)
    print('ARCHIVE_READY',TARGET,'files',len(files),'bytes',sum(r['bytes'] for r in files),flush=True)

if __name__=='__main__':main()
