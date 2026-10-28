"""Package verified figures, provenance and a bounded reproducibility archive."""
import csv,datetime,hashlib,json,os,shutil,tarfile
from pathlib import Path
BASE=Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate')
ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/neutral_m3_sweep_consistent_20261028')
WORK=BASE/'work/M3_sweep_angles_20261027'
DEST=BASE/'outputs/M3_sweep_angles_20261028'
NEW=[30,45,50,55,60,65,70,73]
STAMP=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def js(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n')

def main():
    selection=json.loads((DEST/'selection_manifest.json').read_text());curves=selection['curves']
    audit=json.loads((ROOT/'outputs/final_numerical_audit.json').read_text())
    assert len(curves)==164 and len(audit['curves'])==141 and audit['passed']
    assert all(sha(ROOT/'outputs'/r['curve'])==r['sha256'] for r in audit['curves'])
    coverage=[]
    for angle in NEW+[75]:
        records=[r for r in curves if r['angle']==angle]
        coverage.append(dict(angle=angle,status='retained_original_no_recomputation' if angle==75 else 'computed_from_1p1m_FPP',curves=len(records),frequencies=len(set(r['f_kHz'] for r in records)),points=sum(r['points'] for r in records),sampled_min_f_kHz=min(r['f_kHz'] for r in records),sampled_max_f_kHz=max(r['f_kHz'] for r in records),actual_domain_end_mm=3500 if angle==70 else 1100,note='Frequency extrema are delivered samples, not physical stability cutoffs.'))
    for angle in [0,10,15,63,68,74,74.2,74.5,78,80]:
        coverage.append(dict(angle=angle,status='coarse_FPP_available_excluded_from_dense_comparison' if angle==74 else 'no_complete_selected_FPP_available',curves=0,frequencies=0,points=0,sampled_min_f_kHz=None,sampled_max_f_kHz=None,actual_domain_end_mm=None,note='74°/rough/fpp_output is coarse and is listed separately; not missing.' if angle==74 else 'Original files or derived eigenvalue PLT may exist; no replacement of missing selected FPP by coarse results.'))
    js(DEST/'case_coverage.json',dict(created=STAMP,entries=coverage,inventory='all_fpp_final_inventory.json',classification_rule='All data in 1.1m are dense by user definition; rough means coarse.'))
    with (DEST/'case_coverage.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(coverage[0]));w.writeheader();w.writerows(coverage)
    # Refresh the older pending-input description so it is not mistaken for final coverage.
    pending=json.loads((DEST/'pending_input_status.json').read_text());pending['superseded_by']='case_coverage.json';pending['historical_inventory_only']=True;js(DEST/'pending_input_status.json',pending)
    table='\n'.join('| {angle:g}° | {curves} | {frequencies} | {points} | {sampled_min_f_kHz:g}～{sampled_max_f_kHz:g} |'.format(**r) for r in coverage if r['curves'])
    readme=f'''# 不同后掠角 M3 x–β 中性曲线

本轮交付164条曲线、85124个数据点。新计算30°、45°、50°、55°、60°、65°、70°、73°的141条曲线；75°的23条曲线逐文件SHA核验后原样复用，没有重新求解75°。

## 查看图件和Tecplot数据

- `M3_sweep_angles_comparison.png`：所有角度的共同1.1m窗口与上游细节；SVG同名。
- `M3_sweep_angles_upstream.png`：九角度上游细节，适合先查看。
- `M3_sweep_angles_common_1p1m.png`：九角度统一x窗口；纵轴按各角度范围显示。
- `M3_sweep_angles_shared_beta_scale.png`：x和β均用共同范围，适合比较尺度。
- `M3_<角度>deg_xbeta.png`：各角度完整实际计算范围及上游细节。
- `M3_sweep_angles_all_neutral_curves.plt`：原生Tecplot二进制，164个独立zone、15变量；各角度另有独立 `.plt`。相应 `.dat` 是ASCII版本。
- `curves/`：逐条CSV；75°CSV为原始13列原样副本，新算CSV为15列。

Tecplot使用XY Line，以`x_mm`为横坐标、`beta`为纵坐标。zone名称包含角度、带符号频率、component。相同频率使用相同颜色；负频率虚线、正频率实线、零频率黑线。不同component保留独立zone，没有添入共轭复制曲线。总图图例为各角度频率的并集，个别角度没有交付某频率时，不代表该频率已证实稳定。

## 输入选择和实际覆盖

按用户规则，`1.1m`内均认定为密网格，不按网格文件名或点数改判。新八角度均来自其`1.1m/fpp_output`；`rough`只表示粗网格。保留FPP原400点法向分布，仅作转置/单位换算。

70°的`1.1m/fpp_output`实际含151个沿程截面、x=−35～3500mm，其余新角度含841个截面、x=−35～1100mm。共同图只显示至1100mm，70°独立图、PLT、CSV保留3500mm范围。

| 后掠角 | 曲线数 | 频率数 | 点数 | 已交付采样频率范围/kHz |
|---|---:|---:|---:|---:|
{table}

以上最小/最大频率是已交付采样范围，**不是失稳截止频率**。`frequency_status.csv`逐角度/频率列出状态；`unresolved_no_delivered_curve`表示本轮未求得可交付曲线，不能解释为全域稳定。固定x搜索和允许x移动的频率延拓均为有限分支搜索，不能证明所有独立口袋已经穷尽。

63°、68°、80°当前仍未发现完整的所选密网格FPP输入；不从其单β正增长率PLT伪造中性边界。0°、10°、15°、74.2°、74.5°、78°也未发现完整的所选FPP数据；74°现有`rough/fpp_output`为粗网格，本轮密网格对比未混入，已单独列入`case_coverage.json`。完整目录复查见`all_fpp_final_inventory.json`，其中其他头部半径/构型属于另外的参数变化。

## 模型和比较限制

横坐标为全局x/mm；β、α、ω均以Lref=0.035m、Uref=1810.6307411853716m/s无量纲化。频率标签为带符号的循环物理频率f/kHz，代码换算为`omega_nd = 2*pi*1000*f_kHz*Lref/Uref`，不是把kHz直接作为角频率。密度和温度参考量见`frequency_conversion.json`与`input_selection.json`。

新计算使用冻结局部壁面基底、零显式曲率的空间LST模型，选择`wall_local_consistent`使所有截面统一采用投影法向导数。源码实现与配置保存在复现代码及证据包中。75°按用户要求保留原输入处理和旧法向离散；因此这九角度图不是只改变后掠角的完全受控定量比较。不能将新旧处理差异全部归因于后掠角，也不据此声称CFD网格收敛或已包含全部曲率/非平行效应。

f=0在各自下游计算端仍开放；55°还有非零频率开放边界。开放线未人工闭合。原生站之间的x转折点由相邻实际局部原始变量和法向坐标构建后重新求解；这不是增加CFD网格。70°只将经独立中点误差检查的自适应三次矩阵用于导引，最终原生点与转折点仍独立校正。

## 验证与可追溯记录

- 新141条曲线：74221点，1544次独立START原始输入抽查；最大记录增长率误差约2.206e−10，最大记录右残差约2.464e−11。
- 698个选点以原生或实际剖面矩阵直接复核；原生站上下边界成对且无漏站，已闭合曲线闭合校正通过。见`final_numerical_audit.json`。
- 新法向离散回归与旧模式兼容、坐标舍入扰动、全场/9站导出对照通过，见`consistent_validation.json`、`build_receipt.json`。CMake选中的4项测试通过。
- 二进制PLT由本机Tecplot PREPLOT 2024.1导出，再由独立受限读取器逐值核对：总PLT的1276860个double均与CSV一致，各角度PLT也一致。见`tecplot_validation.json`。读取器依据本机官方`360-data-format.html`的Binary Data File Format节，不修改许可。
- 75°的23份CSV与原交付SHA逐项一致，见`plot_validation.json`；74项受保护输入当前哈希均通过，原始CFD/FPP不改写。
- 70°11个种子的旧频率搜索审计被重启覆盖，索引已明确标为失效历史引用；不补造遗失的搜索历史。其实际种子与最终曲线均有独立START校验。见`provenance_corrections.json`。

以上数值为本轮实际检查记录，不是关于所有特征模式或频率截止值的理论结论。本轮没有新增文献性物理推论。

## 复现文件

`reproducibility/runtime/`为最终运行脚本，`preparation_history/`保存准备与诊断流程，`frozen_solver/START_LST`及SHA保存实际计算二进制。代码备份包含新增可选模式与回归测试；旧global/wall_local行为保留。历史`build_consistent.py`针对修改前基准构建，不可将它当作可对当前HEAD重复应用的补丁。

本机完整计算阶段仍在：
`{ROOT}`

完成结果不会被恢复命令覆盖。需要检查/继续某个新角度时，在上述目录执行：
```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 work/run_batch.py s30 verify_resume
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 work/audit_final_results.py
```
`s30`可替换为其余七个新角度；75°数值入口被明确禁止。此恢复命令对已完成曲线仅检查并跳过，不代表重新独立计算了全部数据。

轻量证据包`reproducibility/numerical_evidence.tar.gz`包含最终/导引JSON、搜索记录、输入哈希、构建测试记录、独立START选点的日志/表格、选定特征向量及70°导引误差回执。为控制体积，未打包原始CFD/FPP、全部矩阵缓存或全部导引向量；全量本地阶段保留。跨机器复算需另行提供清单对应的原始输入或已转换HDF5，并调整脚本中绝对路径与PETSc/SLEPc/MPI环境，不能仅凭轻量包从零重建所有矩阵。特征向量抽样包用于追溯，不是全量断点包。

`evidence_archive_manifest.json`记录压缩包成员及SHA；`delivery_manifest.json`记录最终交付文件SHA。版本提交和远程核验另记在`git_version.json`及工作日志。
'''
    (DEST/'README.md').write_text(readme)
    repro=DEST/'reproducibility';(repro/'runtime').mkdir(parents=True,exist_ok=True)
    for p in (ROOT/'work').glob('*.py'):shutil.copy2(p,repro/'runtime'/p.name)
    for name in ['publish_sweep.py','export_verify_tecplot.py','read_plt.py','stage_delivery.py']:
        shutil.copy2(WORK/name,repro/'runtime'/name)
    (repro/'preparation_history').mkdir(exist_ok=True)
    for name in ['build_consistent.py','verify_consistent.py','prepare_consistent.py','prepare_fpp.py','add_fpp_angles.py','pilot_additional.py','physical_mode_transfer.py','physical_mode_transfer_aligned.py','probe_mixed_flat.py','recover_consistent.py']:
        shutil.copy2(WORK/name,repro/'preparation_history'/name)
    (repro/'frozen_solver').mkdir(exist_ok=True);binary=ROOT/'work/solver_source/bin/START_LST';shutil.copy2(binary,repro/'frozen_solver/START_LST')
    js(repro/'version_provenance.json',dict(created=STAMP,base_commit='b391bde3184540a820e554b153082370e6ef3819',solver_sha256=sha(binary),local_stage=str(ROOT),dependencies=['Ubuntu Python 3.8, numpy, scipy, h5py, matplotlib','PETSc/SLEPc, MPI shared libraries from build_receipt.json'],archive_scope='Lightweight evidence; original fields, all coefficient caches and all continuation vectors are retained locally rather than copied here.'))
    sources={}
    def add(p,name):
        p=Path(p)
        if p.is_file():sources[name]=p
    for p in sorted((ROOT/'outputs').iterdir()):
        if not p.name.startswith('s75_'):add(p,'outputs/'+p.name)
    for p in sorted((ROOT/'work').glob('*.log')):
        if not p.name.startswith('s75_'):add(p,'run_logs/'+p.name)
    runs=set();vectors=set()
    for row in audit['curves']:
        p=ROOT/'outputs'/row['curve'];d=json.loads(p.read_text());native=[r for r in d['points'] if r['kind']==0]
        for r in [native[0],native[len(native)//2],max(native,key=lambda r:abs(r['growth']))]+d['physical_tips']:vectors.add(r['vector_file'])
        for c in d['original_checks']:runs.add(c['receipt']['run'])
        for c in d['physical_tips']:runs.add(c['original_receipt']['run'])
        guide=json.loads(p.with_name(p.name.replace('_final','_v3_guide')).read_text());receipt=guide.get('seed',{}).get('receipt',{})
        if receipt.get('run'):runs.add(receipt['run'])
    for run in sorted(runs):
        for p in sorted((ROOT/'runs'/run).rglob('*')):
            if p.is_file() and p.suffix in ['.json','.ini','.csv','.log','.plt']:add(p,'evidence_runs/'+run+'/'+str(p.relative_to(ROOT/'runs'/run)))
    for p in sorted(vectors):
        path=Path(p);add(path,'selected_vectors/'+str(path.relative_to(ROOT)))
    for p in sorted((ROOT/'physical_guides/s70').glob('i*/receipt.json')):add(p,'physical_guide_receipts/'+str(p.relative_to(ROOT/'physical_guides')))
    for case in NEW:add(ROOT/'cases'/f's{case}'/'config.ini',f'case_configs/s{case}/config.ini')
    original=BASE/'75°/outputs/M3_all_neutral_curves_20261027'
    for name in ['selection_manifest.json','delivery_manifest.json']:add(original/name,'retained75_source/'+name)
    members=[];archive=repro/'numerical_evidence.tar.gz'
    with tarfile.open(archive,'w:gz',compresslevel=5) as tf:
        for name,path in sorted(sources.items()):
            tf.add(path,arcname=name,recursive=False);members.append(dict(path=name,source=str(path),bytes=path.stat().st_size,sha256=sha(path)))
    assert archive.stat().st_size<95_000_000
    # Read every archive member back and verify content rather than trusting creation alone.
    with tarfile.open(archive,'r:gz') as tf:
        for expected,member in zip(members,tf.getmembers()):
            assert member.name==expected['path'] and member.size==expected['bytes']
            assert hashlib.sha256(tf.extractfile(member).read()).hexdigest()==expected['sha256']
    js(repro/'evidence_archive_manifest.json',dict(created=STAMP,archive=archive.name,archive_sha256=sha(archive),archive_bytes=archive.stat().st_size,members=members,all_members_exact_readback=True))
    js(DEST/'delivery_review.json',dict(created=STAMP,numerical_audit_passed=True,original75_exact_csv_copy=True,visual_inspection=['M3_sweep_angles_comparison.png','M3_sweep_angles_upstream.png','M3_70deg_xbeta.png'],plots_have_all_selected_curves=True,csv_directory_has_no_unlisted_curves=set(p.name for p in (DEST/'curves').glob('*.csv'))==set(Path(r['output_csv']).name for r in curves),pending_dense_fpp_angles=[63,68,80],coarse74_excluded=True,frequency_cutoffs_not_claimed=True,remaining_scope='Other angles require selected FPP data or a separately labelled coarse-grid calculation; full modal/frequency completeness is not certified.'))
    print('DELIVERY_STAGED',len(curves),'curves',len(members),'evidence members',archive.stat().st_size,'compressed bytes',flush=True)

if __name__=='__main__':main()
