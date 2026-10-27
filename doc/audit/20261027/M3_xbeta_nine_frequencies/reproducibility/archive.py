"""Archive verified nine-frequency figures, numerical evidence and reproducibility."""
import datetime,hashlib,json,shutil,subprocess,sys
from pathlib import Path
import numpy as np
WIN=Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate/75°')
BASE=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°')
CALC=BASE/'neutral_xbeta_extra_frequencies_20261027'
DEST=WIN/'outputs/M3_xbeta_nine_frequencies_20261027'
REPO=BASE/'bluntplate/neutral_beta_full_20261022/work/solver_source'
TARGET=REPO/'doc/audit/20261027/M3_xbeta_nine_frequencies'
LOG=WIN/'outputs/工作日志_2026-10-27_九频率xβ中性曲线补算.md'
FS=[-100,-75,-50,-25,0,25,50,75,100]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def js(p,r):p.write_text(json.dumps(r,indent=2,ensure_ascii=False)+'\n')
def git(*args):return subprocess.check_output(['git','-C',str(REPO),*args],text=True)
def main():
    sys.path.insert(0,str(CALC/'work'));from common import snapshot
    protected=snapshot()
    selection=json.loads((DEST/'selection_manifest.json').read_text());audit=json.loads((CALC/'outputs/new_frequency_audit.json').read_text());native=json.loads((DEST/'tecplot_validation.json').read_text())
    assert len(selection['source_files'])==18 and len(audit['results'])==8
    for r in selection['source_files']:assert sha(r['path'])==r['sha256']
    for r in native['files']:assert sha(DEST/r['file'])==r['sha256'] and sha((DEST/r['file']).with_suffix('.dat'))==r['source_sha256']
    prefix=str(TARGET.relative_to(REPO))+'/'
    for args in [('diff','--cached','--name-only','-z'),('diff','--name-only','-z'),('ls-files','--others','--exclude-standard','-z')]:
        names=git(*args).split('\0');assert all(p.startswith(prefix) for p in names if p),'Do not mix another task'
    # Preserve all stage logs, including failed attempts, and both failed/successful summary states.
    logs=DEST/'calculation_logs';logs.mkdir(exist_ok=True)
    for source in sorted((CALC/'outputs').glob('*.log')):shutil.copy2(source,logs/source.name)
    for name in ['run_summary_attempt1.json','run_summary.json','new_frequency_audit.json','bootstrap_receipt.json','neutral_seeds.json','input_manifest.json','runtime_versions.json']:
        source=CALC/'outputs'/name
        if source.exists():shutil.copy2(source,DEST/name)
    scripts=DEST/'reproducibility';scripts.mkdir(exist_ok=True)
    for source in sorted((CALC/'work').glob('*.py')):shutil.copy2(source,scripts/source.name)
    for source in sorted(Path(__file__).parent.glob('*.py')):shutil.copy2(source,scripts/source.name)
    shutil.copy2(__file__,DEST/'archive.py')
    # Export point provenance and original solver check receipts including raw vectors.
    receipts=set();point_vectors={};receipt_files=[]
    for r in audit['results']:
        case,f=r['case'],r['f_kHz'];d=json.loads(Path(r['final_file']).read_text());guide=json.loads(Path(r['guide_file']).read_text())
        receipts.add(guide['seed']['receipt']['run'])
        for check in d['original_checks']:receipts.add(check['receipt']['run'])
        for j,p in enumerate(d['physical_tips']):
            receipts.add(p['original_receipt']['run']);point_vectors[f'{case}_f{f:+d}_tip{j}']=np.load(p['vector_file'])
        for repair in d.get('native_crossing_repairs',[]):
            for trace in repair['traces']:receipts.add(trace['original_check']['run'])
        for side in ['lower','upper']:
            p=next(p for p in json.loads((CALC/'outputs/neutral_seeds.json').read_text())['records'] if p['case']==case and p['f_kHz']==f and p['side']==side)
            point_vectors[f'{case}_f{f:+d}_{side}_seed']=np.load(p['vector_file'])
    np.savez_compressed(DEST/'new_neutral_seed_and_tip_vectors.npz',**point_vectors)
    for name in sorted(receipts):
        source=CALC/'runs'/name;target=DEST/'evidence_runs'/name
        for rel in ['config.ini','summary.json','solver.log','out/neutral_field.plt','out/all_vectors.petsc','out/wall_local_profile_metadata.csv']:
            p=source/rel
            if not p.exists():continue
            dest=target/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)
            receipt_files.append(dict(path=str(dest.relative_to(DEST)),source=str(p),sha256=sha(p),bytes=p.stat().st_size))
    js(DEST/'original_run_evidence_manifest.json',dict(runs=len(receipts),files=receipt_files,scope='Frozen START seed, native spot-check, physical-tip and added-native-station receipts; raw files preserved byte for byte'))
    rows=selection['source_files']
    text='''# 两工况 M3：九频率 x–β 中性曲线

新增 **−75、−25、+25、+75 kHz**，与原 −100、−50、0、+50、+100 kHz 合并。两工况各九条曲线，合计18条；每一频率沿用同一个从零频 β=−7 种子追踪得到的代表家族，不额外叠加另一共轭代表。

## 图和 Tecplot 文件

- [全范围与上游对比](M3_xbeta_nine_frequencies_same_branch_comparison.png)
- [上游放大](M3_xbeta_nine_frequencies_same_branch_upstream.png)
- [鼻部放大](M3_xbeta_nine_frequencies_same_branch_nose.png)
- `M3_xbeta_nine_frequencies_same_branch.plt`：合并18 zones。
- `bluntplate_M3_xbeta_nine_frequencies.plt`：9 zones。
- `rough_M3_xbeta_nine_frequencies.plt`：9 zones。

每张图另有 SVG。Tecplot 选 X=`x_mm`、Y=`beta`；附带 `.lay` 可直接打开。布局中的线型只用于辅助区分频率：−75、+25 kHz 用点线，+75 kHz 用虚线，其余用实线；不表示增加了其他失稳带或共轭代表。每个频率的上下中性边界属于同一条曲线。

## 实际计算范围

原五个频率直接复用前一任务已验证数值，并核对源文件哈希；新增四个频率分别直接求解其带符号的频率，没有用曲线之间的几何插值替代计算。新增八条曲线均完成闭合；零频保留至计算域尾端 x=1100 mm 的开放边界，不人为连接两端。

| f (kHz) | bluntplate x 范围 (mm) | rough x 范围 (mm) | 状态 |
|---:|---:|---:|---|
'''
    for f in FS:
        a=next(r for r in rows if r['case']=='bluntplate' and r['f_kHz']==f);b=next(r for r in rows if r['case']=='rough' and r['f_kHz']==f)
        text+=f"|{f:+d}|{a['x_range_mm'][0]:.6f} ～ {a['x_range_mm'][1]:.6f}|{b['x_range_mm'][0]:.6f} ～ {b['x_range_mm'][1]:.6f}|{'计算域尾端开放' if f==0 else '闭合已验证'}|\n"
    total=sum(r['points'] for r in rows);newpoints=sum(r['points'] for r in audit['results'])
    maxgrowth=max(r['maximum_growth_error'] for r in audit['results']);maxres=max(c['absolute_original_residual'] for r in audit['results'] for c in r['matrix_checks'])
    checks=sum(len(r['original_independent_checks']) for r in audit['results']);minMAC=min(r['minimum_adjacent_guide_MAC'] for r in audit['results'])
    text+=f'''
## 方法和验证

空间 LST 与原工况、冻结二进制一致（SHA256：634b45c5d30e076f5b86ccb136e0734c1905beedde108477dceb8a75adc3853c），沿用局部壁面坐标、冻结切平面、零显式曲率模型。物理频率 f 为有符号周频率，单位 kHz；角频率换算 ω_nd=2π f_kHz×1000×0.035/1810.6307411853716。β 为原程序无量纲展向波数。

种子从之前连续频率路径接近目标的点出发，再在精确目标频率下校正 Im(α)=0。原矩阵之间的插值只用于生成追踪导引；最终发布点是原生 CFD 截面的中性根和实际局部剖面插值重组后校正的转折点。实际剖面插值不是新增 CFD 网格计算。

- 新增八条曲线共 {newpoints} 点；与原结果合计 {total} 点。
- 新增点最大 |Im α|={maxgrowth:.4e}；每条16个原矩阵复核点，共128个，最大 ||P(α)q||₂={maxres:.4e}（q 已作欧氏单位范数）。另保留各求解器的原残差，不将不同归一化定义的残差直接混比。
- 最低导引相邻向量 MAC={minMAC:.8f}；闭合返回截面的位置和向量匹配均通过。每个原生截面无缺漏、无重复中性根，物理转折点确定的整数截面范围与发布点完全一致。
- {checks} 份原求解器独立抽查收据通过 α 差异<1e−7、原向量 MAC>1−1e−7 和中性增长误差<1e−8；另包含种子、转折点和新增截面的独立校验。`evidence_runs/` 保留这些运行的配置、日志、中心点采样和原始向量。sample_only 日志中的“no bracketed root”不作为失败判据，实际检查中心点及原特征向量。
- 所取共同截面上下两根的增长率两侧符号符合失稳带方向，伴随梯度与有限差分一致性检查。
- 三份原生 Tecplot 二进制文件均以 Preplot 生成，魔数 #!TDV112；经 Tecplot 360 EX 2024 R1 读入再写出，所有双精度值逐值相等。合并文件18 zones、11026点，bluntplate 9 zones、7991点，rough 9 zones、3035点。
- 41份原输入和旧结果的保护哈希检查通过。本任务没有修改 Fortran 源码、冻结求解器或原五频率的数值。

## 鼻部转折校正与补点

三处原来的 α 消元转折校正停滞，保留原失败日志，改用原 bordered 多项式实/虚残差和 (α,β) Jacobian 子式为零的联合条件，避免在 α 投影敏感处除以 pᴴPαq。该条件由当前残差方程切向关系直接推得，最终仍由原空间矩阵与冻结 START 检查，而不是仅凭非线性优化器返回成功。

rough −75 kHz 的实际剖面转折从导引 xi≈40.378 移到 xi≈41.029，跨过第41原生截面。检查捕获该遗漏后，补算该截面两根 β≈−16.8461365611、−16.3471081448，分别经21个实际局部剖面步从第40截面连接至41截面，并通过原求解器验证。最终只发布修复后的79点曲线，未把原导引点作为物理解交付。

关于数值形式，查阅 [SLEPc 官方 PEP 手册](https://slepc.upv.es/release/documentation/manual/pep.html)，核对原多项式特征问题；联合校正使用有界最小二乘，其接口与终止条件核对 [SciPy 官方 least_squares 文档](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html)。两者于2026-10-27实际查阅，仅用于数值方法背景与接口，不作为本工况物理结果的外部证明。

## 数据、复现与边界

18份 CSV 和三份 DAT 保存每个点的 x、β、f、ω、α、增长率、原网格截面、点类型和残差等13个变量。point_kind=0 原生截面，=1 实际剖面转折，=2 经验证后为绘图重复首点。`selection_manifest.json` 列出数据来源与哈希。

`new_frequency_audit.json`、八份 `*_calculation_evidence.json`、`original_run_evidence_manifest.json`、`calculation_logs/` 和工作日志记录计算证据与失败/修正过程。`new_neutral_seed_and_tip_vectors.npz` 保留新增曲线的种子和转折点向量；完整原生剖面向量仍保留在原环境的独立 WSL 计算目录，其路径与哈希在审核清单中。原始 CFD 和原矩阵不重复上传，复现需要原环境的数据与依赖。

`reproducibility/` 保存实际使用的求解/校正/补点/审核/绘图脚本。原环境步骤为：建立隔离目录 → 目标频率种子及导引 → 原生/物理转折校正 → rough −75 kHz 的 native41 补点 → audit_extra.py → publish.py → Windows native_tecplot.ps1 → publish.py --verify → archive.py。bootstrap 后须使用归档的最终校正脚本（包含 physical_fold_bordered.py 和最终 tip_correction.py）；首次失败流程在日志中保留。

原零频边界至1100 mm仍开放，未由当前网格证明其下游闭合位置。本图只覆盖已跟踪的 M3 代表家族，不是对全部特征模态的排他性全谱声明。不同频率曲线相互交叉，不能据此把同一频率的一条失稳带解释为多个独立机制。
'''
    (DEST/'README_九频率中性曲线.md').write_text(text,encoding='utf-8')
    shutil.copy2(LOG,DEST/LOG.name)
    cleanup=[];receipt=DEST/'cleanup_receipt.json'
    if receipt.exists():cleanup=json.loads(receipt.read_text())['removed']
    for r in native['files']:
        p=DEST/(Path(r['file']).stem+'_roundtrip.dat')
        if p.exists():
            assert p.resolve().parent==DEST.resolve() and sha(p)==r['roundtrip_sha256']
            cleanup.append(dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p),reason='Verified exact read-back; native PLT, source DAT and certificate retained'));p.unlink()
    js(receipt,dict(removed=cleanup))
    # Normalize presentation/reproducibility text; original run receipts remain exact.
    for p in DEST.rglob('*'):
        if p.is_file() and 'evidence_runs' not in p.parts and p.suffix in ['.svg','.csv','.json','.lay','.md','.py','.ps1','.mcr']:
            s=p.read_text(encoding='utf-8-sig');p.write_text('\n'.join(v.rstrip() for v in s.splitlines()).rstrip()+'\n',encoding='utf-8')
    (DEST/'.gitattributes').write_text('*.plt binary\n*.log binary\nevidence_runs/** binary\n',encoding='ascii')
    for r in receipt_files:assert sha(DEST/r['path'])==r['sha256']
    files=[dict(path=str(p.relative_to(DEST)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(DEST.rglob('*')) if p.is_file() and p.name not in ['delivery_manifest.json','git_version.json'] and not p.name.startswith('工作日志_')]
    js(DEST/'delivery_manifest.json',dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),files=files,protected_inputs=protected,execution_metadata_excluded=['工作日志_*.md','git_version.json','delivery_manifest.json']))
    shutil.copytree(DEST,TARGET,dirs_exist_ok=True)
    for p in TARGET.rglob('*'):
        if p.is_file():p.chmod(0o644)
    assert snapshot()==protected and all(sha(TARGET/r['path'])==r['sha256'] for r in files)
    print('ARCHIVE_READY',len(files),'files',sum(r['bytes'] for r in files),'bytes',len(receipts),'original runs',flush=True)
if __name__=='__main__':main()
