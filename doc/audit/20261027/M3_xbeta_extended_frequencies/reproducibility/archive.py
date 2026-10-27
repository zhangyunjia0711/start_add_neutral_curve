"""Archive verified extended-frequency deliverables and frozen-solver evidence."""
import csv,datetime,hashlib,json,platform,shutil,subprocess,sys
from pathlib import Path
import numpy as np,scipy,h5py
WIN=Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate/75°')
BASE=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°')
CALC=BASE/'neutral_xbeta_extended_frequencies_20261027'
DEST=WIN/'outputs/M3_xbeta_extended_frequencies_20261027'
LOG=WIN/'outputs/工作日志_2026-10-27_扩大频率范围中性曲线.md'
REPO=BASE/'bluntplate/neutral_beta_full_20261022/work/solver_source'
TARGET=REPO/'doc/audit/20261027/M3_xbeta_extended_frequencies'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def js(p,d):p.write_text(json.dumps(d,indent=2,ensure_ascii=False)+'\n')
def git(*args):return subprocess.check_output(['git','-C',str(REPO),*args],text=True)
def main():
    sys.path.insert(0,str(CALC/'work'));from common import snapshot
    protected=snapshot();prefix=str(TARGET.relative_to(REPO))+'/'
    for args in [('diff','--cached','--name-only','-z'),('diff','--name-only','-z'),('ls-files','--others','--exclude-standard','-z')]:
        assert all(p.startswith(prefix) for p in git(*args).split('\0') if p),'Another task has changes'
    audit=json.loads((CALC/'outputs/new_frequency_audit.json').read_text());selection=json.loads((DEST/'selection_manifest.json').read_text());validation=json.loads((DEST/'tecplot_validation.json').read_text())
    assert len(audit['results'])==12 and len(selection['source_files'])==30
    for r in selection['source_files']:assert sha(r['path'])==r['sha256']
    for r in validation['files']:assert sha(DEST/r['file'])==r['sha256'] and sha((DEST/r['file']).with_suffix('.dat'))==r['source_sha256']
    logs=DEST/'calculation_logs';logs.mkdir(exist_ok=True)
    for p in sorted((CALC/'outputs').glob('*.log')):shutil.copy2(p,logs/p.name)
    for name in ['run_summary_attempt1.json','run_summary.json','seed_extension_summary.json','new_frequency_audit.json']:
        shutil.copy2(CALC/'outputs'/name,DEST/name)
    for case in ['bluntplate','rough']:
        for sign in ['positive','negative']:
            p=CALC/'outputs'/f'{case}_{sign}_frequency_extension.json';shutil.copy2(p,DEST/p.name)
    scripts=DEST/'reproducibility';scripts.mkdir(exist_ok=True)
    for source in sorted((CALC/'work').glob('*.py')):shutil.copy2(source,scripts/source.name)
    for source in sorted(Path(__file__).parent.glob('*.py')):shutil.copy2(source,scripts/source.name)
    shutil.copy2(__file__,DEST/'archive.py')
    receipts=set();vectors={};point_sources=[];summaries=[]
    for r in audit['results']:
        case,f=r['case'],r['f_kHz'];d=json.loads(Path(r['final_file']).read_text());guide=json.loads(Path(r['guide_file']).read_text());receipts.add(guide['seed']['receipt']['run'])
        for c in d['original_checks']:receipts.add(c['receipt']['run'])
        for k,p in enumerate(d['physical_tips']):receipts.add(p['original_receipt']['run']);vectors[f'{case}_f{f:+d}_tip{k}']=np.load(p['vector_file'])
        for repair in d.get('native_crossing_repairs',[]):
            for trace in repair['traces']:receipts.add(trace['original_check']['run'])
        for p in d['points']:
            if p['kind']==2:continue
            path=Path(p['vector_file']);point_sources.append(dict(case=case,f_kHz=f,kind=p['kind'],xi=p['xi'],x_mm=p['x_mm'],alpha_r=p['alpha_r'],alpha_i=p['alpha_i'],beta=p['beta_r'],vector_path=str(path),vector_sha256=sha(path)))
        summaries.append(dict(case=case,physical_case='baseline_blunt_plate' if case=='bluntplate' else 'continuous_curvature_nose',f_kHz=f,points=len(d['points']),x_min_mm=d['x_range_mm'][0],x_max_mm=d['x_range_mm'][1],beta_min=d['negative_beta_range'][0],beta_max=d['negative_beta_range'][1],maximum_growth_error=r['maximum_growth_error'],maximum_right_residual=r['maximum_right_residual']))
    for p in json.loads((CALC/'outputs/neutral_seeds.json').read_text())['records']:
        assert p['f_kHz']==int(p['f_kHz'])
        vectors[f'{p["case"]}_f{int(p["f_kHz"]):+d}_seed']=np.load(p['vector_file'])
    np.savez_compressed(DEST/'new_seed_and_tip_vectors.npz',**vectors)
    js(DEST/'point_vector_provenance.json',dict(points=point_sources,scope='Native/fold vectors retained in the isolated WSL run; seed/tip vectors and independent frozen START vectors also archived here. Shared original CFD/matrices are referenced by hash.'))
    evidence=[]
    for name in sorted(receipts):
        source=CALC/'runs'/name
        for rel in ['config.ini','summary.json','solver.log','out/neutral_field.plt','out/all_vectors.petsc','out/wall_local_profile_metadata.csv']:
            p=source/rel
            if not p.exists():continue
            target=DEST/'evidence_runs'/name/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target);evidence.append(dict(path=str(target.relative_to(DEST)),source=str(p),sha256=sha(p),bytes=p.stat().st_size))
    js(DEST/'original_run_evidence_manifest.json',dict(runs=len(receipts),files=evidence))
    with (DEST/'new_curve_summary.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(summaries[0]));writer.writeheader();writer.writerows(summaries)
    js(DEST/'runtime_versions.json',dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,h5py=h5py.__version__,platform=platform.platform()))
    text='''# 两构型 M3：扩展频率范围的 x–β 中性曲线

新增 **−112、−110、+125、+150、+175、+200 kHz**。与原九频率合并为15个频率、两构型30条曲线，共13540点；新增12条均经数值验证闭合，零频原曲线保持开放至x=1100 mm。

- [新增范围与原±100 kHz对比](M3_xbeta_extended_frequencies_new_range_comparison.png)
- [15频率总图](M3_xbeta_extended_frequencies_all_frequencies.png)
- [正频率25～200 kHz](M3_xbeta_extended_frequencies_positive_frequency_detail.png)
- [合并原生Tecplot文件](M3_xbeta_extended_frequencies.plt)：30 zones。
- `bluntplate_M3_xbeta_extended_frequencies.plt`、`rough_M3_xbeta_extended_frequencies.plt`：各15 zones。

`bluntplate`是基准钝板；`rough`是用户确认的连续曲率优化头部，文件夹名不表示表面粗糙度。图中不额外叠加另一共轭代表。每个频率沿用从定常β=−7种子连接得到的同一代表家族。

## 数值观察与范围

在新增正频率+125、+150、+175、+200 kHz中，下游最右中性位置随频率上升依次前移；优化构型在相同频率下的最右中性位置仍比基准构型靠后。这是该家族在β上取最右边界的比较，不是固定β扰动的下游截止，也不直接表明增长率、N或转捩位置。

负频率补算−110、−112 kHz，鼻部闭环相较−100 kHz收缩。固定原生截面的负频率续接在约−112.4和−114.4 kHz停下，**并未证明整张网格或全部模态在−125/−150/−200 kHz无失稳**。最初考虑的这三个更负频率没有可交付的完整同家族边界，未画成曲线，也未发布临界频率数值。

沿用当前局部壁面切平面、零显式曲率LST，与原九频率一致。头部形状通过基流和局部壁面方向进入。这里没有实施含显式曲率的鼻部模型、非平行/PSE验证或全部模态搜索。现有原生网格差异仍需在机制研究中另做收敛评估。

## 变量与核验

输入频率是有符号循环频率f，单位kHz。无量纲角频率 `omega_nd=2*pi*(1000*f_kHz)*0.035/1810.6307411853716`；β沿用原程序的无量纲波数。扰动约定下空间增长量为−Im(α)。

CSV/DAT保留13变量：x_mm、beta、f_kHz、omega_nd、alpha_r、alpha_i、growth、x_flat_mm、station_coordinate、point_kind、case_id、right_residual、raw_MAC。point_kind=0表示原生截面，1表示实际局部剖面插值后重组方程的转折点，2为经过验证的闭合重复首点。原生算子插值只用于引导，不将其内部导引点作为新CFD计算点交付。

已核验实际转折范围内原生截面没有遗漏、根不重复、零增长误差、原方程残差、相邻特征向量连续性、β边界两侧增长符号、导数差分，以及冻结START独立抽查。优化构型−112 kHz第24截面由实际剖面续接补足。详见new_frequency_audit.json与evidence_runs。

三份PLT经Tecplot 360 EX 2024 R1/Preplot生成，文件头为#!TDV112；Tecplot重新读取并输出双精度ASCII后逐值一致，详见tecplot_validation.json。临时回读DAT核验后删除，源DAT、PLT、布局、证书与哈希保留。

方法约定的一手资料：局部空间LST扰动表达与曲率度量示例见 [JFM原文§2.4](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/instabilities-in-shockwaveboundarylayer-interactions-at-mach-6/795EAE9D2DAA6834476B70E5D8EB3DC2)。该文工况不同，仅用于方法范围说明，不借用其稳定性结论；本轮曲线来自上列实际求解证据。

## 复现与保存

reproducibility保留频率续接、伪弧长追踪、物理转折校正、原生补点、审查与绘图脚本。独立WSL运行目录、原始CFD/矩阵路径和哈希在input_manifest及各证据文件中；共享的大型原输入/矩阵未重复复制。新频率的原生向量仍保存在独立WSL运行目录，种子、转折点及独立START抽查原始向量另随本目录归档。

必要的失败日志和修正原因保留。原九频率CSV复制前后哈希一致，原CFD、旧数值结果和求解器未改写。
'''
    text+='\n## 新曲线范围\n\n|构型|f (kHz)|x_min (mm)|x_max (mm)|β_min|β_max|\n|---|---:|---:|---:|---:|---:|\n'
    for r in summaries:text+=f'|{r["physical_case"]}|{r["f_kHz"]:+d}|{r["x_min_mm"]:.6f}|{r["x_max_mm"]:.6f}|{r["beta_min"]:.6f}|{r["beta_max"]:.6f}|\n'
    (DEST/'README.md').write_text(text);shutil.copy2(LOG,DEST/LOG.name);(DEST/'.gitattributes').write_text('* -text\n*.md text eol=lf\n*.py text eol=lf\n*.lay -whitespace\n*.log -whitespace\n*.csv -whitespace\n*.json -whitespace\n*.svg -whitespace\nevidence_runs/** -whitespace\n')
    files=[]
    for p in sorted(DEST.rglob('*')):
        if p.is_file() and p.name not in ['delivery_manifest.json','git_version.json']:files.append(dict(path=str(p.relative_to(DEST)),sha256=sha(p),bytes=p.stat().st_size))
    js(DEST/'delivery_manifest.json',dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),files=files,file_count=len(files),total_bytes=sum(p['bytes'] for p in files),protected_inputs=protected,git_scope=prefix))
    TARGET.mkdir(parents=True,exist_ok=True)
    for source in sorted(DEST.rglob('*')):
        if not source.is_file() or source.name=='git_version.json':continue
        target=TARGET/source.relative_to(DEST);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target);assert sha(source)==sha(target)
    print('ARCHIVED',len(files),sum(p['bytes'] for p in files),len(receipts),'frozen runs',str(TARGET),flush=True)
if __name__=='__main__':main()
