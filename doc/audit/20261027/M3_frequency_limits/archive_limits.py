"""Archive verified frequency limits and reproducibility in the user's Git repo."""
import datetime,hashlib,json,platform,shutil,subprocess,sys
from pathlib import Path
import numpy as np,scipy,h5py
WIN=Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate/75°')
BASE=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°')
CALC=BASE/'neutral_m3_frequency_limits_20261027'
DEST=WIN/'outputs/M3_frequency_limits_20261027'
LOG=WIN/'outputs/工作日志_2026-10-27_两构型M3频率截止.md'
REPO=BASE/'bluntplate/neutral_beta_full_20261022/work/solver_source'
TARGET=REPO/'doc/audit/20261027/M3_frequency_limits'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def js(p,d):p.write_text(json.dumps(d,indent=2,ensure_ascii=False)+'\n')
def git(*args):return subprocess.check_output(['git','-C',str(REPO),*args],text=True)
def main():
    sys.path.insert(0,str(CALC/'work'));from common import snapshot
    protected=snapshot();prefix=str(TARGET.relative_to(REPO))+'/'
    for args in [('diff','--cached','--name-only','-z'),('diff','--name-only','-z'),('ls-files','--others','--exclude-standard','-z')]:assert all(p.startswith(prefix) for p in git(*args).split('\0') if p),'Another task has changes'
    selection=json.loads((DEST/'selection_manifest.json').read_text());audit=json.loads((CALC/'outputs/final_numerical_audit.json').read_text());validation=json.loads((DEST/'tecplot_validation.json').read_text());limits=json.loads((DEST/'frequency_limits.json').read_text())['limits'];assert len(selection['curves'])==19 and len(audit['curves'])==15 and len(limits)==4
    for r in selection['curves']:assert sha(r['path'])==r['sha256']
    for r in validation['files']:assert sha(DEST/r['file'])==r['sha256'] and sha((DEST/r['file']).with_suffix('.dat'))==r['source_sha256']
    evidence=DEST/'evidence';evidence.mkdir(exist_ok=True)
    for source in sorted((CALC/'outputs').glob('*.json')):
        # Structured diagnostics retain failures; large shared matrices stay in WSL.
        shutil.copy2(source,evidence/source.name)
    for source in sorted((CALC/'outputs').glob('*.log')):shutil.copy2(source,evidence/source.name)
    for name in ['bootstrap_receipt.json','input_manifest.json','final_numerical_audit.json']:shutil.copy2(CALC/'outputs'/name,DEST/name)
    scripts=DEST/'reproducibility';scripts.mkdir(exist_ok=True)
    for source in sorted((CALC/'work').glob('*.py')):shutil.copy2(source,scripts/source.name)
    for source in sorted(Path(__file__).parent.glob('*.py')):shutil.copy2(source,scripts/source.name)
    shutil.copy2(__file__,DEST/'archive_limits.py');receipts=set();vectors={};vector_sources=[]
    for k,r in enumerate(audit['curves']):
        d=json.loads(Path(r['path']).read_text());points=[p for p in d['points'] if p['kind']!=2]
        for c in d['original_checks']:receipts.add(c['receipt']['run'])
        for p in points:
            if p.get('original_receipt',{}).get('run'):receipts.add(p['original_receipt']['run'])
            vector_sources.append(dict(case=r['case'],f_kHz=r['f_kHz'],region=r['region'],xi=p['xi'],x_mm=p['x_mm'],beta=p['beta_r'],alpha_r=p['alpha_r'],alpha_i=p['alpha_i'],path=p['vector_file'],sha256=sha(p['vector_file'])))
        # Archive all final four small loops, tips, and representative wider-loop vectors.
        selected=range(len(points)) if r['f_kHz'] in [291.6,264.53,-240.2,-237.14] else sorted(set([j for j,p in enumerate(points) if p['kind']==1]+np.linspace(0,len(points)-1,8,dtype=int).tolist()))
        for j in selected:vectors[f'curve{k:02d}_point{j:04d}']=np.load(points[j]['vector_file'])
    for k,r in enumerate(limits):
        d=json.loads(Path(r['source']).read_text());receipts.add(d['original_check']['run']);vectors[f'critical{k}']=np.load(d['vector_file'])
    for case in ['bluntplate','rough']:
        d=json.loads((CALC/'outputs'/f'{case}_negative_secondary_pilot.json').read_text());receipts.add(d['initial_independent_check']['run'])
    np.savez_compressed(DEST/'critical_and_selected_curve_vectors.npz',**vectors);js(DEST/'point_vector_provenance.json',dict(points=vector_sources,scope='All new vectors remain in the isolated WSL run; final small loops, critical vectors and representative wider-loop vectors additionally archived in NPZ. Independent frozen-START vectors preserved as raw files.'))
    raw=[]
    for name in sorted(receipts):
        source=CALC/'runs'/name
        for rel in ['config.ini','summary.json','solver.log','out/neutral_field.plt','out/all_vectors.petsc','out/wall_local_profile_metadata.csv']:
            p=source/rel
            if not p.exists():continue
            target=DEST/'evidence_runs'/name/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target);raw.append(dict(path=str(target.relative_to(DEST)),source=str(p),sha256=sha(p),bytes=p.stat().st_size))
    js(DEST/'original_run_evidence_manifest.json',dict(runs=len(receipts),files=raw));js(DEST/'runtime_versions.json',dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,h5py=h5py.__version__,platform=platform.platform()))
    by={(r['case'],r['sign']):r for r in limits};p0=by['bluntplate',1]['f_kHz'];p1=by['rough',1]['f_kHz'];n0=abs(by['bluntplate',-1]['f_kHz']);n1=abs(by['rough',-1]['f_kHz'])
    report='''# 两构型 M3：延拓至频率截止

本轮将已有代表M3支扩展到正、负频率端，并补入原先未画出的下游失稳口袋。结果为**当前CFD网格、冻结局部壁面坐标及零显式曲率LST模型下的数值截止估计**；不额外镜像叠加共轭代表，不称为全部特征模态的完备稳定性证明。

## 结果

|构型|负端f (kHz)|正端f (kHz)|负端ω_nd|正端ω_nd|
|---|---:|---:|---:|---:|
'''
    for case,name in [('bluntplate','基准钝板'),('rough','连续曲率优化头部')]:
        neg,pos=by[case,-1],by[case,1];report+=f'|{name}|{neg["f_kHz"]:.6f}|{pos["f_kHz"]:.6f}|{neg["omega_nd"]:.8f}|{pos["omega_nd"]:.8f}|\n'
    report+=f'\n在相同计算模型下，优化构型正端频率上限约降低{100*(1-p1/p0):.2f}%，负端绝对值约降低{100*(1-n1/n0):.2f}%。这是频率支持范围的比较，不能据此单独确定N因子、输入扰动响应或转捩位置。\n'
    report+='''
频率沿用此前约定：f为物理循环频率、单位kHz；程序无量纲角频率 `ω_nd=2π(1000 f_kHz)Lref/Uref`，Lref=0.035m、Uref=1810.6307411853716m/s。kHz数值不是角频率rad/s。

## 图与Tecplot

- [频率扩展与两构型对比](M3_frequency_limits_comparison.png)：正端及负端，后者同时画出f=−112kHz的鼻部和下游两片区域。
- [最后闭合小环与临界点](M3_frequency_limits_last_loops.png)。
- [临界位置增长符号核验](M3_frequency_limits_growth_checks.png)。
- [全沿程中性根合并支](M3_frequency_limits_merger_map.png)：灰色虚线为原生算子插值导引。
- `M3_frequency_limits_neutral_curves.plt`：合并曲线与四个临界点；两个构型另有各自PLT/布局。曲线19 zones，临界点4 zones。所有源CSV/DAT与PNG/SVG同时保存。

## 对原图范围的修正

鼻部负频率闭环在约−114.13、−115.96kHz消失后，仍有下游失稳区域。原图只显示先前追踪到的闭环，不能用来判断整个M3家族的负端截止。本轮从定常向量独立续接发现下游正增长，冻结START复核通过；再沿实际中性方程追踪，在投影转折后重新连接到原中性根合并支。这里没有通过取共轭把另一张曲线复制进来。

最终负端估计因此更新到约−240.23、−237.15kHz；主图保留原鼻部参考，并新增下游−112、−150、−200kHz及临界前的小环。中性支和窗口的含义始终限定在所跟踪代表及已解析的连接口袋；复数共轭对称性本身不等于新增独立物理扰动。

## 验证与计算范围

采用原空间多项式的实α、β中性方程及β根合并条件，精确Jacobian/Hessian通过差分核验。逐截面与伪弧长方法联合穿过投影转折；算子插值仅作导引，交付曲线保留原生方程根和实际局部剖面重组的转折/补充点。

四个最终临界点由冻结START独立核验；在各临界位置，将频率向内/向外移动0.01、0.1、1kHz并重新优化β，增长率分别为正/负。该图是临界位置的局部β最大增长，不能冒充全域增长包络或N因子。

15条新曲线1267点通过原矩阵残差、零增长、原生覆盖、相邻向量和原求解器抽查；加4条旧参考共19条1815点。原生Tecplot回读逐值核对，另有四个临界点。详细证书在final_numerical_audit.json、tecplot_validation.json和evidence_runs。

负频率原生临界支联合覆盖基准第4～841截面、优化第2～361截面；正端主支覆盖基准第4～841、优化第5～361。极近附着线处的空间投影退化点及补充β扫描中的两段特征值投影失败均保留为未解析，没有把失败当稳定。正端外侧两构型各14个代表截面扫描β∈[−64,64]，26个完成该窗口、2个部分完成，实际计算点均无正增长；这是补充检查，不是任意β、所有谱分支的严格无根证明。

## 插值与物理精度限制

最后的小环包含实际重组的插值剖面；其尺度可小于原始沿程网格间距。它们表示指定插值模型内的数值延拓，**不表示完成CFD网格加密或网格收敛**。优化头部临界频率在相邻原生剖面间可出现下凹；没有把未通过内部增长检查的264.5kHz等尝试作为完整曲线交付。最终264.53kHz小环的实际内部采样均通过。

当前外形通过基流及局部壁面方向影响LST，显式曲率项仍为零。本轮没有实现含曲率模型或PSE/全局非平行验证；截止数值的小数位是该模型的求根精度，不是物理或CFD离散误差的置信区间。

## 复现和资料

reproducibility保存最终求解、修正、审查、绘图及归档脚本。evidence保留阶段结果和失败诊断；原输入、冻结二进制和原矩阵缓存按路径/SHA256引用，未重复上传大型共享数据。全部新向量保留在独立WSL目录，最终小环/临界向量及代表性向量另外归档为NPZ；冻结START独立运行的配置、日志、采样表和原始向量逐字节保存。

方法核验资料：[SLEPc PEP官方手册](https://slepc.upv.es/release/documentation/manual/pep.html)、[SciPy least_squares官方说明](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html)；空间扰动与局部曲率LST的适用范围见[JFM原文§2.4](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/instabilities-in-shockwaveboundarylayer-interactions-at-mach-6/795EAE9D2DAA6834476B70E5D8EB3DC2)。折点方程和导数由本任务的原方程推导并数值核验，不冒充文献已给出的本工况结果。
'''
    (DEST/'README.md').write_text(report);shutil.copy2(LOG,DEST/LOG.name);(DEST/'.gitattributes').write_text('* -text\n*.md text eol=lf\n*.py text eol=lf\n*.lay -whitespace\n*.log -whitespace\n*.csv -whitespace\n*.json -whitespace\n*.svg -whitespace\nevidence_runs/** -whitespace\n')
    files=[dict(path=str(p.relative_to(DEST)),sha256=sha(p),bytes=p.stat().st_size) for p in sorted(DEST.rglob('*')) if p.is_file() and p.name not in ['delivery_manifest.json','git_version.json']]
    js(DEST/'delivery_manifest.json',dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),files=files,file_count=len(files),total_bytes=sum(p['bytes'] for p in files),git_scope=prefix,protected_inputs=protected))
    TARGET.mkdir(parents=True,exist_ok=True)
    for source in sorted(DEST.rglob('*')):
        if not source.is_file() or source.name=='git_version.json':continue
        target=TARGET/source.relative_to(DEST);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target);assert sha(source)==sha(target)
    print('ARCHIVED',len(files),sum(p['bytes'] for p in files),len(receipts),'frozen runs',flush=True)
if __name__=='__main__':main()
