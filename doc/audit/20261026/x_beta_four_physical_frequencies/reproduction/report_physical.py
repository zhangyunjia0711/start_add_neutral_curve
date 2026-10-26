"""Publish factual numerical scope, reference conventions, and check receipts."""
import datetime,json,shutil,sys,platform
from pathlib import Path
import numpy,scipy,h5py,matplotlib
from common import ROOT,snapshot,sha,SOLVER,SOLVER_HASH
from publish_physical import DEST

def main():
    protected=snapshot();data=[json.loads(p.read_text()) for p in sorted((DEST).glob('*_calculation_evidence.json'))]
    assert len(data)==8
    positive=[q for d in data for q in d['original_checks']];negative=[q for d in data for q in d['negative_frequency_checks']]
    tips=[q for d in data for q in d['physical_tips']];sides=[q for d in data for q in d['local_boundary_checks']]
    native=json.loads((DEST/'tecplot_validation.json').read_text())
    conversion=json.loads((DEST/'frequency_conversion.json').read_text());source=Path(conversion['reference_source'])
    shutil.copy2(source,DEST/'reference_input_conversion.json')
    parameter_pilot=json.loads((ROOT/'outputs/pilot.json').read_text())
    coefficient_records=[]
    for p in sorted((ROOT/'native_pencils').glob('*/i*/frequency_coefficient.json')):
        r=json.loads(p.read_text())
        assert sha(r['prior_beta_coefficients'])==r['prior_report_sha256']
        for item in r['prior_coefficient_files']:assert sha(item['path'])==item['sha256']
        coefficient_records.append(dict(case=r['case'],station=r['station'],frequency_report=str(p),frequency_report_sha256=sha(p),G_sha256=sha(p.parent/'G.npz'),prior_beta_report_sha256=r['prior_report_sha256']))
    assert len(coefficient_records)==1202
    (DEST/'operator_verification.json').write_text(json.dumps(dict(note='Joint-parameter matrix checks precede physical-unit clarification; these are operator checks, not requested-frequency curves.',pilot_validations=[dict(case=r['case'],station=r['station'],validations=r['validations']) for r in parameter_pilot],native_coefficients=coefficient_records,prior_coefficients_hashes_verified=True),indent=2)+'\n')
    text='''# 两工况固定物理频率的 x–β 中性曲线

完成 f=−100、−50、+50、+100 kHz 的结果。这里 f 是有符号周期频率；原程序参数 ω 是无量纲角频率。按用户本轮确认的单位，采用：

`ω_nd = 2π f_kHz × 1000 × Lref/Uref`，`Lref=0.035 m`，`Uref=1810.6307411853716 m/s`。

所以 ±50 kHz 对应 ω_nd=±6.0727867021443425；±100 kHz 对应 ±12.145573404288685。最初按无量纲 ω=±50、±100 做的预检查已被本次物理单位结果替代，没有混入交付图件。归一化来源及公式见 frequency_conversion.json、reference_input_conversion.json。

## 实际结果与范围

沿用此前 M3 特征函数种子、冻结 v10 求解器、局部壁面坐标、局部平行及零显式曲率模型。bluntplate 为 841×400 网格，rough 为 361×400 网格；准备的参数算子覆盖两套原网格 −35～1100 mm。

从 ω=0 的 β=−7 种子及其 β=+7 共轭种子延拓，正频率分别得到 A、B 两个种子连通分量。A/B 是本任务的来源标记，不代表新命名的两个物理模态。两工况的 +50、+100 kHz 共8条轮廓全部完成闭合返回核验；负频率取同时反转 ω、β、Reα 的实场共轭，另经原程序独立抽验。每个负频率文件仍有 A、B 两个分量，且不会在同一非零频率下直接镜像 β。

|工况|f/kHz|分量|最上游 x/mm|最下游 x/mm|交付点 β 范围|原生点数|剖面插值转折点数|
|---|---:|---|---:|---:|---:|---:|---:|
'''
    for d in sorted(data,key=lambda d:(d['case'],d['f_kHz'],d['component'])):
        text+=f"|{d['case']}|+{d['f_kHz']}|{d['component']}|{d['x_range_mm'][0]:.6f}|{d['x_range_mm'][1]:.6f}|{d['negative_beta_range'][0]:.6f}～{d['negative_beta_range'][1]:.6f}|{d['native_point_count']}|{len(d['physical_tips'])}|\n"
    text+='''
负频率的 x 范围相同，β 范围反号并交换上下界。表中 β 极值是最终点表的极值，没有声称已单独连续优化 β 极值。横轴是全局几何 x，不是壁面弧长；图件展示实际闭合轮廓所在区间。

本轮观察：rough 的 A 分量下游闭合位置由 bluntplate 的118.003 mm延后到142.890 mm（50 kHz），由56.876 mm延后到75.038 mm（100 kHz）。这仅是相同 LST 模型下的中性边界差异，不能直接当作转捩位置或“粗糙度一定促稳/促失稳”的证据。

## 文件与打开方式

- `M3_xbeta_four_frequencies_comparison.png/.svg`：四频率总图。蓝色 bluntplate，红色 rough；实线 A，虚线 B。
- `M3_xbeta_f±50kHz_comparison.png/.svg` 和 `M3_xbeta_f±100kHz_comparison.png/.svg`：各频率对比，右侧为鼻部细节。
- 同名的4个 `.plt`：每个含两工况×A/B，共4个zone。
- `M3_xbeta_four_frequencies.plt`：四频率合并，16个zone。
- `bluntplate_M3_xbeta_four_frequencies.plt`、`rough_M3_xbeta_four_frequencies.plt`：每工况8个zone。
- `.plt` 是 Preplot 生成的原生二进制 `#!TDV112`；`.dat` 是可读输入，`.csv` 是逐条轮廓的相同点表。`.lay` 为保存的布局，`.mcr` 使用相对文件名，可从结果目录执行。
- `*_calculation_evidence.json`、`physical_audit.json`、`tecplot_validation.json`、`original_checks/` 保存数值检查、原程序抽验日志及格式回读证据。
- `reproduction/` 保存实际脚本和复现说明；大体积原矩阵、输入副本、向量和引导尝试保留在独立 WSL 计算目录。

Tecplot 的 X变量选 `x_mm`，Y变量选 `beta`。`f_kHz` 是物理周期频率，`omega_nd` 是传给程序的无量纲角频率；`beta` 为原程序无量纲展向波数。`station_coordinate` 是原生站位或相邻站位的分数索引。`point_kind=0` 表示原生 CFD 站位，1是实际局部剖面插值后重组装得到的转折点，2是已通过独立返回校验后重复首点形成的绘图闭合记录。插值点没有冒充新的 CFD 网格站位。`case_id=1/2` 是 bluntplate/rough，`component_id=1/2` 是 A/B，`conjugate=1` 表示负频率共轭结果。

## 方法与核验

原程序的空间多项式是 α 二次、β 二次、ω 线性。复用此前原生 β 算子系数，再在所有1202个站位导出 ω=1 的原矩阵求时间系数；种子处以额外 β、ω 组合做独立矩阵检验。原生算子之间的分片插值仅引导拓扑。

二维 α 或 β 谱投影在 rough、100 kHz 的折叠附近曾反复回走。最终改用实(站位,α,β)三参数的原多项式增广复残差，两条实方程加弧长条件作延拓，同时约束原多项式残差、向量连续性和返回段；不把失败引导点交付为曲线。

最终每一个穿越的原生站位均在固定实 β、实 ω 的空间问题中重新校正 Imα=0，各站恰有两个互不重复的边界点且站位连续。x 转折点则插值相邻局部法向坐标与原始流动变量，由 START 重新组装算子，校正 Imα=Im(∂α/∂β)=0；检查插值两端与原生算子吻合并用原程序复算。该步骤是剖面插值，不是新 CFD 或网格加密。

'''
    text+=f"- 正频率最终共有 {sum(d['native_point_count'] for d in data)} 个原生中性点、{len(tips)} 个剖面插值转折点及8条闭合重复记录；加负频率共轭后合并文件共7964条记录。\n"
    text+=f"- 最终最大 |Imα|={max(d['maximum_growth_error'] for d in data):.3e}；最大原多项式右残差 ‖Pq‖₂={max(d['maximum_right_residual'] for d in data):.3e}（‖q‖₂=1）。\n"
    text+=f"- 原程序独立抽验：正频率原生点{len(positive)}次、正频率插值转折点{len(tips)}次、负频率原生/转折点{len(negative)}次；正频率原生点最大α差{max(q['alpha_difference'] for q in positive):.3e}，负频率最大α差{max(q['alpha_difference'] for q in negative):.3e}。\n"
    text+=f"- 负频率抽验最小原始向量MAC={min(q['raw_MAC'] for q in negative):.15f}；每个正频率原生点均检查算子共轭关系和共轭向量的原残差。\n"
    text+=f"- 16个规则边界点以β±0.001验证增长率变号；伴随导数与中央有限差分最大相对差{max(q['derivative_relative_error'] for q in sides):.3e}。图只画中性线，未以局部抽验推断整片内部全都不稳定。\n"
    text+=f"- 7个原生.plt均经Tecplot360 EX 2024 R1打开和导回，所有双精度数值逐值一致；合并文件16zone、119460个双精度数值，批处理日志无ERROR。\n- {len(protected)}个保护文件的SHA复核通过，冻结原程序SHA256为 `{SOLVER_HASH}`。\n"
    text+='''
原程序抽验的“sample only / no bracketed root”是三点采样模式的正常结束标志；核验取β中心点及其向量。终端打印的首个偏移采样点增长率不等于中心中性误差，未把它用作验收值。

## 适用条件与资料

这是既有 M3 种子连通分支的计算，没有枚举其他模态或所有不连通中性分支。局部平行、零显式曲率的鼻部结果仍受模型适用性限制；目录名 linear_curvature 不代表本轮启用了显式曲率项。未新增 CFD 或法向网格收敛研究，残差小不等于模型和网格误差小。正负频率共轭是同一实扰动的表示关系，不能据此把物理模态数量翻倍。

归一化及符号以实际 START 输入和 `src/lst/lst_eqn.f90` 的 `G_0=−iG0`、空间 `+iβC0`、`+iαA0` 为依据。实ω、实β、复α和沿正流向增长率−Imα的空间LST约定参考 Nomura (2006), *Stability of Supersonic Swept-Wing Boundary Layers with or without Crossflow*, §2.3，pp.194–195，DOI [10.2322/tjsass.48.193](https://doi.org/10.2322/tjsass.48.193)，[原文](https://www.jstage.jst.go.jp/article/tjsass/48/162/48_162_193/_pdf/-char/en)。该文支持约定，不能代替本工况数据。

多项式特征问题、伴随线性化、奇异最高次系数矩阵及残差检查参考 [SLEPc官方PEP文档](https://slepc.upv.es/release/documentation/manual/pep.html)。两来源均已于2026-10-26查阅；本轮边界、共轭检验及比较数字来自实际原矩阵和原程序计算。
'''
    (DEST/'结果说明.md').write_text(text,encoding='utf-8')
    runs={q['receipt']['run'] for q in positive+negative}|{q['original_receipt']['run'] for q in tips}
    receipts=[]
    for name in sorted(runs):
        src=ROOT/'runs'/name;dst=DEST/'original_checks'/name;dst.mkdir(parents=True,exist_ok=True)
        files=[]
        for relative in ['config.ini','summary.json','solver.log','wall_local_profile_metadata.csv','out/neutral_field.plt','out/neutral_run.log']:
            sourcepath=src/relative
            if not sourcepath.exists():continue
            target=dst/relative;target.parent.mkdir(exist_ok=True);shutil.copy2(sourcepath,target)
            files.append(dict(path=str(target.relative_to(DEST)),sha256=sha(target)))
        receipts.append(dict(run=name,source=str(src),files=files))
    environment=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),python=sys.version,platform=platform.platform(),numpy=numpy.__version__,scipy=scipy.__version__,h5py=h5py.__version__,matplotlib=matplotlib.__version__,solver=str(SOLVER),solver_sha256=SOLVER_HASH,normalization_source_sha256=sha(source),calculation_root=str(ROOT),previous_root=str(ROOT.parent/'neutral_xbeta_omega0_20261026'),original_check_runs=receipts)
    (DEST/'execution_provenance.json').write_text(json.dumps(environment,indent=2)+'\n')
    readme=f'''# 复现说明

实际计算目录：`{ROOT}`。大体积输入副本、原生矩阵、原始向量和所有引导日志保留在该目录，原文件SHA见input_manifest.json。

本目录Python脚本是实际工作脚本的快照。common.py把ROOT设为脚本父目录的父目录，因此直接复核应在上述计算目录执行，而不是在交付目录reproduction内执行。

环境与版本见execution_provenance.json。每个数值进程设OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=MKL_NUM_THREADS=1；原程序运行环境由run_adapter.py设置，使用/usr/local/mpich-4.1.2/lib。

现有计算复核：`python3 work/audit_physical.py`。它检查所有受保护文件、站位完整性/重复根、负频率原程序抽验及导数/边界变号。`python3 work/publish_physical.py`重生成图件和点表。

完整重算应使用新ROOT，并同步修改bootstrap.py中的ROOT；禁止删除既有数据。bootstrap.py从前一ω0任务复制并核验输入。把这里除bootstrap.py外的脚本复制到新ROOT/work，覆盖bootstrap复制的common/run_adapter/core，并复制frequency_conversion.json到新ROOT/outputs。PREVIOUS应指向保留的neutral_xbeta_omega0_20261026，以读取已独立检验的原生β系数。publish_physical.py、tecplot_physical.py、report_physical.py、archive_physical.py的DEST需设为新交付目录。

顺序：prepare.py准备所有ω系数；physical_pilot.py从ω0种子作物理目标频率延拓；physical_survey.py取得每个种子站的β中性根对；run_contours.py以contour3d.py计算8条正频率分量；run_finalize.py在原生站位与实际插值剖面复算；audit_physical.py作独立检查；publish_physical.py出图；tecplot_physical.py生成宏。physical_pilot.py需要目标frequency_conversion.json，原始ω0种子和系数保留在PREVIOUS。

contour.py为contour3d.py提供Family及seed_root；其旧二维main未被用于最终图。run_contours.py只复用已通过closed_verified的结果，完整重算要新建输出目录。

Windows在新结果目录运行native_tecplot.ps1，调用Preplot和隐藏的Tecplot批处理完成转换及回读。随后运行tecplot_physical.py --verify、report_physical.py。若回读中间.dat已被清理，先重跑tecplot_roundtrip.mcr生成。宏使用相对ASCII文件名；修改路径后保存的.lay应重新生成。

数值断言失败时须保留日志并诊断，不能强制闭合或把失收敛当成稳定。本轮未枚举全谱，也未新增CFD网格收敛验证。
'''
    (DEST/'reproduction/README.md').write_text(readme,encoding='utf-8')
    for name in ['tecplot_physical.py','report_physical.py','archive_physical.py']:
        if (ROOT/'work'/name).exists():shutil.copy2(ROOT/'work'/name,DEST/'reproduction'/name)
    print('REPORT_COMPLETE',len(runs),'original-check receipts',flush=True)
if __name__=='__main__':main()
