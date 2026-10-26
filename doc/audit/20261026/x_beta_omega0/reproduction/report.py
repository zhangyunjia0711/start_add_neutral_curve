"""Write a traceable Chinese delivery note after all numerical checks pass."""
import json,shutil,datetime
from pathlib import Path
from common import ROOT,snapshot
DEST=Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate/75°/outputs/M3_omega0_xbeta_20261026')
s=json.loads((DEST/'summary.json').read_text());audit=json.loads((ROOT/'outputs/band_audit.json').read_text());positive=json.loads((ROOT/'outputs/positive_beta_checks.json').read_text());native=json.loads((DEST/'tecplot_validation.json').read_text());snapshot()
for name in ['band_audit.json','positive_beta_checks.json']:shutil.copy2(ROOT/'outputs'/name,DEST/name)
text='''# 两工况固定 ω=0 的 x–β 中性曲线

已完成此前物理特征函数对应的 **M3 分支** 的固定实 ω=0 中性边界，使用与既有 β=−7 计算相同的冻结 v10 程序、局部壁面坐标、局部平行及零显式曲率模型。β 是原程序的无量纲波数，Lref=0.035 m。没有把文件夹名中的 linear_curvature 当作已经启用显式曲率项的证据。

两条边界在鼻部相接，但在 **x=1100 mm 计算域尾端仍然开口**。图和数据均没有人为封闭。该结果描述种子连通的 M3 分支，不是对全部特征谱或其他不连通分支的枚举。

|工况|有中性点的原生站位数|原生中性点|剖面插值鼻部点|最上游 x / mm|对应 β|x=1100 mm 的不稳定β绝对值区间|
|---|---:|---:|---:|---:|---:|---:|
'''
for r in s:text+=f"|{r['case']}|{r['station_count']}|{r['native_points']}|{r['profile_tips']}|{r['x_min_mm']:.6f}|{r['beta_at_xmin']:.6f}|{r['tail_abs_beta_band'][0]:.6f}～{r['tail_abs_beta_band'][1]:.6f}|\n"
text+='''
主图给出与既有 β=−7 的 M3 方向一致的负 β 支，浅色区是这条分支界定的不稳定带。signed 图和全部 `.plt` 文件还提供正 β 共轭支；已对每个原生站位检查矩阵对称及共轭向量的原多项式残差，另在原程序中独立复算种子、鼻部点和尾端两侧边界。正负支是同一静止实扰动的共轭表示，不应据此将物理模态数量翻倍。

## 文件

- `M3_omega0_xbeta_comparison.png/.svg`：全域叠加及鼻部局部图。
- `bluntplate_M3_omega0_xbeta.png/.svg`、`rough_M3_omega0_xbeta.png/.svg`：各自的单独曲线。
- `M3_omega0_xbeta_signed.png/.svg`：正负 β 的共轭表示。
- `M3_omega0_xbeta_flat_origin.png/.svg`：使用距各自平板起点的坐标比较。bluntplate 的平板起点 x≈0，rough 为14.858163 mm。主图横轴均是全局几何 x，不是鼻部弧长。
- `M3_omega0_xbeta_comparison.plt`：原生 Tecplot 二进制，4个 zone，依次为 bluntplate 负β、rough 负β、bluntplate 正β共轭、rough 正β共轭。
- 两个工况的独立 `.plt`：各有负β及正β共轭2个 zone。`.lay` 和 `.mcr` 默认显示负β支，正β数据仍在文件内。
- `.dat` 为可读的 Tecplot 输入，`.csv` 是相同物理点的表格。
- `*_calculation_evidence.json`、`band_audit.json`、`positive_beta_checks.json`、`tecplot_validation.json` 保留计算与验证证据。
- `reproduction/`：实际使用的最终脚本；完整原矩阵、原始向量和原程序输出保留在独立 WSL 计算目录。

CSV/Tecplot 的 `point_kind=0` 是原生 CFD 站位，`point_kind=1` 是相邻局部剖面的插值转折点；`conjugate=1` 表示由已验证的共轭变换生成。插值点不是新的 CFD 网格计算点。`station_coordinate` 对插值点是相邻原生站位的分数坐标；`case_id=1/2` 分别代表 bluntplate/rough。所有点的 ω 均严格为0，β为实数。

## 数值方法与验证

从原程序在 β=−16,0,16 导出空间二次多项式的矩阵，对其 β 依赖作二次重建；在两个种子另取 β=−7,0.37,−35 作独立矩阵检验。固定实 α 后解复 β 二次谱，令 Imβ=0，以伪弧长追踪拓扑；最终每一个原生站位均回到固定实β的空间问题，校正 Imα=0。中间的分片算子插值点不进入交付点表。鼻部转折点另对实际局部剖面插值重新组装原程序矩阵、验证两端与原生算子一致，并独立调用原程序复算。

'''
growth=max(r['maximum_growth_error'] for r in s);res=max(r['maximum_right_residual'] for r in s)
derivative=max(q['relative_derivative_difference'] for c in audit['cases'] for st in c['cross_sections'] for q in st['boundaries'])
checks=[q for case in ['bluntplate','rough'] for q in json.loads((ROOT/'outputs'/f'{case}_final.json').read_text())['original_checks']]
text+=f"- 最终点最大 |Imα|={growth:.3e}；最大原多项式右残差={res:.3e}。\n- 24个原生边界点的独立原程序抽验，最大α差={max(q['alpha_difference'] for q in checks):.3e}；两鼻部点另行原程序复算。\n- 16个横截面检查边界内外增长率符号，带内中点增长率均为正；从边界跟踪至中点的最小相邻原始向量 MAC={min(st['minimum_continuation_MAC'] for c in audit['cases'] for st in c['cross_sections']):.9f}。\n- 左向量解析导数与有限差分的最大相对差={derivative:.3e}。\n- 三个原生 `.plt` 经 Tecplot 360 EX 2024 R1 回读，{sum(q['double_values'] for q in native['files'])}个双精度数值逐值完全一致；二进制标识为 `#!TDV112`，批处理日志无ERROR。\n"
text+='''
## 适用范围与资料

数值残差小不等于 CFD 网格、法向网格或模型误差小。这里比较的是两个既有输入和同一局部 LST 模型的结果，没有额外作网格收敛或全局稳定性验证；鼻部最上游位置的微小差别不宜单凭这两张不同网格认定为显著物理变化。ω=0、β=−7 是不同二维切片；固定β−7在下游稳定，不能推出其他β的静止扰动也都稳定。

空间扰动的实 ω、实 β、复 α 以及增长率 −Imα 的约定参考 Nomura (2006), *Stability of Supersonic Swept-Wing Boundary Layers with or without Crossflow*, §2.3，pp.194–195，DOI:[10.2322/tjsass.48.193](https://doi.org/10.2322/tjsass.48.193)，[原文](https://www.jstage.jst.go.jp/article/tjsass/48/162/48_162_193/_pdf/-char/en)。该文支持约定，不是本工况的实测数据来源。

二次多项式谱、伴随线性化、奇异最高次矩阵及残差检查参考 [SLEPc 官方 PEP 文档](https://slepc.upv.es/release/documentation/manual/pep.html)，2026-10-26查阅。边界位置、共轭性质和比较数字来自本轮实际原矩阵与原程序计算，不将文献结果移作本工况结论。
'''
(DEST/'结果说明.md').write_text(text,encoding='utf-8')
readme='''# 复现说明

实际计算目录：`/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/neutral_xbeta_omega0_20261026`。
该目录保存两工况输入副本、全部原矩阵、特征向量和原程序输出；原输入路径及SHA见input_manifest.json。冻结程序SHA为634b45c5d30e076f5b86ccb136e0734c1905beedde108477dceb8a75adc3853c。

环境：Ubuntu WSL，Python3.8.10，NumPy1.24.4，SciPy1.10.1，h5py3.11.0，Matplotlib3.7.5；原程序使用/usr/local/mpich-4.1.2/lib；每个数值进程的OPENBLAS_NUM_THREADS、OMP_NUM_THREADS、MKL_NUM_THREADS均设为1。

现有计算可在原目录直接运行 `python3 work/audit.py`、`python3 work/positive_checks.py`，读取缓存并复核原输入；`python3 work/tecplot.py --verify` 需要原生回读文件，若已清理可在Windows输出目录重跑tecplot_roundtrip.mcr生成。

完整重算应使用新目录，禁止删除或覆盖现有结果：先在bootstrap.py中将ROOT末级目录名改为新名称，运行bootstrap.py；它会复制并校验原输入，同时生成common.py和run_adapter.py。把这里其他Python脚本复制到新ROOT/work。publish.py、tecplot.py、report.py中的DEST也应改为新的Windows交付目录。

按依赖顺序运行：solver.py（种子/β扫描）、prepare.py（1202站原矩阵）、contour.py bluntplate、contour.py rough、finalize.py bluntplate、finalize.py rough、audit.py、positive_checks.py、publish.py。两个contour以及两个finalize可在各自依赖完成后独立并行。

最后运行tecplot.py生成宏；在Windows使用native_tecplot.ps1（从75°任务目录执行，若重算需调整目标目录），它调用Preplot转换原生.plt并用Tecplot进行回读和布局保存。再运行tecplot.py --verify及report.py。宏使用相对ASCII文件名，避免目录名75°经ASCII转码导致路径错误。

contour.py --resume 只用于日志状态未完成的轮廓。finalize.py枚举全部实际站位穿越，不依赖预测步恰好落在网格上；每个穿越都通过原空间多项式求根。point_kind和conjugate含义见结果说明。最终原生站位必须各有两侧边界点；任何模式MAC、残差、稳定侧符号或完整性断言失败都应保留记录并诊断。

大体积输入、矩阵和向量保留在WSL，不随Git提交重复上传；这里提供最终结果、哈希、计算证据和复现代码。本次不是CFD网格加密，也不是全谱支枚举。
'''
(DEST/'reproduction/README.md').write_text(readme,encoding='utf-8')
for name in ['solver.py','tecplot.py','positive_checks.py','report.py']:shutil.copy2(ROOT/'work'/name,DEST/'reproduction'/name)
print('REPORT_COMPLETE',flush=True)
