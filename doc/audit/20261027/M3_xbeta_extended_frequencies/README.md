# 两构型 M3：扩展频率范围的 x–β 中性曲线

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

## 新曲线范围

|构型|f (kHz)|x_min (mm)|x_max (mm)|β_min|β_max|
|---|---:|---:|---:|---:|---:|
|baseline_blunt_plate|-112|-30.977423|-25.691200|-22.893328|-20.842535|
|baseline_blunt_plate|-110|-31.749890|-24.359324|-23.032504|-20.166672|
|baseline_blunt_plate|+125|-28.022505|42.348018|-1.381670|9.890923|
|baseline_blunt_plate|+150|-26.275977|32.402404|1.498846|12.067668|
|baseline_blunt_plate|+175|-24.405581|24.583868|4.488356|14.212320|
|baseline_blunt_plate|+200|-22.343512|18.115258|7.571480|16.316640|
|continuous_curvature_nose|-112|-31.916879|-24.977563|-23.360180|-20.510488|
|continuous_curvature_nose|-110|-32.434569|-23.916805|-23.407750|-19.905601|
|continuous_curvature_nose|+125|-28.413022|57.163280|-1.317077|9.409043|
|continuous_curvature_nose|+150|-26.706283|43.548126|1.618148|11.478406|
|continuous_curvature_nose|+175|-24.808340|32.453276|4.681406|13.496109|
|continuous_curvature_nose|+200|-22.622247|22.932289|7.853822|15.435403|
