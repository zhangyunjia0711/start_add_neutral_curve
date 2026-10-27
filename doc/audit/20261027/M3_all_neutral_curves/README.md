# 全部已交付M3中性曲线合图

本轮把此前的固定物理频率x–β中性曲线与频率截止延拓结果合在一起。没有重新求解、拟合曲线、插值生成新频率或额外镜像生成共轭曲线。

- [两工况并排比较，含全域和上游放大](M3_all_neutral_curves_comparison.png)：正负频率在每个坐标轴中同时绘制；两个工况用同一套频率颜色和线型。黑星是数值截止点。
- [两个工况直接叠在同一坐标轴中](M3_all_neutral_curves_both_geometries_overlay.png)：颜色表示频率；实线为基准钝板，虚线为连续曲率优化头部。左侧全域、右侧上游放大。
- 两图另存SVG。图例B、O分别表示该频率值只有基准、优化工况的已计算曲线。

共45条曲线、14807个曲线点、25个不同的物理频率值；基准23条、优化22条。f=−112kHz在每个工况有鼻部和下游两个已解析口袋，各自保留独立曲线，未用线段将它们拼接。f=0的两条曲线均在x=1100mm处开放，没有人为闭合。

物理循环频率f沿用kHz单位。无量纲角频率为 `omega_nd=2*pi*1000*f_kHz*0.035/1810.6307411853716`；数据中的beta为无量纲值，x使用原全局坐标mm。

## 数据与Tecplot

- `M3_all_neutral_curves.plt/.dat/.lay`：全部45条曲线和4个单独的临界点，共49 zones、14811点。
- `bluntplate_M3_all_neutral_curves.plt/.dat/.lay`：基准钝板。
- `rough_M3_all_neutral_curves.plt/.dat/.lay`：连续曲率优化头部，rough是原目录名称，不表示表面粗糙度。
- `curves/`：45份源CSV的逐字节副本。每条曲线保留为独立zone，频率、复波数、增长、原截面坐标、点类型、残差等列原样保存。
- `selection_manifest.json`：逐条来源、SHA256、点数、范围及4次重复引用去除记录。去重前先逐值核对重复曲线一致。
- `plot_validation.json`：各坐标轴实际绘制的曲线ID及坐标范围，全域图包含全部点。
- `tecplot_validation.json`：原生Tecplot 360 EX 2024 R1/Preplot写入后再回读，各变量逐值核对证书。
- `publish_all.py`、`native_tecplot.ps1`及`tecplot_roundtrip.mcr`：汇总、绘图、原生格式生成及回读复现文件。

来源为同工作区此前交付的 `M3_xbeta_extended_frequencies_20261027` 和 `M3_frequency_limits_20261027`。CSV先与各自交付清单的SHA256核对，再直接拷贝和读入；原数值结果不作修改。本轮只整理已有计算，不引入新的理论性结论。

“全部”指此前已交付的已验证曲线，并非全部特征谱的完备枚举。极近附着线退化、部分模态投影失败等原计算限制继续有效；尚未算出的其他口袋不能被合图证明不存在。截止附近的小环沿用此前的局部剖面插值，未做CFD网格收敛。本轮模型仍是零显式曲率的局部壁面LST。
