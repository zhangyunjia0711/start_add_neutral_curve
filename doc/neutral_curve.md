# 中性曲线模块

该模块借鉴 `LST_LJS` 中“沿参数逐步扫描、由增长率变号括区、迭代校正中性点”的工作流程，调用 START 自己的基流与特征值求解接口。它输出数值中性点，再用 Python 或 Tecplot 绘图。中性条件为时间模式的 `Im(omega)=0`，或空间模式的 `-Im(alpha)=0`；两种模式均以增长率为正表示不稳定。

## 运行条件

- 使用一个 MPI 进程；伴随问题不支持此模块。
- 扫描参数名为 `iloc`、`Re`、`alpha`、`beta`、`omega`，大小写不敏感；`x` 只作为求解状态输出，不能用作扫描轴。`iloc` 仅能作 outer 轴，并按整数站位处理。
- 扫描 `alpha`、`beta`、`omega` 时改变其输入实部；增长率来自求解得到的复特征值。时间问题的特征值为 `omega`，因此不能将 `omega` 选为扫描轴；空间问题相应禁止扫描 `alpha`。
- outer 与 inner 参数必须不同。固定站位的参数图应使用两个允许的波数、频率或 Reynolds 数作为扫描轴。
- 选择的参数范围应覆盖同一物理模态。特征值误差超限或相邻求解跳变过大的点被记为失败，不参与中性根和增长率图插值。该检查不能代替对特征函数、分辨率和模态连续性的物理判断。

## 菜单与配置

启动后使用中性曲线菜单，或在 INI 文件中启用 `[neutral]`。菜单编号保留原参考程序的对应入口，但实现调用 START 的求解器：

| 编号 | 功能 |
|---|---|
| 31 | 沿 outer 轴逐行扫描 inner，增长率变号括区后用带括区约束的割线法校正 |
| 32 | 同样逐行求根，使用二分法校正 |
| 33 | 扫描增长率场，同时输出中性根，供绘制增长率等高图 |
| 36 | 从中性种子点沿二维参数平面作伪弧长零线跟踪，可经过相对于单一参数的转折点 |
| 37 | 固定基流站位，在选定参数平面逐行扫描并求中性根 |

可直接运行的完整配置位于 `lst_test/lst_test/neutral_curve.ini`，包含基流文件、站位和 LST 模态设置。下面给出其中的中性扫描参数及可调选项。本次已验证 `beta=-3.25..-3.15`、`alpha=2.8..4.5` 的三个 outer 行，每行得到两个中性根；更换基流、站位或模态后，应重新确认范围与模态连续性。

```ini
[neutral]
enabled = T
task = 37
outer_parameter = beta
inner_parameter = alpha
outer_min = -3.25
outer_max = -3.15
outer_points = 3
inner_min = 2.8
inner_max = 4.5
scan_points = 31
growth_tol = 1.0e-7
parameter_tol = 1.0e-9
max_iterations = 60
max_points = 300
eigen_error_tol = 1.0e-7
max_eigen_jump = 0.5
continuation_step = 0.02
ds = 0.01
ds_min = 1.0e-5
ds_max = 0.05
direction = 1
output_prefix = out/neutral
```

数值参数的含义：

| 参数 | 含义 |
|---|---|
| `outer_min/max/points` | 外层参数取值范围及行数 |
| `inner_min/max`、`scan_points` | 每行内层参数范围及初始扫描点数 |
| `growth_tol`、`parameter_tol` | 增长率与参数校正精度 |
| `max_iterations` | 单个中性点或校正步骤允许的最大迭代数 |
| `max_points` | 弧长跟踪的最大点数 |
| `eigen_error_tol` | 特征值误差阈值 |
| `max_eigen_jump` | 相邻延拓求解允许的复特征值跳变上限 |
| `continuation_step` | 波数或频率参数的延拓步进上限，默认 `0.02`；Re 与 iloc 的专用规则见下文 |
| `start_outer`、`start_inner` | 弧长初始种子点的两轴坐标 |
| `start_branch` | 初始行存在多个根时选择的中性根分支 |
| `tangent_angle` | 初始切向角，单位为度，用于选择弧长跟踪方向 |
| `ds`、`ds_min/max` | 弧长初始、最小与最大步长；两轴量纲和数量级应适合相应步长 |
| `direction` | 弧长跟踪方向，取 `1` 或 `-1` |
| `output_prefix` | 输出文件前缀，含路径时该目录须事先存在 |

程序提供默认 `out/` 目录。若修改 `output_prefix`，应先创建父目录。

计算新参数点前，程序从已接受的模态状态出发作小步延拓。波数/频率每步变化受 `continuation_step` 限制；Re 每步使用基准 Re 的 5% 作为上限，`iloc` 按 1 个整数站位前进。首次用较小步长估计特征值变化，之后使用割线预测复特征值，再选择最接近预测目标的特征值。离散换站先在固定 inner 参数下逐站求解，不对站位作差分外推；随后在固定站位沿 inner 参数使用割线预测。按最近目标选择模态只在中性曲线任务期间开启，任务结束后恢复原有单点求解规则。延拓能降低跨模态跳变的风险，但仍应检查特征函数与网格收敛。

本次原版基准求解使用 `Ma=6`、`Re=79100`、`alpha=3.5`、`beta=-3.2059`、`iloc=270`，时间特征值约为 `omega=(-0.03933479, 0.01802792)`。该点的时间增长率大于零，是定位模态的参考点。以同一模态在 `beta=-3.2059` 扫描，已得到 `alpha≈2.986692664285` 和 `3.957146023035` 两个中性根，增长率绝对值分别约为 `8.3e-10` 和 `2.6e-11`，该次求解没有拒绝点。应以实际运行日志中的收敛结果为准。

## 输出文件

给定前缀 `out/neutral`：

- `out/neutral_roots.plt`：成功校正的中性点以及程序写出的状态信息。
- `out/neutral_field.plt`：扫描及根校正的求解评估样本，含重复点和非规则间隔；弧长任务是否输出此文件取决于是否进行了初始扫描。
- `out/neutral_ridge.plt`：每个 outer 行的增长率候选最大值。任务 31/32 在最优扫描点的邻域用最多 24 步黄金分割进一步细化；任务 33/37 只取扫描样本最优点。该文件不构成全参数区间的全局最大值证明。
- `out/neutral_run.log`：扫描参数、求解及终止状态。无根、求解失败、达到点数上限与正常完成应依据日志区分。

程序会保留失败前已接受的有效点，因此有 `.plt` 输出并不代表任务完整成功。出现求解失败时，批处理以退出码 1 返回；检查 shell 退出状态以及 `_run.log` 中各行状态、`rejected` 与最终 `result`。绘图只显示成功样本，不能代替这些检查。

Tecplot 数据为 POINT 格式，列顺序固定如下：

```text
VARIABLES = "outer" "inner" "iloc" "x" "Re" "alpha_r" "alpha_i" "beta_r" "beta_i" "omega_r" "omega_i" "growth" "branch" "status"
ZONE T="row_1"
...
ZONE T="row_2"
...
```

`outer/inner` 是所选扫描轴的实际值；`x`、`Re` 以及其余复数量是该次求解状态。`growth` 为上述统一增长率。`branch` 用于区分中性根分支，增长率场中的值为 0。`status=0` 表示求解成功，非零表示失败，不能作为有效场值绘图。`status=2` 表示求解评估失败，数值列中可能包含 NaN。运行日志的终止码另作解释：0 成功、1 未找到根、2 求解失败、3 达到上限、4 到达边界、5 闭合、6 步长降至下限、7 输入无效、8 退化；应结合日志上下文判读。

逐行扫描的 roots、field 和 ridge 使用每个 outer 行一个 `ZONE T="row_N"`；弧长任务的 roots 使用单个 `ZONE T="arclength"`，按记录顺序连线。空行或失败行应保留原有编号，绘图程序据此识别中断，不将两侧中性点直接相连。分支编号用于图形分组，不应仅凭该编号断定两个中性点属于同一物理模态。

## Python 绘图

依赖 Python 3、NumPy 和 Matplotlib，不需要 SciPy。在仓库根目录执行以下命令，可使用随程序提供的完整配置完成计算并绘图：

```bash
cd lst_test/lst_test
../../bin/START_LST -f neutral_curve.ini
# 求解正常结束后，绘制中性根与增长率背景：
python3 ../../tools/plot_neutral_curve.py out/neutral_roots.plt \
  --field out/neutral_field.plt \
  --x-label beta --y-label alpha --output out/neutral_preview.png
```

仅绘制中性根时省略 `--field`：

```bash
python3 ../../tools/plot_neutral_curve.py out/neutral_roots.plt \
  --x-label beta --y-label alpha --output out/neutral_roots.png
```

不传 `--output` 时打开 Matplotlib 窗口。使用 `--output` 时采用无界面的 Agg 后端，适合 WSL 或计算节点；输出图的父目录须存在。

脚本读取 VARIABLES 与多条 ZONE，支持 Fortran 的 D 指数，保留失败点和空行信息。中性点按 `branch` 分段，并在缺失行或失败点处断线；弧长点保持文件中的顺序。增长率场优先使用带掩码的直角网格，整格屏蔽触及失败样本的单元；非规则数据使用三角剖分，合并重复坐标；只要重复点中有失败记录，该坐标即被屏蔽。同一行按 inner 坐标排序建立邻接，跨行只连接原始编号相邻的行，并屏蔽失败节点及过大间隔；无有限坐标的失败评估使对应行保守留白。颜色为正负增长率，黑线为求得的中性点。掩码可能留下空白，这是数值数据缺失或不满足局部邻接条件的区域，不宜把空白解释为稳定。

空文件或只有表头的文件可解析；若中性点与增长率场均没有成功数据，脚本以退出码 2 报告原因，不生成误导性的空结果图。只有场数据成功时仍可绘制背景，并标注没有成功中性根。

## 验证与限制

本模块的验证包含数值求根/弧长方法的解析函数测试、11 项绘图输入解析与失败缺口测试，以及 START 实际算例的编译和求解测试。

实际时间模式任务 37 在 `Ma=6`、`Re=79100`、`iloc=270`、`beta=[-3.25,-3.20,-3.15]`、`alpha=2.8..4.5` 上完成：三个 outer 行各得到两个中性根，全部六个根的增长率绝对值小于 `6.5e-9`，140 次评估、545 次内部求解、0 个拒绝点，最终状态为 `converged`。对应结果为审计目录 `real_cases/task37/out/neutral_roots.plt`，预览为同目录 `neutral_preview.png`，已检查图中两条中性分支、增长率背景和坐标标签。所有实际执行记录以审计目录的过程日志为准。

每行扫描收集范围内发现的所有根，数量不限于两个。`ridge` 只代表扫描样本最优值或最优样本邻域的局部细化结果，不保证全局最大增长率。

扫描只能发现离散样本之间的变号根，或满足增长率容差的采样点；若中性曲线与扫描线相切而不变号，可能漏检。应根据增长率场增加采样密度，必要时使用任务 36 作弧长延拓。发现根后，应检查网格收敛、模态连续性和参数无量纲定义，再用于物理结论。

任务 36 要求 outer 与 inner 两个区间均有正宽度；单行固定 outer 配置只适用于逐行扫描任务。
