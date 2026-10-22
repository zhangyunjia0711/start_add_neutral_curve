# 局部壁面坐标、谱支跟踪与有限 β 域扫描

本文补充 [中性曲线基础说明](neutral_curve.md)。下面选项保留旧配置的默认行为；只有显式启用时才生效。已验证环境为单 MPI 进程、复数 PETSc 3.22.1 / SLEPc 3.22.0。

## 坐标模型

`[lst] profile_coordinates` 默认为 `global`。`wall_local` 按站位计算壁面切向与向外法向，将二维速度转到该基底，投影法向坐标，并重新计算法向导数。平直、对齐的恒等变换保留原剖面及导数；坐标元数据写入运行目录的 `wall_local_profile_metadata.csv`。

`wall_local` 只支持正问题、空间模式、`LST_local`、`COMOVING_CS`、一个 MPI 进程和完整法向域；拒绝伴随、时间模式、LPSE 以及非零显式曲率。法向投影非单调等无效几何会停止计算。

这是每站冻结基底、零显式曲率的局部切平面模型。缺失的切向基流导数置零；旋转坐标没有补上曲率与非平行项。因此它修正剖面方向，并不构成完整鼻部稳定性模型。输出 `x` 仍是原网格壁面全局坐标；空间特征值 `alpha` 对应局部切向，不能把鼻部 `-Im(alpha)` 直接当成沿全局 x 的增长率。

## 跟踪选项

这些选项均位于 `[neutral]`：

| 选项 | 默认 | 实际作用 |
|---|---|---|
| `adaptive_steps` | `F` | 连续参数步失败时减小步长；不放宽残差、跳跃或向量阈值。 |
| `minimum_step` | `1e-5` | 自适应连续步长下限。离散 `iloc` 不能无限细分。 |
| `vector_tracking` | `F` | 用独立保存的原始右特征向量跟踪谱支，在本次已收敛候选中选择满足门槛且重叠最大的根。 |
| `minimum_overlap` | `0.95` | 归一化欧氏平方重叠下限，范围 `(0,1]`。这是跟踪参数，不是通用物理判据。 |
| `retry_eigenpairs` | `0` | `0` 或 `1` 不扩充；大于 1 时，在有参考向量但本步未通过向量选择后，同一点只重算一次。只支持空间 PEP 且须 `vector_tracking=T`。 |
| `sample_only` | `F` | 只允许 `task=33`；输出指定采样点，不求中性根，也不细化峰值。 |
| `reverse_outer` | `F` | 逐行扫描任务按 outer 从大到小执行；上下界仍按 `min ≤ max` 填写。 |
| `reverse_inner` | `F` | `sample_only` 时，inner 从大到小采样。 |
| `anchor_row_seed` | `F` | `sample_only` 时，使用本行首个采样点的特征值作为下一行初值；该点失败则没有新锚点。默认用本行有效采样中增长率最高点的特征值。 |

当前重叠定义为 `|q_refᴴ q|² / (‖q_ref‖² ‖q‖²)`，对缩放与整体复相位不敏感。这里没有物理能量权重或跨网格插值，不能直接比较不同法向网格的原始向量。仅有残差合格也不足以保证非正规问题的特征值前向误差或物理模态身份；仍需分辨率、模态形状和正反向跟踪检查。[SLEPc 可靠性说明](https://slepc.upv.es/release/documentation/manual/eps.html#reliability-of-the-computed-solution)

普通求解可显式请求一个候选：

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  ./bin/START_LST -f config.ini \
  -pep_tol 1e-12 -pep_max_it 1000 -pep_nev 1 -pep_ncv 8 \
  -st_pc_factor_mat_solver_type superlu
```

示例使用的串行 SuperLU 须由 PETSc 启用；也可使用原默认 SuperLU_DIST。INI 中 `retry_eigenpairs=8` 会在下一次重试中、读取命令行选项之后覆盖 `nev=8`、`ncv=max(24,2*8+8)=24`，并使用参考向量的独立副本作为初始空间。求解返回即清理一次性引用；再下一次求解恢复原维数。扩充不保证全谱覆盖，`eigen_error_tol`、`max_eigen_jump` 和 `minimum_overlap` 始终保留。候选序号只属于本次求解，不是持久谱支编号。[SLEPc PEP 说明](https://slepc.upv.es/release/documentation/manual/pep.html)

## 扫描片段与结果解释

[配置片段](neutral_tracking_example.ini) 必须合并到已有、可运行的工况配置中；它没有基流路径、全局参数和有效特征值种子，不能单独运行。修改坐标模型后应重新求有效种子。

空间模式保持实 `omega` 和实 `beta`，增长率定义为 `g=-Im(alpha)`。对声明的有限实 β 区间 `B` 和实际跟踪谱支集合 `J`，可在外部归约中定义 `G_B(x,omega)=max_{beta∈B,j∈J} g_j`。绘制的是 `G_B>0` 区域边界；多个不稳定频率区间之间的稳定空隙应保留，不能只连最小/最大频率。

二维扫描模块本身不会自动完成“沿程 × 频率 × β × 全谱”的穷尽搜索。需逐站、多种子计算，保留所有竞争峰，并核验 β/频率范围扩展、网格加密和谱支切换。报告中必须写明有限域、实际谱支与未解决区；β 最优值触及边界时，曲线可能被搜索范围截断。不同站位分别优化 β 得到的局部包络，也不等同于单个固定 β 扰动的 N 因子。

`*_field.plt` 中 `status=0` 才是有效样本。任何非零状态、缺失点或失败恢复都是 **unknown**，不能填成零增长、稳定点或跨过去插值。`sample_only` 成功采样后仍可能以 `no bracketed root`、退出码 1 结束，因为没有执行求根；须联合检查 field 状态、`rejected`、日志及 PETSc 错误，不能只按退出码分类。`ridge` 在这种模式下只是有效采样中的最高值，不是已优化的 β 极值。

## 仓库内短回归

```bash
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build --target test_slepc_candidate_retry test_raw_vector_resources -j4
ctest --test-dir build --output-on-failure \
  -R 'slepc_candidate_retry|raw_vector_resources'
```

测试无需工况文件：小对角问题检查候选/向量选择、独立副本、阈值拒绝、一次性维数恢复、重复求解和归一化的内存释放，并检查 EPS 明确拒绝扩充且不触发信号。CTest 使用 `-malloc_debug -malloc_dump -options_left`，不把未消费选项或退出分配残留当成成功。这些测试验证程序机制，不能代替真实流动的物理与网格收敛核验。

在线 SLEPc 文档随发行版更新；本文行为以本仓库实现和上述已测 3.22 环境为准。查阅日期：2026-10-22。
