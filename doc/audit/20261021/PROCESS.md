# 中性曲线模块移植过程记录

日期：2026-10-21（Asia/Shanghai）
参考程序：E:\work\code\LST_LJS
目标程序：/home/zaza/start/START_ALL_IN_ONE/START_ALL_IN_ONE

初步只读检查：目标已有 Git 仓库，HEAD=60445c4；origin=https://gitee.com/ltlx_at_tju/START_ALL_IN_ONE.git。
本地原有修改：CMakeLists.txt、src/atnpse/main.f90、src/lst/alst_eqn.f90、src/lst/lst_eqn.f90、src/npse/main.f90；未跟踪 lst_test/、stdlib/。
WSL Git 2.25.1 可用；Windows Codex 随附 Git 2.53.0.windows.3 可用；未发现常规 Git for Windows 安装。
未发现已检查祖先目录的 AGENTS.md。已并行安排两个只读代码分析任务。
用户要求先上传原版再修改；因 origin 为原作者仓库，已请求用户给出有权限的上传目的地；在上传完成前不修改目标程序源码。
从本记录建立起，Linux 操作经 run_logged.py 记录命令、输出、退出码。

## 原版保存结果

完整归档：original_full.tar.gz（包含 .git、原有未提交修改、算例、stdlib、原可执行文件及构建目录）。
SHA256：84f9d2f50c734f9b7cd9da83fcda57eb7db72f28c719d0900b5f8829b56517ff
Git 原版分支：backup/original-before-neutral-20261021
Git 原版提交：7b12c3c37a47ee9caa664c6e08815c5dd976a335
此提交保存项目已跟踪文件的本地修改和 lst_test；stdlib 自带独立 Git，历史和本地补丁分别保存为 stdlib_history.bundle、stdlib_local.patch，完整工作树位于 full 归档。
独立 Git 备份：original_history.bundle（原有历史），original_with_local_changes.bundle（含原版分支）。
文件校验清单：original_manifest.json。
当前 develop 分支及用户暂存区未切换、未改写，目标源码尚未修改。

## Git 下载/安装

依据 Git 官方安装页 https://git-scm.com/install/windows 的 winget 命令下载 Git.Git；仓库本身使用现有 WSL Git 2.25.1。
winget 下载并验证 Git 2.55.0.5 安装包；系统安装提示管理员权限，最终返回 “You cancelled the installation”，exit code 2。未宣称 Windows 安装成功，未重试管理员安装。将安装包另外保留到 E:\work\code\Git-download-20261021。
详见 git_install.log、git_download.log。

## 远程上传

已检查连接的 GitHub 账户 zhangyunjia0711：当前仅发现公开 FPP 仓库，与本任务无关，未向其写入。
用户尚未给出上传目的仓库；未对原作者 Gitee origin 执行 push。为遵守先上传再修改的明确顺序，等待仓库信息期间仅备份、分析与隔离验证。

## 目标程序接口与移植方案

入口 src/lst/module_lst.f90；新增独立 neutral 模块，复用 lst%solve(ome,alf,bat)，无需复制参考程序的 O-S 方程矩阵。
时间模式增长率 sigma=Im(omega)，空间模式 sigma=-Im(alpha)。
建议保留 31/32/33/36/37 五个功能入口，统一中性点括区、有限次割线/二分校正、续接预测、状态报告与 Tecplot 数据输出；另附绘图脚本和独立配置。
沿站位计算使用 assemble_lst_grid_profile/baseflow_profile 后 set_grid/set_baseflow/setup；Re 改变后 set_parameter('Re',Re)+setup。
时间模式二维曲线在 (beta_r,alpha_r) 平面，空间模式在 (beta_r,omega_r) 平面，文档与轴名应明确区分。
参考源码算法详细分析见 REFERENCE_ANALYSIS.md。

待处理风险：SLEPc adapter 未处理 nconv=0；重复调用 Vec 生命周期需要检查；已有 MPI 参数广播被注释；应明确单进程行为并避免新功能扩大既有缺陷。中性点收敛必须同时验证特征值残差和增长率容差，保留完整复数 seed 并限制搜索范围/步数。

## 原版隔离验证

复制原有 bin/START_LST、算例 config.ini 和 in/ 到 baseline_test/，避免原程序 cfg_loader 改写目标算例。
环境使用现有 PETSc 3.22.1、SLEPc 3.22.0、MPICH。
单点与连续两次 alpha 小步进均退出 0，未重现原源码注释中的主对角线零错误。
alpha=3.500: omega=-0.03933479054543185+0.01802792265614591i，绝对残差3.4827e-15。
alpha=3.501: omega=-0.03888732105721095+0.01802726681876940i，绝对残差5.4586e-15。
alpha=3.502: omega=-0.03843983939042378+0.01802635641223358i，绝对残差9.6023e-15。
原程序均有 MPICH yaksa 4 leaked handle pool objects 警告，已保留原始输出。
详见 commands.log、baseline_test/single.log、baseline_test/repeated_alpha.log。

## 当前状态

已完成 Git 安装包下载、完整原版归档、Git 原版分支与 bundle、源码算法分析、原版隔离运行和日志。
尚未执行远程上传，也未修改目标程序源码或添加功能。下一步需要用户提供有写权限的远程仓库地址，以遵守先上传原版再修改的顺序。

## 用户提供仓库后的继续执行

用户指定仓库：https://github.com/zhangyunjia0711/start_add_neutral_curve.git。
连接器可读取仓库，但写入请求连续出现传输错误；本机 Git 最初无凭据，执行 dry-run 时明确失败，未宣称上传成功。通过 Git Credential Manager 设备登录，经用户完成授权后，用 Git CLI 推送原版分支。
远程 refs/heads/backup/original-before-neutral-20261021 已核实为 7b12c3c37a47ee9caa664c6e08815c5dd976a335；original_push.log 和 original_remote_verification.log 保存证据。目标源码修改从此核验之后才开始。

创建 feature/neutral-curves 工作分支。直接切换原版分支因原有未提交/未跟踪文件被 Git 拒绝；在确认原版清单完全一致后，新建 feature 分支，用 read-tree 更新索引、update-ref 指向已经备份的原版提交，工作树文件不变。原 develop 分支和备份分支均保留。
为保留用户刚创建仓库的 README 和初始历史，拉取并合并 neutral-backup/main 的初始提交，无强制推送。

## 最终实现

- src/lst/module_neutral_numerics.f90：独立有界多根扫描、割线/二分校正、伪弧长续接、闭合/范围/有限性/残差检查。
- src/lst/module_neutral_curve.f90：调用 START 原有 LST 求解器，菜单31/32/33/36/37，独立 [neutral] 配置，时间/空间增长率，完整复数种子，换站/Re 更新与恢复，Tecplot 数据/日志输出。
- src/lst/module_lst.f90、main.f90：交互菜单、INI自动批处理、EOF安全退出、失败退出码、单进程批处理限制。
- src/lst/lst_eqn.f90：暴露所选特征值误差和可切换的目标跟踪状态。
- src/public/module_slepc_adapter.f90：零收敛返回NaN/空数组，旧Vec销毁，已填特征向量复用，绝对残差last_error；中性任务可选nearest-target，普通求解保留原增长率选择规则。
- tools/plot_neutral_curve.py：绘制根分支和增长率背景，保留失败/缺失区间，不跨区连接。
- doc/neutral_curve.md：中文说明；lst_test/lst_test/neutral_curve.ini：验证过的可运行示例。
- CMake及tests：数值核心与SLEPc adapter的CTest；Python绘图测试。

## 调试中发现并修正的问题

1. 首轮编译：global_parameters的size变量遮蔽Fortran size内建函数；以use关联重命名修复。
2. 扫描返回部分有效根及失败状态时，不能按根数覆盖为成功；现保留失败优先级，并检查峰值精化及交互态恢复。
3. 首次从配置点直接跳到扫描边界会跳模态；增加有界小步续接，首个微小步估计斜率，之后割线预测复特征值。
4. 实际算例揭示原adapter在多个收敛特征值之间按最大增长率切换，导致追踪错支；新增仅在中性任务期间开启的nearest-target选择，并测试其与默认规则不同。
5. 离散换站与连续波数变化混用斜率会错支；先在固定inner参数下逐个整数站位过渡，再沿连续参数预测校正，跨多个站位也逐站求解。
6. 独立测试中验证失败/NaN/过窄区间/无变号切触不虚报根；绘图测试修正了重复评估及非规则refine采样的遮罩处理。
以上失败、修正命令和运行摘要均在commands.log；最终各算例的配置、求解原始输出、场/中性点和终止状态位于real_cases/。

## 最终验证结果

- START_LST 编译成功，数值核心31项解析检查通过。
- CTest neutral_numerics、slepc_adapter_safety 两项均通过；adapter覆盖空结果、EPS/PEP重复求解、所选残差、默认与目标跟踪选择。
- Python绘图11项测试通过；真实中性曲线图已生成并目视检查。
- 14个实际场景通过：菜单31/32/33各跨3个站位（各6根）；37在3个beta上6根；36弧长13点并正常到达参数边界；空间37两根；交互任务结束后输出特征函数成功；Re扫描6根；跨站位步长2的扫描4根；错误轴、严格跳变限制、无根区域、MPI双进程分别按预期返回非零。
- 所有成功中性点均满足abs(growth)<=1e-7，真实时间模式37六根最大残差6.431e-9，0拒绝点；空间模式两根最大残差2.621e-8。
- beta=-3.2059、iloc=270的时间模式中性alpha约为2.986692664285及3.957146023035。
- 原版已存在的MPICH yaksa 4 leaked handle pool objects警告仍可能出现；新增adapter独立测试无Vec生命周期错误。未将此次验证视为物理网格收敛验证。

## 备份与交付范围

完整本地归档及依赖原工作树保存于此审计目录的original_full.tar.gz；GitHub原版分支保存START项目源码、原有本地修改及算例。stdlib的独立历史/本地修改另存本地bundle/patch及完整归档，不将其518MB构建/venv目录纳入项目Git提交。
新增功能在feature/neutral-curves上提交并同步到用户仓库main；原作者Gitee origin未被推送或更改。最终远程提交及推送核验结果将附于外部日志。

最终只读复核发现弧长零宽outer区间会被单行扫描padding意外扩大；现强制任务36要求outer_min<outer_max。新增zero_width_arc实际回归验证退出1、提示清晰且不生成越界曲线。累计14个实际场景通过。
