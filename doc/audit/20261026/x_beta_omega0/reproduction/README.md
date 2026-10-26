# 复现说明

实际计算目录：`/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/neutral_xbeta_omega0_20261026`。
该目录保存两工况输入副本、全部原矩阵、特征向量和原程序输出；原输入路径及SHA见input_manifest.json。冻结程序SHA为634b45c5d30e076f5b86ccb136e0734c1905beedde108477dceb8a75adc3853c。

环境：Ubuntu WSL，Python3.8.10，NumPy1.24.4，SciPy1.10.1，h5py3.11.0，Matplotlib3.7.5；原程序使用/usr/local/mpich-4.1.2/lib；每个数值进程的OPENBLAS_NUM_THREADS、OMP_NUM_THREADS、MKL_NUM_THREADS均设为1。

现有计算可在原目录直接运行 `python3 work/audit.py`、`python3 work/positive_checks.py`，读取缓存并复核原输入；`python3 work/tecplot.py --verify` 需要原生回读文件，若已清理可在Windows输出目录重跑tecplot_roundtrip.mcr生成。

完整重算应使用新目录，禁止删除或覆盖现有结果：先在bootstrap.py中将ROOT末级目录名改为新名称，运行bootstrap.py；它会复制并校验原输入，同时生成common.py和run_adapter.py。把这里其他Python脚本复制到新ROOT/work。publish.py、tecplot.py、report.py中的DEST也应改为新的Windows交付目录。

按依赖顺序运行：solver.py（种子/β扫描）、prepare.py（1202站原矩阵）、contour.py bluntplate、contour.py rough、finalize.py bluntplate、finalize.py rough、audit.py、positive_checks.py、publish.py。两个contour以及两个finalize可在各自依赖完成后独立并行。

最后运行tecplot.py生成宏；在Windows使用native_tecplot.ps1（从75°任务目录执行，若重算需调整目标目录），它调用Preplot转换原生.plt并用Tecplot进行回读和布局保存。再运行tecplot.py --verify及report.py。宏使用相对ASCII文件名，避免目录名75°经ASCII转码导致路径错误。

contour.py --resume 只用于日志状态未完成的轮廓。finalize.py枚举全部实际站位穿越，不依赖预测步恰好落在网格上；每个穿越都通过原空间多项式求根。point_kind和conjugate含义见结果说明。最终原生站位必须各有两侧边界点；任何模式MAC、残差、稳定侧符号或完整性断言失败都应保留记录并诊断。

大体积输入、矩阵和向量保留在WSL，不随Git提交重复上传；这里提供最终结果、哈希、计算证据和复现代码。本次不是CFD网格加密，也不是全谱支枚举。
