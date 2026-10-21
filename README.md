# start_add_neutral_curve
20261021


已添加中性曲线菜单 **31 / 32 / 33 / 36 / 37**，以及 INI 批处理与 Python 绘图。

- [中文使用说明](doc/neutral_curve.md)
- [完整配置示例](lst_test/lst_test/neutral_curve.ini)
- [过程日志和验证记录](doc/audit/20261021/PROCESS.md)
- [真实算例曲线预览](doc/audit/20261021/real_cases/task37/neutral_preview.png)
- 原版备份分支：`backup/original-before-neutral-20261021`（提交 `7b12c3c`）。

在已配置 PETSc / SLEPc / Fortran stdlib 的 WSL 环境中：

```bash
cmake -S . -B build-neutral
cmake --build build-neutral --target START_LST test_neutral_numerics test_slepc_adapter -j4
ctest --test-dir build-neutral --output-on-failure
cd lst_test/lst_test
../../bin/START_LST -f neutral_curve.ini
python3 ../../tools/plot_neutral_curve.py out/neutral_roots.plt \
  --field out/neutral_field.plt --x-label beta --y-label alpha \
  --output out/neutral.png
```

当前中性曲线模块支持单 MPI 进程的正问题 LST，时间/空间模式均已验证。
