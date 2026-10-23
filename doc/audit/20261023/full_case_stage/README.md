# 75°工况阶段快照

这是未完成鼻尖/有限域验收时的版本检查点，不是完整认证结果。查看 outputs/本工况可变β计算说明.md 与 case_stage_status.json。

work/保存实际使用的脚本快照，依赖原工况目录、配置、H5及先前归档结果；不能把本目录当独立可运行示例。原完整计算目录和每文件SHA见archive_manifest.json。需要复现时保持原ROOT/work与ROOT/outputs结构和原CASE/M1/M2/M3输入，使用冻结v10（核心fbf78b1）；每个真实run的config.ini/solver.log/summary在原ROOT/runs。未将大型原始CFD或全部求解日志复制入Git。

不同数值/来源模式明确区分：实际中性点、正增长见证、共同共轭派生点、仅辅助剖面点、失败/unknown。此提交只归档证据与脚本，不修改Fortran核心。
