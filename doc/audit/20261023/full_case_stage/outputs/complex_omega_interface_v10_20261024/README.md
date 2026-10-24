# v10 复频率接口核验

未改核心。原交互入口一次进程共 3 次空间 PEP 求解，输入位于 runs/interactive_omega_i_three_points/stdin.txt：
7 1 0.005、7 1 0.005、-1 0 0。

native10、β=12.13790098329066，从已验 Ω=40 实轴点，沿 Ω_i=0 →0.005→0.01。
3 个选中本征对都有同求解组日志、原始 Vec 和残差绑定。每步 raw MAC ≥0.99999998279，最大残差 4.21e-12。实际耗时约 1.013 s。

接口事实（附 source_interface_snippets.json 与原始源码匹配片段）：

- cfgio 读两个 omega 分量；module_lst 构造 omega(1)+i*omega(2)。原交互 k=7 改变虚部，空间矩阵保留 complex omega。
- neutral 入口在进入扫描前拒绝非实 omega/beta（阈值 1e-12）。扫描 omega 时另会将虚部设 0。
- full_driver.py 也把 seed omega 写成 real,0；不能原样用于复频率。
- 此探针走原交互默认选根策略，外部再以真实向量 MAC/残差核验。当前 neutral 未提供以 omega_i 为轴的 MAC 自适应续算；可在不改 core 的独立驱动中读导出本征对、限制步长并接受/拒绝。
- 3 个邻点只验证接口和局部向量连续性，不验证全局解析支拓扑、pinching、物理 absolute 或下游因果传播。

complex_omega_probe.json 保存完整参数、残差、Vec 记录索引、scalar join 和 SHA；complex_omega_points.csv 为简表。
全部 v10 构建源文件 SHA 与冻结构建清单一致；原始输入、原程序与冻结 v10 SHA 不变。
probe_complex_omega.py 保留复现，已有结果目录时主动拒绝覆盖；本次只执行一次进程，不含扩频或全域扫描。
