# 空间根方向与复alpha鞍点：有限核验

2026-10-24T12:05:08.880135+08:00

本次从10月24日恢复时刻建立新日志，未把quota/暂停期间写成持续计算。沿用冻结v10与原H5，复用已验空间PEP三系数矩阵，不修改core、边界条件或物理模型。两个站分别固定各自实际beta，只沿同一temporal简单支有界求复alpha下omega_alpha=0。

|原生站|固定beta|复alpha0|复omega0|迭代|
|---|---:|---:|---:|---:|
|10|12.137900983290656|−22.055166411716655 +0.06120268031775814i|8.902789246227368 −1.02547876047e−5i|3|
|11|12.34948078116693|−22.811316916796155 −0.13839233297358722i|9.045241670115805 −4.63671472950e−5i|4|

两处左右特征向量的同PEP回代残差约3–4e−13，原生来源q/MAC保持，omega_alpha绝对值分别2.86e−10与5.82e−14。native10的omega_alphalpha≈0.002306214−0.001224215i，native11≈0.002107695−0.001279647i；实/虚alpha方向h/h2的导数核验一致。它们是局部鞍点候选，不是已完成因果分类。

另实际运行8次原v10 complex-omega空间PEP，保留append原生Vec并逐index绑定打印高精度alpha。全部通过原残差1e−7与MAC.95门限；最大原PEP残差2.24473e−8、最大同矩阵回代残差2.24473e−8、最小对鞍点raw MAC0.999991737。这里不把最差残差说成机器精度。8次原生短调用合计4.157秒，另有单线程稀疏矩阵特征求解。

native10在omega0+i*epsilon、epsilon=1e−4/1e−5/1e−6处分别得到两不同的真实合格空间根；间距为0.553500705、0.175040747、0.055353006，除sqrt(epsilon)后为55.3500705、55.3527445、55.3530060，呈平方根收拢。两根中点到鞍点的误差约1.73e−3、1.73e−4、1.74e−5。native11也用两个原生邻根复核。上述是局部空间双支/二阶退化的数值证据，不是精确代数重数证明，更不是pinch两半平面起源证明。

资料已核对：Huerre与Monkewitz的正确DOI为[10.1146/annurev.fl.22.010190.002353](https://www.annualreviews.org/content/journals/10.1146/annurev.fl.22.010190.002353)，不是此前的002151。实际逐页读了[原文扫描件](https://indico.ictp.it/event/a02272/contribution/25/material/0/0.pdf)的原印刷pp477–485，重点pp480–484、图2、式9–17；其封面是1994 ICTP讲义，后附1990原文复印。出版社与扫描件页码上限存在不一致，未假称已解决。Briggs1964仅核实[出版社信息](https://mitpress.mit.edu/9780262523509/electron-stream-interaction-with-plasmas/)，全文未获取。Bers1983原章未获取；补读了[Ram与Bers1991原作者报告的文字镜像](https://studylib.net/doc/11312294/pfc-ja-91-4-absolute--versus-space--plasmas)第II节pp2–5，明确注明镜像来源。

Huerre–Monkewitz pp483–484强调：omega_alpha=0及omega0的符号单独不够，必须核对空间两支在足够高的因果Laplace轮廓上来自Fourier轮廓的两侧；同半平面支的合并不构成该pinch。Ram–Bers1991第II节pp3–5同样要求随Laplace轮廓下降作解析延拓，并指出一般传播速度不能由任意模态群速代替。因此Re(omega_alpha)的局部符号、MAC连续和大空间g均不能自行决定下游因果支。

本次没有严格的高Imomega上界。要完成因果检验，需证明初始Laplace线高于固定beta、全部实alpha相关temporal谱/奇点，并从那里跟踪两支及允许的Fourier轮廓变形直到夹逼，同时排除其他控制奇点。少数target temporal根、实omega=40或beta到±40都不是这个上界。当前M含约束且非Hermitian，不能套用正定Rayleigh商虚构界。

因此两个负Imomega0只能作为已计算局部鞍点的性质；在起源、主导pinch与其他谱支未核验前，整体local算子的absolute/convective类别仍未知，更不能由此给全鼻部物理分类、下游不稳定认定或无限包络结论。现有鼻部候选继续保留物理未知及parallel/zero-curvature模型限制。本阶段按有界任务到此结束，不继续无穷高Imomega探测。

保留脚本、原系数矩阵来源、配置、append原生向量、PEP日志、鞍点JSON、原文PDF与实际查看页图。已删除本任务独占的临时pypdf vendor及无文本提取文件，删除前验证绝对路径在本RUN/work目录。原H5/bin及只读9个core源hash均保持。
