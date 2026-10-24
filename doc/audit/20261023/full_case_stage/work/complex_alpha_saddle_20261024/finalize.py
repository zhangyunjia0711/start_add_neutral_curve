"""Qualification, primary-source record, and bounded saddle-stage final audit."""
import hashlib,json,shutil
from pathlib import Path
import saddle_probe as sa
tp=sa.tp;RUN=sa.RUN;HERE=sa.HERE
OUT=RUN/'outputs/spatial_root_causality_saddle_review_20261024.json'
MD=OUT.with_suffix('.md')

def main():
    assert not OUT.exists()
    inputs=[RUN/'outputs'/name for name in ['complex_alpha_saddle_native10_20261024.json','complex_alpha_saddle_native10_neighbors_20261024.json','complex_alpha_saddle_native11_20261024.json']]
    ten,pairs,eleven=[json.loads(p.read_text()) for p in inputs]
    rootpairs=[p for r in pairs['records'] for p in r['pairs']]+eleven['native_neighbor_pairs']
    original_hashes=RUN/'outputs/temporal_frozen_global_source_hashes_20261023.json';hashes=json.loads(original_hashes.read_text())
    assert all(tp.driver.sha(Path(p['path']))==p['sha256'] for p in hashes['files'])
    tp.driver.snapshot();assert tp.driver.sha(tp.BIN)==tp.SHA
    ref=HERE/'references/Huerre_Monkewitz_1990_ICTP.pdf'
    bibliography=[dict(authors='P. Huerre and P. A. Monkewitz',year=1990,title='Local and Global Instabilities in Spatially Developing Flows',
        correct_doi='10.1146/annurev.fl.22.010190.002353',previously_supplied_doi='10.1146/annurev.fl.22.010190.002151',
        publisher='https://www.annualreviews.org/content/journals/10.1146/annurev.fl.22.010190.002353',
        original_text_scan='https://indico.ictp.it/event/a02272/contribution/25/material/0/0.pdf',
        actual_read='Original article pp477–485, especially pp480–484 Eq9–17 and Fig2, as scanned reproduction after ICTP March1994 workshop cover; visually inspected PDF pages4–8',
        local_pdf=str(ref),local_sha256=tp.driver.sha(ref),
        limit='Publisher metadata currently prints pp473–537; the scan first article page prints473–517. No claim made that this pagination discrepancy or erratum has been resolved.'),
        dict(authors='A. K. Ram and A. Bers',year=1991,title='Absolute versus convective analysis of instabilities in space plasmas',report='MIT PFC/JA-91-4',
        original_author_text_mirror='https://studylib.net/doc/11312294/pfc-ja-91-4-absolute--versus-space--plasmas',
        actual_read='Original-author report text mirrored online, SectionII pp2–5 Eq4–6: Fourier/Laplace contour deformation and opposite-side branch origin; p5 distinguishes propagation from generic modal group velocity',
        limit='Mirror original text was read; no institutional PDF of this 1991 report or full Bers1983 chapter obtained in this session.'),
        dict(authors='R. J. Briggs',year=1964,title='Electron-Stream Interaction with Plasmas',
        publisher='https://mitpress.mit.edu/9780262523509/electron-stream-interaction-with-plasmas/',actual_read='Publisher metadata and description only',
        limit='Full book not retrieved or read; its technical criterion is not presented as independently checked from that book.'),
        dict(authors='A. K. Ram and A. Bers',year=1992,title='Comments on absolute and convective instabilities',doi='10.1029/91GL02953',
        publisher='https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/91GL02953',
        institutional_record='https://dspace.mit.edu/entities/publication/d9a8e3c2-7761-4a14-935f-ca4b3552ff0e',
        actual_read='Metadata/abstract only; attempted original PDF failed with429/403/timeout',limit='Not used for unverified full-text technical claims.'),
        dict(authors='A. Bers',year=1983,title='Space-time evolution of plasma instabilities—absolute and convective',
        citation_source='Ram/Bers1991 and Huerre/Monkewitz1990 references',actual_read='Full original chapter not obtained',
        limit='Bers2016 publisher abstract was accessible but is a later book chapter, not substituted for the unread1983 original.')]
    requirements=[
        'At fixed real beta, establish the relevant temporal dispersion branches and all singularities of the initial-value resolvent, including constraints and any continuous spectrum/truncation issues.',
        'Choose a causal Laplace line Imomega=gamma above ALL relevant temporal singularities for real alpha, with a proved finite bound or a justified limiting argument; a few k=2 target eigenvalues do not establish this.',
        'Track the two spatial branches from that causal line while lowering/deforming the omega contour; verify origins on opposite sides of the real-alpha Fourier contour and preserve admissible deformation without unaccounted singularities.',
        'Demonstrate that those particular branches pinch the Fourier contour at the candidate saddle; same-half-plane coalescence or a local square-root pair is insufficient.',
        'For classification of the local operator, identify the controlling admissible pinch frequency, not merely one finite candidate; global/physical flow conclusions require additional model and boundary-condition justification.'
    ]
    cleanup=[]
    vendor=(HERE/'vendor').resolve();expected=(RUN/'work/complex_alpha_saddle_20261024/vendor').resolve()
    assert vendor==expected and str(vendor).startswith(str(RUN.resolve())+'/work/')
    if vendor.exists():
        shutil.rmtree(vendor);cleanup.append(dict(path=str(vendor),reason='Temporary pypdf5.5.0/typing_extensions4.13.2 installed only to extract scanned PDF page images; bundled Windows pypdf is available and numeric reproduction does not import these packages'))
    empty=HERE/'references/Huerre_Monkewitz_1990_ICTP.txt'
    if empty.exists():
        assert empty.resolve().parent==(HERE/'references').resolve()
        empty.unlink();cleanup.append(dict(path=str(empty),reason='pypdf extraction contained page labels only: PDF is scanned; no usable text'))
    report=dict(created=tp.stamp(),status='bounded_local_saddle_and_native_neighbor_validation_complete',
        source_paths=[dict(path=str(p),sha256=tp.driver.sha(p)) for p in inputs],bibliography=bibliography,
        model='Frozen v10 local parallel wall-coordinate spatial PEP, original CFD H5 and boundary conditions; fixed real beta for each station',
        solver=str(tp.BIN),solver_sha256=tp.SHA,cases=[dict(station=10,point=ten['result']['point'],iterations=ten['result']['iterations'],h_h2=ten['saddle_h_h2_verification']),
            dict(station=11,point=eleven['result']['point'],iterations=eleven['result']['iterations'],h_h2=eleven['h_h2_verification'])],
        local_neighbor_result=dict(actual_native_calls=len(rootpairs),actual_native_seconds=sum(p['seconds'] for p in rootpairs),
            largest_reported_original_PEP_residual=max(p['reported_polynomial_residual'] for p in rootpairs),
            largest_same_matrix_residual=max(p['exported_same_PEP_matrix_residual'] for p in rootpairs),minimum_MAC_to_saddle=min(p['source_MAC'] for p in rootpairs),
            original_residual_gate=1e-7,original_MAC_gate=.95,
            native10_separation_over_sqrt_epsilon=[r['separation_over_sqrt_epsilon'] for r in pairs['records']],
            qualification='All eight pairs pass the original1e-7 PEP gate; worst2.24e-8 is explicitly retained, not presented as machine precision. Spatial root conditioning worsens near coalescence.'),
        causal_pinch_origin_verified=False,absolute_or_convective_classification='unknown',
        required_causal_validation=requirements,proved_high_Imomega_upper_bound=None,
        upper_bound_qualification='No rigorous gamma above all temporal singularities is available. Omega40 and beta±40 are sampled real parameters, not an Imomega spectral bound. A singular non-Hermitian M cannot be treated as positive definite to invent an energy or Rayleigh-quotient bound.',
        physical_downstream_instability_inferred=False,
        local_observation='Two finite simple-temporal-branch saddles with negative Imomega0 and genuinely residual-qualified nearby spatial roots. This does not establish causal branch selection, leading pinch, physical convective/absolute instability, or an infinite neutral envelope.',
        integrity=dict(original_inputs_and_binary_unchanged=True,read_only_core_manifest=str(original_hashes),read_only_core_hashes_unchanged=True),
        retained_reproduction=[dict(path=str(p),sha256=tp.driver.sha(p)) for p in HERE.glob('*.py')],cleanup=cleanup,
        failures=['Supplied Huerre-Monkewitz DOI002151 could not be verified; corrected against publisher002353.',
            'WSL pdftotext absent; bundled fitz absent. Temporary pypdf extraction found no text; visual page-image inspection used instead.',
            'Briggs full book and Bers1983 full chapter not obtained; MIT1991 comments PDF attempts returned429/403/timeouts.',
            'No high-Imomega branch-origin continuation or all-spectrum bound attempted under the bounded task; no causal classification inferred.'])
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    text=['# 空间根方向与复alpha鞍点：有限核验',tp.stamp(),
        '本次从10月24日恢复时刻建立新日志，未把quota/暂停期间写成持续计算。沿用冻结v10与原H5，复用已验空间PEP三系数矩阵，不修改core、边界条件或物理模型。两个站分别固定各自实际beta，只沿同一temporal简单支有界求复alpha下omega_alpha=0。',
        '|原生站|固定beta|复alpha0|复omega0|迭代|',
        '|---|---:|---:|---:|---:|',
        '|10|12.137900983290656|−22.055166411716655 +0.06120268031775814i|8.902789246227368 −1.02547876047e−5i|3|',
        '|11|12.34948078116693|−22.811316916796155 −0.13839233297358722i|9.045241670115805 −4.63671472950e−5i|4|',
        '两处左右特征向量的同PEP回代残差约3–4e−13，原生来源q/MAC保持，omega_alpha绝对值分别2.86e−10与5.82e−14。native10的omega_alphalpha≈0.002306214−0.001224215i，native11≈0.002107695−0.001279647i；实/虚alpha方向h/h2的导数核验一致。它们是局部鞍点候选，不是已完成因果分类。',
        '另实际运行8次原v10 complex-omega空间PEP，保留append原生Vec并逐index绑定打印高精度alpha。全部通过原残差1e−7与MAC.95门限；最大原PEP残差2.24473e−8、最大同矩阵回代残差2.24473e−8、最小对鞍点raw MAC0.999991737。这里不把最差残差说成机器精度。8次原生短调用合计4.157秒，另有单线程稀疏矩阵特征求解。',
        'native10在omega0+i*epsilon、epsilon=1e−4/1e−5/1e−6处分别得到两不同的真实合格空间根；间距为0.553500705、0.175040747、0.055353006，除sqrt(epsilon)后为55.3500705、55.3527445、55.3530060，呈平方根收拢。两根中点到鞍点的误差约1.73e−3、1.73e−4、1.74e−5。native11也用两个原生邻根复核。上述是局部空间双支/二阶退化的数值证据，不是精确代数重数证明，更不是pinch两半平面起源证明。',
        '资料已核对：Huerre与Monkewitz的正确DOI为[10.1146/annurev.fl.22.010190.002353](https://www.annualreviews.org/content/journals/10.1146/annurev.fl.22.010190.002353)，不是此前的002151。实际逐页读了[原文扫描件](https://indico.ictp.it/event/a02272/contribution/25/material/0/0.pdf)的原印刷pp477–485，重点pp480–484、图2、式9–17；其封面是1994 ICTP讲义，后附1990原文复印。出版社与扫描件页码上限存在不一致，未假称已解决。Briggs1964仅核实[出版社信息](https://mitpress.mit.edu/9780262523509/electron-stream-interaction-with-plasmas/)，全文未获取。Bers1983原章未获取；补读了[Ram与Bers1991原作者报告的文字镜像](https://studylib.net/doc/11312294/pfc-ja-91-4-absolute--versus-space--plasmas)第II节pp2–5，明确注明镜像来源。',
        'Huerre–Monkewitz pp483–484强调：omega_alpha=0及omega0的符号单独不够，必须核对空间两支在足够高的因果Laplace轮廓上来自Fourier轮廓的两侧；同半平面支的合并不构成该pinch。Ram–Bers1991第II节pp3–5同样要求随Laplace轮廓下降作解析延拓，并指出一般传播速度不能由任意模态群速代替。因此Re(omega_alpha)的局部符号、MAC连续和大空间g均不能自行决定下游因果支。',
        '本次没有严格的高Imomega上界。要完成因果检验，需证明初始Laplace线高于固定beta、全部实alpha相关temporal谱/奇点，并从那里跟踪两支及允许的Fourier轮廓变形直到夹逼，同时排除其他控制奇点。少数target temporal根、实omega=40或beta到±40都不是这个上界。当前M含约束且非Hermitian，不能套用正定Rayleigh商虚构界。',
        '因此两个负Imomega0只能作为已计算局部鞍点的性质；在起源、主导pinch与其他谱支未核验前，整体local算子的absolute/convective类别仍未知，更不能由此给全鼻部物理分类、下游不稳定认定或无限包络结论。现有鼻部候选继续保留物理未知及parallel/zero-curvature模型限制。本阶段按有界任务到此结束，不继续无穷高Imomega探测。',
        '保留脚本、原系数矩阵来源、配置、append原生向量、PEP日志、鞍点JSON、原文PDF与实际查看页图。已删除本任务独占的临时pypdf vendor及无文本提取文件，删除前验证绝对路径在本RUN/work目录。原H5/bin及只读9个core源hash均保持。']
    blocks=[];rows=[]
    for line in text:
        if line.startswith('|'):rows.append(line);continue
        if rows:blocks.append('\n'.join(rows));rows=[]
        blocks.append(line)
    if rows:blocks.append('\n'.join(rows))
    MD.write_text('\n\n'.join(blocks)+'\n')
    sa.log('原文/依赖操作补录：WSL无pdftotext，bundled无fitz；临时pypdf5.5.0得到扫描页图，正文无文本层。Bers/Briggs完整原文获取失败、HTTP429/403/timeout都保留为未获取，不作已读声明。已读HM原印刷pp477–485和Ram/Bers1991报告镜像II节。')
    sa.log('结束审查：两局部saddle、8原生邻根有界验证完成；全部来源/core/H5/bin hash保持。无严格高Imomega界、无因果pinch起源，absolute/convective类别unknown。清理 '+json.dumps(cleanup)+'；保留必要复现证据与文献。主报告可引用 '+str(MD))
    print(str(OUT));print(str(MD))

if __name__=='__main__':main()
