"""Finalize saved metadata without any solver invocation."""
import csv,datetime,hashlib,json,pathlib,subprocess
RUN=pathlib.Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/bluntplate/neutral_beta_full_20261022')
A=RUN/'outputs/complex_omega_interface_v10_20261024';C=RUN/'outputs/native10_beta_cuts_v10_20261023_2233'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def write(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
p=A/'complex_omega_probe.json';d=json.loads(p.read_text())
for x in d['points']:
    x['printed_beta_r']=x['beta_r'];x['beta_r']=d['reference_point']['beta_r']
write(p,d)
with (A/'complex_omega_points.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(d['points'][0]));w.writeheader();w.writerows(d['points'])
probe=A/'probe_complex_omega.py';text=probe.read_text(encoding='utf-8-sig')
text=text.replace("beta_r=row[9],beta_i=row[10],","beta_r=point['beta_r'],printed_beta_r=row[9],beta_i=row[10],")
probe.write_text(text)
S=RUN/'work/solver_source';original=pathlib.Path('/home/zaza/start/START_ALL_IN_ONE/START_ALL_IN_ONE')
specs=[(S/'src/lst/cfgio_adapter.f90',263,285),(S/'src/lst/module_lst.f90',113,117),
    (S/'src/lst/module_lst.f90',136,155),(S/'src/lst/module_lst.f90',192,226),
    (S/'src/lst/module_neutral_curve.f90',70,95),(S/'src/lst/module_neutral_curve.f90',378,402),
    (S/'src/lst/lst_eqn.f90',923,933)]
originals=[]
for rel in ['src/lst/cfgio_adapter.f90','src/lst/module_lst.f90','src/lst/lst_eqn.f90']:
    p=original/rel;z=p.read_text().splitlines()
    terms=('omega','ome =','case (7)','pk*ZI','complex(REAL64)','ome*G_0')
    originals.append(dict(path=str(p),sha256=sha(p),matches=[dict(line=i+1,text=l) for i,l in enumerate(z) if any(t in l for t in terms)]))
write(A/'source_original_snippets.json',originals)
snippets=[dict(path=str(p),sha256=sha(p),start_line=lo,end_line=hi,
    lines=[dict(line=i+1,text=l) for i,l in enumerate(p.read_text().splitlines()) if lo<=i+1<=hi]) for p,lo,hi in specs]
write(A/'source_interface_snippets.json',snippets)
build=json.loads((RUN/'outputs/v10_build_validation.json').read_text())
mismatches=[str(S/rel) for rel,h in build['source_hashes'].items() if sha(S/rel)!=h]
assert not mismatches,mismatches
git_status=subprocess.check_output(['git','-C',str(S),'status','--short'],text=True)
write(A/'final_source_integrity.json',dict(created=now,
    all_v10_build_source_hashes_match=True,compared_files=len(build['source_hashes']),mismatches=mismatches,
    solver_sha256=sha(RUN/'outputs/bin/START_LST_wall_local_v10'),
    git_status=git_status,qualification='Only pre-existing parent audit-cleanup deletions in doc/audit are present; this task made no core edits.'))
readme="""# v10 复频率接口核验

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
"""
(A/'README.md').write_text(readme)
report=json.loads((C/'native10_cut_report.json').read_text())
write(C/'final_source_integrity.json',json.loads((A/'final_source_integrity.json').read_text()))
for area,names in [(A,['probe_complex_omega.py','finalize_metadata.py','complex_omega_probe.json','complex_omega_points.csv','complex_omega_integrity.json','final_source_integrity.json','source_interface_snippets.json','source_original_snippets.json','README.md']),
    (C,['reduce_cuts.py','native10_cut_report.json','native10_cut_samples.csv','native10_cut_summary.csv','native10_cut_witnesses.csv','native10_beta_cuts.png','reduction_integrity.json','final_source_integrity.json','README.md','cut_progress.json','cut_plan.json','source_fresh_join.json','integrity.json','run_cuts.py'])]:
    write(area/('manifest.json' if area==A else 'reduction_manifest.json'),dict(created=now,files=[dict(path=str(area/n),sha256=sha(area/n)) for n in names],
        raw_bindings=report['raw_bindings'] if area==C else d['points']))
print(json.dumps(dict(created=now,sourcefiles_compared=len(build['source_hashes']),core_unchanged=True,
    cut_report_sha256=sha(C/'native10_cut_report.json'),cut_manifest_sha256=sha(C/'reduction_manifest.json'),
    complex_report_sha256=sha(A/'complex_omega_probe.json'),complex_manifest_sha256=sha(A/'manifest.json'))))

