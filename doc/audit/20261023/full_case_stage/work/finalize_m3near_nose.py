"""Preserve bounded-jump nose recovery evidence and its genuine MAC stop."""
import json,re
from pathlib import Path
import full_driver as d
from axial_wall_local import stamp,NUMBER

R=d.ROOT
out=R/'outputs'
state=json.loads((out/'axial_wall_local_v10_near_j64.json').read_text())
family=state['families']['M3near']
names={b['name'] for b in family['blocks'] if re.match(r'axial_local_j(?:8|16|32|64)\.0_',b['name'])}
probes=[]
for key in ['37_36','24_23','17_16','10_9']:
    path=out/f'local_mode_audit_M3near_nose_probe_{key}.json'
    raw=json.loads(path.read_text())
    assert raw['all_requested_probes_passed_same_root']
    probes.append(dict(path=str(path),sha256=d.sha(path),result=raw))
    for p in raw['probes']:
        for k in ['forward_block','back_block']:names.add(p[k]['name'])
last=family['blocks'][-1]
log=(R/'runs'/last['name']/'solver.log').read_text()
rx=r'MAC best index=\s*(\d+)\s+eigen=\s*\(\s*('+NUMBER+r')\s*,\s*('+NUMBER+r')\s*\)\s+overlap=\s*('+NUMBER+r')\s+residual=\s*('+NUMBER+r')\s+jump=\s*('+NUMBER+r')'
all_mac_best=[dict(index=int(m[1]),alpha_r=float(m[2]),alpha_i=float(m[3]),MAC=float(m[4]),
               residual=float(m[5]),jump=float(m[6])) for m in re.finditer(rx,log)]
rejected=[v for v in all_mac_best if v['MAC']<.95]
assert rejected and max(v['MAC'] for v in rejected)<.95
assert all(v['residual']<=1e-7 and v['jump']<64 for v in rejected)
summaries=[json.loads((R/'runs'/n/'summary.json').read_text()) for n in sorted(names)]
records=[r for r in state['records'] if r['family']=='M3near' and r['status']==0]
new=[r for r in records if r['iloc']<37]
original=out/'axial_wall_local_v10_near.json'
original_hash=json.loads((out/'axial_wall_local_v10_near_j8.json').read_text())['inherited_seed_file_sha256']
assert d.sha(original)==original_hash
assert len(records)==838 and {r['iloc'] for r in records}==set(range(4,842))
assert len(new)==33 and all(r['printed_pep_error']<=1e-7 for r in new)
report=dict(created=stamp(),status='stopped_on_actual_MAC_gate',source=str(out/'axial_wall_local_v10_near_j64.json'),
    source_sha256=d.sha(out/'axial_wall_local_v10_near_j64.json'),solver=state['solver'],solver_sha256=state['solver_sha256'],
    source_original_preserved=d.sha(original)==original_hash,accepted_stations=838,range=[4,841],new_accepted_stations=33,
    unknown_stations=[1,2,3],fixed_omega=0,fixed_beta=-7,final_anchor=family['directions']['up']['anchor'],
    stopping_block=last,rejected_MAC_candidates=rejected,probe_results=probes,
    actual_new_calls=len(names),actual_new_solves=sum(s['solves'] for s in summaries),
    actual_new_seconds=sum(s['seconds'] for s in summaries),
    maximum_new_accepted_PEP_residual=max(r['printed_pep_error'] for r in new),
    qualification='Unclassified same-origin continuation; no physical identity or all-beta stability claim. No i1 seed exported.',
    call_summary_paths=[str(R/'runs'/n/'summary.json') for n in sorted(names)])
d.snapshot()
target=out/'local_mode_audit_M3near_nose_recovery.json'
target.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
rows=[]
for p in probes:
    for v in p['result']['probes']:
        rows.append('|%s|%g|%.12g|%.12g|%.5g|'%(Path(p['path']).stem.split('probe_')[-1],v['jump_limit'],v['delta_alpha'],v['min_actual_mac'],v['roundtrip_alpha_difference']))
md='\n\n'.join(['# M3near 鼻尖固定参数恢复结果',stamp(),
    '固定 ω=0、β=−7，实际接受范围由 i37..841 扩展为 i4..841，共838站。i1..3仍未知。最终独立输出 `'+report['source']+'`；旧原输出和各级j8/j16/j32输出均保留，旧文件hash未变。',
    '相邻根变化增大曾在4、8、16、32绝对jump门限处停止。每次升门限前均独立正向求解并返回原站，保留实际MAC/原PEP残差、不同门限同根比较；未放宽MAC0.95或残差1e−7。',
    '\n'.join(['|相邻步|门限|实际复α差|最小实际MAC|返回α误差|','|---|---:|---:|---:|---:|']+rows),
    '最终i4实际 α=230.9455858284939+448.9590432636435i，x=−34.99325341687227 mm。i4→3最佳候选 α≈220.32083598369+468.24718835501i，原多项式残差≈2.46e−12，实际复α差≈22.02085<64；但最佳MAC≈0.940228075821<0.95。nev1、retry8及多个shift均未通过实际MAC。因此此处是向量归属门限失败，不能继续把相近特征值标为同族接受点。未扩至更大门限，未导出未经证实的i1种子。',
    '所有新增实际调用共%d次、原始EVP solves%d、累计求解墙时%.3f秒。33个新增接受站的最大匹配原多项式残差为%.5g。完整逐调用来源、原始日志、向量检查、门限探针与失败候选见同名JSON。'%(report['actual_new_calls'],report['actual_new_solves'],report['actual_new_seconds'],report['maximum_new_accepted_PEP_residual']),
    '本报告只整理实际操作与求解证据，未新增理论性结论。M3near是搜索起点/来源标签；上述模态物理身份仍未分类。固定β上的负增长不能推出可变β或全谱稳定性。局部切平面模型近鼻尖的适用性也不能仅凭残差和MAC证明。',
    '结束审查：授权的有界门限恢复与失效定位完成，i1目标未达到，原因与unknown区已明确保留。冻结v10与原配置/H5的snapshot核验通过，旧输出hash保持。保留必要复现脚本、各级结果、配置、原始日志和探针；没有归属明确而不再需要的尝试代码，未删除用户或其他任务文件。'])+'\n'
(out/'local_mode_audit_M3near_nose_recovery_说明.md').write_text(md)
with (out/'工作日志_2026-10-23_M3near鼻尖固定参数恢复.md').open('a') as f:
    f.write('\n### '+stamp()+' 有界恢复结束审查\n\n')
    f.write('- 补录：j32接受i16..10；独立64探针10→9及回10通过，|δα|35.20802305741382，MAC0.9803889491764841，回程误差3.6270652172114494e−9。执行后原日志追加命令因末尾字面量\\n发生SyntaxError；该记录命令未写入，本段如实补录，不影响求解。\n')
    f.write('- j64实际再接受i9..4，最终838站；i3真实最佳MAC0.940228075821低于0.95，jump22.02085<64且PEP误差通过。按约定停止；i1..3未知，未导出i1种子。\n')
    f.write('- 新实际调用%d次 / solves%d / %.3f秒，新增33站。详见独立恢复报告及JSON。旧输出、冻结binary、原输入hash核验通过；阶段目标是有界核验，i1目标未达但失败证据完整。无可删除的不必要尝试代码。\n'%(report['actual_new_calls'],report['actual_new_solves'],report['actual_new_seconds']))
print(json.dumps({k:report[k] for k in ['status','accepted_stations','range','actual_new_calls','actual_new_solves','actual_new_seconds']}))
