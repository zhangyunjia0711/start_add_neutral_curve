"""Read-only reduction of the six already finished native10 finite-beta cuts."""
import csv, datetime, hashlib, json, pathlib, struct, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
AREA=pathlib.Path(__file__).resolve().parent
RUN=AREA.parent.parent
COMP=RUN/'outputs/upper_boundary_139_v10_20261023/compensation_20261023_2148'
sys.path.insert(0,str(COMP)); import compensate as c
def sha(p): return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def stamp(): return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
def save(name,data): (AREA/name).write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def csvsave(name,rows):
    cols=list(dict.fromkeys(k for row in rows for k in row))
    with (AREA/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(rows)
def arrays_and_pairs(path):
    data=(path/'out/eigenvectors.petsc').read_bytes(); qs=[]; off=0
    while off<len(data):
        cls,n=struct.unpack_from('>ii',data,off); off+=8
        assert cls==1211214 and 0<n<1000000
        end=off+16*n; assert end<=len(data)
        qs.append(np.frombuffer(data[off:end],dtype='>c16').astype(np.complex128));off=end
    pairs=[];groups=[]
    for line in (path/'solver.log').read_text(errors='replace').splitlines():
        if 'ncov=' in line:groups.append([])
        m=c.PAIR.search(line)
        if m:
            p=dict(index=int(m[1]),alpha_r=float(m[2].replace('D','E')),
                   alpha_i=float(m[3].replace('D','E')),residual=float(m[4].replace('D','E')))
            assert groups; groups[-1].append(p);pairs.append(p)
    assert len(qs)==len(pairs)>0
    return qs,pairs,[len(x) for x in groups]
progress=json.loads((AREA/'cut_progress.json').read_text())
assert progress['finished'] and progress['pending_tasks']==0 and len(progress['full_results'])==6
rows=[];summary=[];witness=[];bindings=[];curves={}
for item in progress['full_results']:
    file=pathlib.Path(item['file']);d=json.loads(file.read_text());r=d['run']
    assert r['complete'] and d['source_identity_verified']
    p=AREA/'runs'/r['name'];qs,pairs,groups=arrays_and_pairs(p)
    rawsha=sha(p/'out/eigenvectors.petsc'); logsha=sha(p/'solver.log')
    points=r['points']; actual=c.a.driver.table(p/'out/neutral_field.plt')
    assert points==actual
    valid=[x for x in points if c.a.valid(x)]
    assert len(valid)==len(points)>0
    runrows=[]
    for i,x in enumerate(points):
        matches=[j for j,z in enumerate(pairs) if abs(c.point_eigen(z)-c.point_eigen(x))<1e-8 and z['residual']<1e-7]
        assert matches,(r['name'],i)
        j=min(matches,key=lambda k:abs(c.point_eigen(pairs[k])-c.point_eigen(x)))
        q=qs[j]; assert np.isfinite(q).all() and np.vdot(q,q).real>0
        assert all(c.mac(q,qs[k])>=1-1e-8 for k in matches)
        z=pairs[j]
        row=dict(tag=d['tag'],kind=d['kind'],side=d.get('side') or '',run=r['name'],field_index=i,**x,
            polynomial_residual=z['residual'],alpha_vec_delta=abs(c.point_eigen(z)-c.point_eigen(x)),
            raw_vector_index=j,raw_vector_norm=float(np.linalg.norm(q)),raw_vector_count=len(qs),
            raw_vector_component_sha256=hashlib.sha256(q.astype('>c16').tobytes()).hexdigest(),
            raw_vector_file_sha256=rawsha)
        rows.append(row);runrows.append(row)
    end=min(runrows,key=lambda x:abs(x['beta_r']-d['endpoint'])) if d['kind']=='arm' else runrows[-1]
    peak=max(runrows,key=lambda x:x['growth'])
    if d['kind']=='arm': assert abs(end['beta_r']-d['endpoint'])<1e-8;witness.append(end)
    else:witness.append(end)
    summary.append(dict(tag=d['tag'],kind=d['kind'],side=d.get('side') or '',run=r['name'],
        omega_r=d['target_omega'],complete=True,source_identity_verified=True,samples=len(runrows),
        positive_samples=sum(x['growth']>0 for x in runrows),endpoint_beta=end['beta_r'],
        endpoint_alpha_r=end['alpha_r'],endpoint_alpha_i=end['alpha_i'],endpoint_growth=end['growth'],
        sampled_max_growth=peak['growth'],sampled_max_beta=peak['beta_r'],
        minimum_accepted_MAC=r['metrics']['minimum_accepted_mac'],
        maximum_accepted_residual=r['metrics']['maximum_accepted_residual'],
        maximum_output_pair_residual=max(x['polynomial_residual'] for x in runrows),
        vector_records=len(qs),all_output_points_vector_bound=True))
    bindings.append(dict(run=r['name'],result_file=str(file),result_sha256=sha(file),
        field_sha256=sha(p/'out/neutral_field.plt'),config_sha256=sha(p/'config.ini'),
        solver_log_sha256=logsha,vector_file=str(p/'out/eigenvectors.petsc'),
        vector_sha256=rawsha,vector_records=len(qs),solve_group_sizes=groups,
        field_points=len(runrows),all_field_points_matched=True))
    if d['kind']=='arm':curves[(d['tag'],d['side'])]=runrows
csvsave('native10_cut_samples.csv',rows);csvsave('native10_cut_summary.csv',summary);csvsave('native10_cut_witnesses.csv',witness)
fig,axes=plt.subplots(1,2,figsize=(11,4.2),sharey=True)
for ax,tag in zip(axes,['wplus005','w40']):
    for side in ['left','right']:
        rr=sorted(curves[tag,side],key=lambda x:x['beta_r'])
        ax.plot([x['beta_r'] for x in rr],[x['growth'] for x in rr],'.-',ms=3,lw=1,label=side+' path')
    w=next(x['omega_r'] for x in summary if x['tag']==tag)
    ax.axhline(0,c='black',lw=.8);ax.axvline(12.137900983290656,c='.65',ls=':',lw=.8)
    ax.set(xlim=(-41,41),xlabel=r'$\beta$',title=r'$\Omega=%.8g$'%w)
    ax.grid(alpha=.2);ax.legend(fontsize=8)
axes[0].set_ylabel(r'Spatial growth $g=-\mathrm{Im}(\alpha)$')
fig.suptitle('Native station 10: finite beta cuts; edge-positive, window unclosed',fontsize=11)
fig.tight_layout();fig.savefig(AREA/'native10_beta_cuts.png',dpi=190);plt.close(fig)
integrity=json.loads((AREA/'integrity.json').read_text())
now=[dict(path=x['path'],sha256=sha(x['path'])) for x in integrity['before']]
assert now==integrity['before']
solver=RUN/'outputs/bin/START_LST_wall_local_v10';frozen=sha(solver)
assert frozen=='634b45c5d30e076f5b86ccb136e0734c1905beedde108477dceb8a75adc3853c'
report=dict(reduced_at=stamp(),calculation_finished_at=progress['updated'],finished=True,complete_calls=6,
    total_output_samples=len(rows),vector_bound_output_samples=len(rows),positive_output_samples=sum(x['growth']>0 for x in rows),
    source_classification='beta_minimum_not_upper_maximum',source_fresh_join=json.loads((AREA/'source_fresh_join.json').read_text()),
    beta_domain=[-40,40],omega_cuts=[8.952811251576883,40.0],summary=summary,raw_bindings=bindings,
    positive_beta_endpoints=[x for x in witness if x['kind']=='arm'],unknown_calls=[],
    classification='finite_window_unclosed_edge_positive',
    qualification='Finite spatial PEP paths of this actual seed only. Positive g at beta boundaries and Omega=40 shows these sampling windows are not closed; it establishes neither physical absolute instability, downstream causality, nor exhaustive spectrum coverage.',
    original_inputs_and_original_binary_unchanged=now==integrity['before'],solver_sha256=frozen)
save('native10_cut_report.json',report)
save('reduction_integrity.json',dict(created=stamp(),current=now,unchanged=True,frozen_solver_sha256=frozen))
readme=['# native10 已完成切片归约','',
'六次调用均完整，所有 field 输出点已与 append Vec、日志本征值及残差逐点绑定。',
'源点是 β minimum，未当作上边界峰。两条有限频率切片为 Ω=8.952811251576883 和 40，β∈[-40,40]。',
'四个 β 边界均有实际正增长；Ω=40 也有正增长。这只说明当前有限窗口未闭合，不能判定物理 absolute 或下游因果传播，不能宣称全谱穷尽。','',
'计算完成于 2026-10-23 22:40:39+08:00；2026-10-24 依据保存原始结果归约，未重新求解。',
'详见 native10_cut_report.json、native10_cut_samples.csv、native10_cut_witnesses.csv 与 native10_beta_cuts.png。']
(AREA/'README.md').write_text('\n'.join(readme)+'\n')
owned=['reduce_cuts.py','native10_cut_report.json','native10_cut_samples.csv','native10_cut_summary.csv','native10_cut_witnesses.csv','native10_beta_cuts.png','reduction_integrity.json','README.md','cut_progress.json','cut_plan.json','source_fresh_join.json','integrity.json','run_cuts.py']
save('reduction_manifest.json',dict(created=stamp(),files=[dict(path=str(AREA/n),sha256=sha(AREA/n)) for n in owned],raw_bindings=bindings))
print(json.dumps(dict(complete_calls=6,samples=len(rows),all_points_vector_bound=True,classification=report['classification'],summary=summary),indent=2))

