"""Reduce saved finite i841 gap paths; no solver calls."""
import collections,csv,hashlib,json,pathlib,struct,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
AREA=pathlib.Path(__file__).resolve().parent;RUN=AREA.parent.parent
COMP=RUN/'outputs/upper_boundary_139_v10_20261023/compensation_20261023_2148'
sys.path.insert(0,str(COMP));import compensate as c
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def write(n,d): (AREA/n).write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def csvsave(n,rows):
    fields=list(dict.fromkeys(k for x in rows for k in x)) or ['run','iloc','omega_r','beta_r','alpha_r','alpha_i','growth','vector_bound','polynomial_residual']
    with (AREA/n).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
def parse_run(d):
    summary=json.loads((d/'summary.json').read_text());points=c.a.driver.table(d/'out/neutral_field.plt')
    log=(d/'solver.log').read_text(errors='replace');pairs=[];groups=[]
    for line in log.splitlines():
        if 'ncov=' in line:groups.append([])
        m=c.PAIR.search(line)
        if m:
            z=dict(index=int(m[1]),alpha_r=float(m[2].replace('D','E')),alpha_i=float(m[3].replace('D','E')),residual=float(m[4].replace('D','E')))
            assert groups;groups[-1].append(z);pairs.append(z)
    path=d/'out/eigenvectors.petsc';qs=[];errors=[]
    try:
        data=path.read_bytes();off=0
        while off<len(data):
            cls,n=struct.unpack_from('>ii',data,off);off+=8
            assert cls==1211214 and 0<n<1000000
            end=off+16*n;assert end<=len(data)
            qs.append(np.frombuffer(data[off:end],dtype='>c16').astype(np.complex128));off=end
        assert len(pairs)==len(qs)>0,(len(pairs),len(qs))
    except Exception as e:errors.append(repr(e))
    rows=[]
    for k,p in enumerate(points):
        z=dict(run=d.name,field_index=k,**p,vector_bound=False)
        if not errors and c.a.valid(p):
            matches=[j for j,a in enumerate(pairs) if abs(c.point_eigen(a)-c.point_eigen(p))<1e-8 and a['residual']<1e-7]
            if matches:
                j=min(matches,key=lambda j:abs(c.point_eigen(pairs[j])-c.point_eigen(p)))
                q=qs[j]
                consistent=np.isfinite(q).all() and np.linalg.norm(q)>0 and all(c.mac(q,qs[l])>=1-1e-8 for l in matches)
                if consistent:z.update(vector_bound=True,raw_vector_index=j,polynomial_residual=pairs[j]['residual'],
                    alpha_vec_delta=abs(c.point_eigen(pairs[j])-c.point_eigen(p)),raw_vector_norm=float(np.linalg.norm(q)),
                    raw_vector_component_sha256=hashlib.sha256(q.astype('>c16').tobytes()).hexdigest(),raw_vector_file_sha256=sha(path))
        rows.append(z)
    return dict(name=d.name,spec=summary['spec'],summary=summary,field_count=len(points),
        vector_records=len(qs),printed_pairs=len(pairs),group_sizes=[len(x) for x in groups],
        vector_bound_count=sum(x['vector_bound'] for x in rows),binding_errors=errors,
        raw_files=[dict(path=str(p),sha256=sha(p)) for p in [d/'config.ini',d/'solver.log',d/'out/neutral_field.plt',path] if p.exists()]),rows
def main():
    progress=json.loads((AREA/'audit_progress.json').read_text());assert progress['finished'] and progress['completed']==6
    plan=json.loads((AREA/'audit_plan.json').read_text());summaries=[];rows=[];runs=[];byname={}
    for d in sorted((AREA/'runs').iterdir()):
        if not (d/'summary.json').exists():continue
        raw,rs=parse_run(d);runs.append(raw);rows.extend(rs);byname[d.name]=rs
    byrun={r['name']:r for r in runs}
    positive=[p for p in rows if p['vector_bound'] and p['status']==0 and abs(p['omega_r']-plan['target_omega'])<1e-10 and p['growth']>1e-8]
    for r in progress['results']:
        arms=r.get('arms',[]);midnames=[a['name'] for a in arms]
        if r.get('center'):midnames.append(r['center']['name'])
        midrows=[p for name in midnames for p in byname.get(name,[]) if abs(p['omega_r']-plan['target_omega'])<1e-10 and p['status']==0]
        bound=bool(midrows) and all(p['vector_bound'] for p in midrows)
        pos=[p for p in midrows if p['vector_bound'] and p['growth']>1e-8]
        final='positive_actual_sample' if pos else 'sampled_nonpositive_known_path' if r['state']=='sampled_nonpositive_known_path' and bound else 'unknown'
        peak=max(midrows,key=lambda p:p['growth']) if midrows else None
        origin=r['job']['point'];fresh=r.get('fresh',{})
        summaries.append(dict(label=r['label'],state=final,source_omega=origin['omega_r'],source_beta=origin['beta_r'],
            source_alpha_r=origin['alpha_r'],source_alpha_i=origin['alpha_i'],target_omega=plan['target_omega'],
            source_fresh_original_alpha_delta=fresh.get('original_alpha_delta'),
            source_vector_sha256=fresh.get('vector',{}).get('sha256'),
            source_vector_index=fresh.get('vector',{}).get('vector_index'),
            center_complete=r.get('center',{}).get('passed',False),completed_arms=sum(a['complete'] for a in arms),
            all_output_vectors_bound=bound,midfrequency_samples=len(midrows),
            maximum_sampled_growth=peak['growth'] if peak else None,maximum_sampled_beta=peak['beta_r'] if peak else None,
            seconds=r.get('seconds'),reason=r.get('reason')))
    # Accepted path gates come from actual solver logs, separately from endpoint Vec joins.
    for s,r in zip(summaries,progress['results']):
        metrics=[stage['run']['metrics'] for stage in r.get('arms',[])]
        if r.get('center'):metrics.append(r['center']['run']['metrics'])
        s['minimum_accepted_MAC']=min((x['minimum_accepted_mac'] for x in metrics if x['minimum_accepted_mac'] is not None),default=None)
        s['maximum_accepted_residual']=max((x['maximum_accepted_residual'] for x in metrics if x['maximum_accepted_residual'] is not None),default=None)
    csvsave('gap_samples.csv',rows);csvsave('gap_source_summary.csv',summaries);csvsave('gap_positive_witnesses.csv',positive)
    fig,axes=plt.subplots(1,2,figsize=(11,4.6))
    for r in progress['results']:
        rr=[p for a in r.get('arms',[]) for p in byname.get(a['name'],[]) if p['status']==0 and p['vector_bound']]
        # Retain every independent seed path even if scalar eigenvalues later coincide.
        rr=sorted(rr,key=lambda x:x['beta_r'])
        if rr:
            for ax in axes:ax.plot([x['beta_r'] for x in rr],[x['growth'] for x in rr],'.-',lw=.8,ms=2.7,label=r['label'])
    for ax in axes:
        ax.axhline(0,color='k',lw=.8);ax.set(xlabel=r'$\beta$',ylabel=r'$g=-\mathrm{Im}(\alpha)$');ax.grid(alpha=.2)
    axes[0].set(xlim=(-41,41),title='Finite domain');axes[0].legend(fontsize=7)
    axes[1].set(xlim=(-5,10),ylim=(-.15,.005),title='Enlarged view near zero growth')
    fig.suptitle('i841 known-seed samples at Omega=%.10g; one source unknown'%plan['target_omega'],fontsize=11)
    fig.tight_layout();fig.savefig(AREA/'gap_beta_cuts.png',dpi=190);plt.close(fig)
    checks=json.loads((AREA/'audit_integrity.json').read_text())
    current=[dict(path=x['path'],sha256=sha(x['path'])) for x in checks['before']]
    assert current==checks['before']
    build=json.loads((RUN/'outputs/v10_build_validation.json').read_text());S=RUN/'work/solver_source'
    assert all(sha(S/rel)==h for rel,h in build['source_hashes'].items())
    report=dict(created=c.a.stamp(),finished=True,station=841,x_mm=1100.,target_omega=plan['target_omega'],
        overall_state='inconclusive_gap_known_paths_nonpositive_one_unknown',
        beta_domain=plan['beta_domain'],known_seed_count=6,actual_beta_arm_calls=sum(len(x.get('arms',[])) for x in progress['results']),
        state_counts=dict(collections.Counter(s['state'] for s in summaries)),source_summary=summaries,
        total_actual_calls=len(runs),total_output_points=len(rows),vector_bound_points=sum(x['vector_bound'] for x in rows),
        robust_positive_midfrequency_points=positive,results=progress['results'],raw_runs=runs,plan=plan,
        original_inputs_and_original_binary_unchanged=True,all_63_frozen_build_sources_unchanged=True,solver_sha256=sha(c.a.SOLVER),
        qualification='Only beta in [-40,40] sampled from six already known actual seed paths. No positive sample on a completed path is a finite sampled observation, not a full-beta or full-spectrum stability proof. Unknown stays unknown. Source tags do not establish distinct physical modes; no semantic merge was made without a same-parameter raw-vector join.')
    write('gap_audit_report.json',report)
    (AREA/'README.md').write_text('# i841 中频孔洞有限路径检查\n\n'+
        '目标频率 Ω='+str(plan['target_omega'])+'，β∈[-40,40]，6 个已知实际来源、最多 12 条 β arms。\n'+
        '结果：5 个来源完成中心和10条β arms，未采到正增长；M2_upper_peak_a中心120s超时、0个中频field输出，保留unknown。\n'+
        '最接近零的实际样本为M1high来源 β=3.048257434757353，g=-0.01180631477770426。粗步长采样不能排除点间窄峰。\n'+
        '每来源 fresh 实际根/Vec→实ω中心输运→左右β路径，跨进程同点以原始Vec MAC及标量join核验。\n'+
        '完成/失败状态见 gap_audit_report.json；gap_samples.csv 逐 field 点与 raw Vec/残差绑定；原始 config/log/Vec 保留。\n'+
        '这是有限已知路径采样，不能称全β或全谱稳定，也不能仅按来源标签将同根视为不同物理模式。\n')
    names=['run_gap_audit.py','reduce_gap_audit.py','audit_plan.json','audit_progress.json','audit_integrity.json',
        'gap_audit_report.json','gap_source_summary.csv','gap_samples.csv','gap_positive_witnesses.csv','gap_beta_cuts.png','README.md']
    write('manifest.json',dict(created=c.a.stamp(),files=[dict(path=str(AREA/n),sha256=sha(AREA/n)) for n in names],raw_runs=runs))
    print(json.dumps(dict(state_counts=report['state_counts'],output_points=len(rows),bound_points=report['vector_bound_points'],
        positives=len(positive),summaries=summaries,report_sha256=sha(AREA/'gap_audit_report.json')),indent=2))
if __name__=='__main__':main()

