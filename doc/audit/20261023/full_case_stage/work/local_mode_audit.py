"""Small immutable-v5 local-root frequency/MAC diagnostic; no mode sign filter."""
import concurrent.futures as cf
import argparse, configparser, csv, datetime, json, re, struct, subprocess, time
from pathlib import Path
import numpy as np
import h5py
import full_driver as driver

ROOT=driver.ROOT
BIN=ROOT/'outputs/bin/START_LST_wall_local_v5'
EXPECTED='528d6d6a1b7aab34c73a70b09721295089ba293b2ef2c9bec173aba8e589c266'
OUTPUT=ROOT/'outputs/local_mode_audit.json'
LOG=ROOT/'outputs/local_mode_audit_工作日志_2026-10-22.md'
FLOAT=r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?'
CENTERS=[
    ('local119_stable',119,0.,-7.,37.4096026990,.1089840811),
    ('local119_growing',119,0.,-7.,52.6648433525,-.0899795234),
    ('local1_largegrowth',1,0.,-7.,-4.0473993233,-5.7873998359),
    ('flat317_M1',317,8.17586,-.98903,30.87103096,-.8416318),
    ('flat312_M2seed',312,15.,15.,6.3395684955,-.5184689106),
]

def profile(station):
    with h5py.File(driver.CASE/'M1/in/grid_vertical2d_N.h5','r') as f:g=f['grid'][:,station-1,:,0]
    with h5py.File(driver.CASE/'M1/in/flow_vertical2d_N.h5','r') as f:q=f['baseflow'][:,station-1,:,0]
    normal=g[0]/np.linalg.norm(g[0]) if station<=241 else np.array([0.,1.])
    n=(g-g[0])@normal
    w99=np.flatnonzero(q[:,3]>=.99*q[-1,3]);j99=int(w99[0]) if len(w99) else None
    gradients=np.column_stack([np.gradient(q[:,k],n,edge_order=2) for k in range(6)])
    exterior_start=max(j99 if j99 is not None else 0,300)
    k=exterior_start+int(np.argmax(abs(gradients[exterior_start:,0])))
    return n,q,dict(normal=normal.tolist(),nmax_N=float(n[-1]),nmax_mm=float(35*n[-1]),
                    spanwise_W99_j=None if j99 is None else j99+1,
                    spanwise_W99_n_N=None if j99 is None else float(n[j99]),
                    spanwise_W99_qualification='First |W|/outer W crossing at 0.99 is a descriptive spanwise-velocity proxy, not a uniquely defined three-dimensional boundary-layer thickness.',
                    outer_density_gradient_max_j=k+1,outer_density_gradient_max_n_N=float(n[k]),
                    outer_gradient_search_begin_j=exterior_start+1,
                    outer_density_gradient_peak=float(gradients[k,0]))

def structure(vec,station):
    qv=vec.reshape((400,5));a=np.sum(abs(qv)**2,axis=1);total=float(a.sum())
    n,q,geo=profile(station);k=int(np.argmax(a));cum=np.cumsum(a)/total
    result=dict(geo,metric='unweighted squared primitive-vector amplitude; not physical energy',
                component_fractions=dict(zip(['rho','u_t','v_n','w','T'],(np.sum(abs(qv)**2,axis=0)/total).tolist())),
                peak_j=k+1,peak_n_N=float(n[k]),peak_n_mm=float(35*n[k]),peak_to_outer_boundary_N=float(n[-1]-n[k]),
                peak_baseflow=dict(zip(['rho','u_global','v_global','w','T','p'],q[k].tolist())),
                j_containing_50_percent=int(np.searchsorted(cum,.5))+1,
                j_containing_90_percent=int(np.searchsorted(cum,.9))+1,
                outer_10_percent_j_amplitude_fraction=float(a[360:].sum()/total),
                outer_10_percent_normal_distance_amplitude_fraction=float(a[n>=.9*n[-1]].sum()/total),
                interpretation='Index localization is not physical far-field localization on a stretched grid. Requires wall-normal grid/domain and outer-flow/shock sensitivity checks; no mode-identity assertion.')
    if geo['spanwise_W99_n_N'] is not None:
        result['peak_n_over_spanwise_W99']=float(n[k]/geo['spanwise_W99_n_N'])
        result['amplitude_fraction_below_spanwise_W99']=float(a[n<=geo['spanwise_W99_n_N']].sum()/total)
    return result

def log(message):
    stamp=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
    with LOG.open('a',encoding='utf-8') as f:f.write('\n'+stamp+' '+message+'\n')

def vectors(path):
    data=path.read_bytes();offset=0;arrays=[]
    while offset<len(data):
        cls,n=struct.unpack_from('>ii',data,offset);offset+=8
        assert cls==1211214 and n==2000,(cls,n,offset)
        end=offset+16*n;assert end<=len(data)
        arrays.append(np.frombuffer(data[offset:end],dtype='>c16').astype(np.complex128));offset=end
    return arrays

def raw_vectors(d,center,station):
    # PEP view files are overwritten after each solve; the final restore solve
    # is at the initial center. Parse only the last printed PEP eigenpair block.
    path=d/'out/all_vectors.petsc'
    if not path.exists():return dict(available=False,reason='no native vector output')
    vv=vectors(path)
    blocks=[];current=[]
    for line in (d/'solver.log').read_text(errors='replace').splitlines():
        m=re.search(r'^\s*(\d+)\s+\(\s*('+FLOAT+r')\s*,\s*('+FLOAT+r')\s*\)\s+The error is\s+('+FLOAT+r')',line)
        if not m:continue
        index=int(m[1])
        if index==0 and current:blocks.append(current);current=[]
        current.append(dict(index=index,alpha_r=float(m[2].replace('D','E')),alpha_i=float(m[3].replace('D','E')),residual=float(m[4].replace('D','E'))))
    if current:blocks.append(current)
    pairs=blocks[-1] if blocks else []
    if len(pairs)!=len(vv):return dict(available=False,reason='last raw pair/vector count mismatch',pair_count=len(pairs),vector_count=len(vv))
    a=complex(center['alpha_r'],center['alpha_i'])
    k=min(range(len(pairs)),key=lambda j:abs(complex(pairs[j]['alpha_r'],pairs[j]['alpha_i'])-a))
    error=abs(complex(pairs[k]['alpha_r'],pairs[k]['alpha_i'])-a)
    if error>1e-7:return dict(available=False,reason='final restore raw alpha does not match stencil center',nearest=pairs[k],alpha_difference=error)
    # Vector entries follow the five primitive variables at each j in the
    # verified solver ordering. These are unweighted amplitude measures.
    return dict(available=True,native_vector_file=str(path),native_values_file=str(d/'out/all_values.txt'),
                raw_index=k,raw_pair=pairs[k],center_alpha_difference=error,
                **structure(vv[k],station))

def probe(center,pep_tol='1e-12',nev=8,ncv=24,seed_shift=1e-6):
    label,station,omega,beta,ar,ai=center
    name='admissibility_'+label+'_v5'
    spec=dict(name=name,station=station,mode='M1',coordinates='wall_local',
              seed=dict(omega_r=omega,beta_r=beta,alpha_r=ar,alpha_i=ai),
              solver=str(BIN),solver_sha256=EXPECTED,task=33,sample_only=True,
              omega_min=omega-.01,omega_max=omega+.01,omega_points=5,
              beta_min=beta-1e-8,beta_max=beta+1e-8,beta_points=2,
              step=.005,seed_shift=seed_shift,eigen_tol=1e-7,jump=2,timeout=300,
              adaptive=True,minimum_step=1e-6,vector_tracking=True,minimum_overlap=.95,
              options=['-pep_tol',pep_tol,'-pep_max_it','1000','-pep_nev',str(nev),'-pep_ncv',str(ncv),
                       '-st_pc_factor_mat_solver_type','superlu',
                       '-pep_view_values','ascii:out/all_values.txt',
                       '-pep_view_vectors','binary:out/all_vectors.petsc','-options_left'])
    summary=driver.run(spec);d=ROOT/'runs'/name
    fields=driver.table(d/'out/neutral_field.plt')
    overlaps=[float(v.replace('D','E')) for v in re.findall(r'vector_overlap=\s*('+FLOAT+r')',(d/'out/neutral_run.log').read_text())]
    grouped=[]
    for om in [omega+delta for delta in (-.01,-.005,0.,.005,.01)]:
        rows=[r for r in fields if abs(r['omega_r']-om)<1e-10]
        accepted=[r for r in rows if r['status']==0]
        record=dict(requested_omega=om,rows=rows,accepted=len(accepted)==2)
        if len(accepted)==2:
            a,b=accepted
            record.update(alpha_r=(a['alpha_r']+b['alpha_r'])/2,alpha_i=(a['alpha_i']+b['alpha_i'])/2,
                          tiny_beta_alpha_difference=abs(complex(a['alpha_r'],a['alpha_i'])-complex(b['alpha_r'],b['alpha_i'])))
        grouped.append(record)
    result=dict(label=label,station=station,omega_center=omega,beta=beta,coordinates='wall_local',
                initial_guess=dict(alpha_r=ar,alpha_i=ai),summary=summary,stencil=grouped,
                actual_vector_overlap_count=len(overlaps),actual_minimum_vector_overlap=min(overlaps) if overlaps else None,
                physical_mode_identity='unclassified new local-model root' if station<242 else 'reference family seed; diagnostic only')
    if all(r['accepted'] for r in grouped):
        aa=[complex(r['alpha_r'],r['alpha_i']) for r in grouped]
        d01=(aa[4]-aa[0])/.02;d005=(aa[3]-aa[1])/.01
        d4=(aa[0]-8*aa[1]+8*aa[3]-aa[4])/.06
        def z(v):return dict(real=v.real,imag=v.imag)
        result['dispersion']=dict(dalpha_domega_h001=z(d01),dalpha_domega_h0005=z(d005),dalpha_domega_fourth_order=z(d4),
                                   derivative_relative_difference=abs(d01-d005)/max(abs(d005),1e-30),
                                   reciprocal_derivative=z(1/d4),
                                   diagnostic_only='Reciprocal complex derivative along this fixed real-beta branch; its real part is not a standalone causal/radiation admissibility filter.')
        result['native_center_structure']=raw_vectors(d,grouped[2],station)
    return result

def nearwall():
    assert driver.sha(BIN)==EXPECTED
    log('新阶段：主任务要求对近壁候选35.9340816817−4.0232078978i执行tol1e−12/1e−13短频率窗口，优先nev1/ncv8、sigma偏移1e−4，输出独立M3near初值；原M3外层根保留独立。此检查不证明物理身份或全程内存无泄漏。')
    r12=probe(('local119_nearwall_tol12',119,0.,-7.,35.9340816817,-4.0232078978),pep_tol='1e-12',nev=1,ncv=8,seed_shift=1e-4)
    a12=r12.get('native_center_structure',{}).get('raw_pair')
    assert a12 and all(x['accepted'] for x in r12['stencil']), 'First near-wall stencil failed'
    r13=probe(('local119_nearwall_tol13',119,0.,-7.,a12['alpha_r'],a12['alpha_i']),pep_tol='1e-13',nev=1,ncv=8,seed_shift=1e-4)
    a13=r13.get('native_center_structure',{}).get('raw_pair')
    good=bool(a13 and all(x['accepted'] for x in r13['stencil']))
    report=dict(solver=str(BIN),solver_sha256=EXPECTED,runs=[r12,r13],passed=good,
                physical_role='near-wall candidate, not final validation',
                caveats=['This is a short-call comparison, not a long-running memory-leak certification.',
                         'Near-wall shape, residual and local MAC continuity do not establish causal/radiation conditions, physical mode identity or normal-grid convergence.'])
    if good:
        v12=vectors(Path(r12['native_center_structure']['native_vector_file']))[r12['native_center_structure']['raw_index']]
        v13=vectors(Path(r13['native_center_structure']['native_vector_file']))[r13['native_center_structure']['raw_index']]
        mac=float(abs(np.vdot(v12,v13))**2/(np.vdot(v12,v12).real*np.vdot(v13,v13).real))
        differences=[]
        for x,y in zip(r12['stencil'],r13['stencil']):
            differences.append(abs(complex(x['alpha_r'],x['alpha_i'])-complex(y['alpha_r'],y['alpha_i'])))
        report['tightening_comparison']=dict(center_alpha_difference=abs(complex(a12['alpha_r'],a12['alpha_i'])-complex(a13['alpha_r'],a13['alpha_i'])),
                    max_stencil_alpha_difference=max(differences),center_native_vector_MAC=mac,
                    normal_structure_peak_j=[r['native_center_structure']['peak_j'] for r in (r12,r13)],
                    normal_structure_peak_n_N=[r['native_center_structure']['peak_n_N'] for r in (r12,r13)])
        passed=all(r['summary']['rejected']==0 and r['summary']['petsc_errors']==0 and r['summary']['unused_options']==0 and
                   r['actual_minimum_vector_overlap']>=.95 and r['native_center_structure']['raw_pair']['residual']<=1e-7 for r in (r12,r13))
        passed &= report['tightening_comparison']['max_stencil_alpha_difference']<1e-6 and mac>=.9999
        report['passed']=bool(passed)
        if passed:
            # The native final-restore pair is an exact fixed-beta, center-omega
            # solve, unlike the average of the two tiny beta perturbations.
            n,q,geo=profile(119)
            with h5py.File(driver.CASE/'M1/in/grid_vertical2d_N.h5','r') as f:x=float(f['grid'][0,118,0,0])
            seed=dict(family='M3near',physical_role='near-wall candidate, not final validation',
                      physical_mode_identity='unclassified new local-model near-wall root',
                      iloc=119,station=119,x=x,x_N=x,x_mm=35*x,Re=79100.,omega_r=0.,omega_i=0.,beta_r=-7.,beta_i=0.,
                      alpha_r=a13['alpha_r'],alpha_i=a13['alpha_i'],growth=-a13['alpha_i'],status=0,
                      coordinates='wall_local',mode='M1',source_mode='M1 main841',
                      solver=str(BIN),solver_sha256=EXPECTED,source_solver_sha256=EXPECTED,
                      source_run=r13['summary']['spec']['name'],source_path=str(ROOT/'runs'/r13['summary']['spec']['name']/'solver.log'),
                      source_kind='final restore raw PEP pair at exactly omega=0 beta=-7',raw_PEP_index=a13['index'],
                      absolute_polynomial_residual=a13['residual'],
                      native_vector_file=r13['native_center_structure']['native_vector_file'],
                      raw_vector_index=r13['native_center_structure']['raw_index'],
                      probe_report=str(ROOT/'outputs/local_mode_audit_nearwall.json'),
                      checked_frequency_offsets=[-.01,-.005,0.,.005,.01],minimum_MAC=min(r['actual_minimum_vector_overlap'] for r in (r12,r13)),
                      normal_structure=r13['native_center_structure'])
            path=ROOT/'outputs/local_mode_audit_nearwall_seed.json';path.write_text(json.dumps(seed,indent=2,allow_nan=False)+'\n')
            report['initial_seed_record']=str(path)
    report['solver_unchanged']=driver.sha(BIN)==EXPECTED
    (ROOT/'outputs/local_mode_audit_nearwall.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    log('近壁候选收紧诊断结束：'+json.dumps({k:report.get(k) for k in ('passed','tightening_comparison','initial_seed_record','solver_unchanged')},ensure_ascii=False))
    print(json.dumps({k:report.get(k) for k in ('passed','tightening_comparison','initial_seed_record','solver_unchanged')},indent=2),flush=True)

def target_probe(target):
    name=f'admissibility_local119_target{target:03d}_v5';d=ROOT/'runs'/name;d.mkdir(exist_ok=True)
    out=d/'diagnostic.json'
    if out.exists():return json.loads(out.read_text())
    assert not (d/'solver.log').exists(),'incomplete existing target diagnostic'
    cfg=configparser.ConfigParser();cfg.read(driver.CASE/'M1/config.ini')
    cfg['filenames']['dir']=str(driver.CASE/'M1/in')+'/'
    cfg['domain']['iloc']='119';cfg['lst']['profile_coordinates']='wall_local'
    cfg['lst']['alpha']=f'{target}, -0.3';cfg['lst']['omega']='0, 0';cfg['lst']['beta']='-7, 0'
    cfg['neutral']={'enabled':'F'}
    with (d/'config.ini').open('w') as f:cfg.write(f)
    options=['-pep_tol','1e-12','-pep_max_it','1000','-st_pc_factor_mat_solver_type','superlu','-pep_nev','8','-pep_ncv','24',
             '-pep_view_values','ascii:out/all_values.txt','-pep_view_vectors','binary:out/all_vectors.petsc','-options_left']
    cmd=[str(BIN),'-f','config.ini']+options;t=time.time()
    with (d/'solver.log').open('w') as f:
        f.write('COMMAND '+json.dumps(cmd)+'\n');f.flush()
        try:p=subprocess.run(cmd,cwd=d,env=driver.ENV,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,timeout=180);code=p.returncode
        except subprocess.TimeoutExpired:code=124
    text=(d/'solver.log').read_text(errors='replace');pairs=[]
    for line in text.splitlines():
        m=re.search(r'^\s*(\d+)\s+\(\s*('+FLOAT+r')\s*,\s*('+FLOAT+r')\s*\)\s+The error is\s+('+FLOAT+r')',line)
        if m:pairs.append(dict(index=int(m[1]),alpha_r=float(m[2].replace('D','E')),alpha_i=float(m[3].replace('D','E')),absolute_polynomial_residual=float(m[4].replace('D','E'))))
    vv=vectors(d/'out/all_vectors.petsc') if (d/'out/all_vectors.petsc').exists() else []
    matching=bool(len(vv)==len(pairs) and len(vv)>0 and [p['index'] for p in pairs]==list(range(len(pairs))))
    for k,pair in enumerate(pairs):
        pair['growth']=-pair['alpha_i'];pair['residual_accepted']=pair['absolute_polynomial_residual']<=1e-7
        pair['physical_mode_identity']='unclassified local-model eigenpair'
        if matching:pair['structure']=structure(vv[k],119)
    result=dict(name=name,target_alpha_r=target,target_alpha_i=-.3,station=119,omega=0,beta=-7,
                solver=str(BIN),solver_sha256=driver.sha(BIN),seconds=round(time.time()-t,3),exit=code,
                petsc_errors=text.count('PETSC ERROR'),unused_options=text.count('Option left:'),
                raw_vector_index_association_verified=matching,command=cmd,pairs=pairs)
    out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('TARGET '+json.dumps({k:result[k] for k in ('name','seconds','exit','petsc_errors','unused_options','raw_vector_index_association_verified')}),flush=True)
    return result

def targets():
    assert driver.sha(BIN)==EXPECTED
    log('阶段审查：前五频率窗口均0拒绝；i119两根按索引在j383/385集中。原网格强拉伸，不能据末40索引点直接认定贴外边界或伪谱，补充实际normal距离、W99描述尺度、外层rho梯度峰及8个alpha shift候选搜索。')
    report=dict(solver=str(BIN),solver_sha256=EXPECTED,scope='Eight finite shifts at i119 omega=0 beta=-7; not a full-spectrum search.',targets=[])
    out=ROOT/'outputs/local_mode_audit_targets.json'
    with cf.ThreadPoolExecutor(max_workers=3) as ex:
        jobs={ex.submit(target_probe,target):target for target in (10,20,30,40,50,60,80,100)}
        for job in cf.as_completed(jobs):
            try:r=job.result()
            except Exception as e:r=dict(target_alpha_r=jobs[job],failed=True,error=repr(e))
            report['targets'].append(r);out.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    allpairs=[]
    for result in report['targets']:
        for pair in result.get('pairs',[]):
            row=dict(target=result['target_alpha_r'],source_run=result['name'],**pair)
            alpha=complex(row['alpha_r'],row['alpha_i'])
            close=next((r for r in allpairs if abs(alpha-complex(r['alpha_r'],r['alpha_i']))<1e-5),None)
            if close:close['seen_targets'].append(result['target_alpha_r'])
            else:allpairs.append(dict(row,seen_targets=[result['target_alpha_r']]))
    report['deduplication_absolute_alpha_tolerance']=1e-5;report['unique_pairs']=sorted(allpairs,key=lambda r:r['alpha_r'])
    report['solver_unchanged']=driver.sha(BIN)==EXPECTED
    out.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    log('八个有限alpha shift诊断完成；详细结果 '+str(out)+'；不以growth或alpha_r符号删除候选。')

def report_files():
    report=json.loads(OUTPUT.read_text());targets=json.loads((ROOT/'outputs/local_mode_audit_targets.json').read_text())
    assert report['solver_unchanged'] and targets['solver_unchanged'] and driver.sha(BIN)==EXPECTED
    compact=[];profiles=[]
    for r in report['results']:
        s=r['native_center_structure'];d=r['dispersion'];a=s['raw_pair']
        compact.append(dict(label=r['label'],i=r['station'],omega=r['omega_center'],beta=r['beta'],
                            alpha_r=a['alpha_r'],alpha_i=a['alpha_i'],growth=-a['alpha_i'],residual=a['residual'],
                            min_MAC=r['actual_minimum_vector_overlap'],derivative_relative_change=d['derivative_relative_difference'],
                            reciprocal_derivative_real=d['reciprocal_derivative']['real'],reciprocal_derivative_imag=d['reciprocal_derivative']['imag'],
                            peak_j=s['peak_j'],peak_n_N=s['peak_n_N'],nmax_N=s['nmax_N'],W99_n_N=s['spanwise_W99_n_N'],
                            peak_n_over_W99=s['peak_n_over_spanwise_W99'],outer_density_gradient_peak_n_N=s['outer_density_gradient_max_n_N'],
                            W99_inner_amplitude_fraction=s['amplitude_fraction_below_spanwise_W99']))
        v=vectors(Path(s['native_vector_file']))[s['raw_index']];profiles.append((r['label'],r['station'],v))
    with (ROOT/'outputs/local_mode_audit_centers.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(compact[0]));w.writeheader();w.writerows(compact)
    candidates=[]
    for r in targets['unique_pairs']:
        s=r['structure'];candidates.append(dict(alpha_r=r['alpha_r'],alpha_i=r['alpha_i'],growth=r['growth'],residual=r['absolute_polynomial_residual'],
                    residual_accepted=r['residual_accepted'],seen_targets=';'.join(map(str,r['seen_targets'])),source_run=r['source_run'],raw_index=r['index'],
                    peak_j=s['peak_j'],peak_n_N=s['peak_n_N'],peak_n_over_W99=s['peak_n_over_spanwise_W99'],W99_inner_amplitude_fraction=s['amplitude_fraction_below_spanwise_W99']))
    with (ROOT/'outputs/local_mode_audit_candidates.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(candidates[0]));w.writeheader();w.writerows(candidates)
    new=min(targets['unique_pairs'],key=lambda r:abs(complex(r['alpha_r'],r['alpha_i'])-complex(35.9340816817,-4.0232078978)))
    v=vectors(ROOT/'runs'/new['source_run']/'out/all_vectors.petsc')[new['index']]
    profiles.append(('local119_new35.934_candidate',119,v))
    with (ROOT/'outputs/local_mode_audit_profiles.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['label','i','j','n_N','n_over_W99','total_primitive_amplitude_fraction','rho_fraction','ut_fraction','vn_fraction','w_fraction','T_fraction'])
        for label,station,v in profiles:
            n,q,geo=profile(station);a=abs(v.reshape((400,5)))**2;a/=a.sum()
            for j in range(400):w.writerow([label,station,j+1,n[j],n[j]/geo['spanwise_W99_n_N'],a[j].sum(),*a[j]])
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(12,4.8))
    for label,station,v in profiles:
        n,q,geo=profile(station);a=np.sum(abs(v.reshape((400,5)))**2,axis=1);a/=a.max();x=n/geo['spanwise_W99_n_N']
        for ax in axes:ax.plot(x,np.maximum(a,1e-50),label=label,lw=1.5)
    for ax in axes:
        ax.set_yscale('log');ax.set_ylim(1e-10,1.5);ax.set_xlabel('n / spanwise W99 distance');ax.grid(alpha=.25)
        ax.axvline(1,color='gray',ls='--',lw=1);ax.set_ylabel('Primitive-vector amplitude squared / peak')
    axes[0].set_xscale('log');axes[0].set_xlim(.03,180);axes[0].set_title('Entire normal range (log distance)')
    axes[1].set_xlim(0,5);axes[1].set_title('Near-wall comparison')
    handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=3,fontsize=8)
    fig.suptitle('Frozen v5 raw-vector diagnostics; identities unclassified; W99 is a descriptive proxy',fontsize=10)
    fig.tight_layout(rect=(0,.13,1,.95));fig.savefig(ROOT/'outputs/local_mode_audit_profiles.png',dpi=180);plt.close(fig)
    lines=['# 局部模型根的频率与法向结构小诊断','',
           '2026-10-22，Asia/Shanghai。使用冻结 v5 SHA-256 `'+EXPECTED+'`；原程序、原输入和源码均未修改。','',
           '五个中心各实际计算 omega−0.01、−0.005、0、+0.005、+0.01。每个频率的 beta 为目标±1e−8，用两行均值逼近固定 beta；所有 50 行 status=0，五个窗口均无拒绝，PETSc 错误和未使用参数为零。`sample_only` 扫描的 no bracketed root 是正常的采样终态，不是失败。每一步归一化欧氏向量 overlap≥0.9998602；这不等价于物理能量或因果性验证。','',
           '原生 PEP 二进制向量逐个与原始 PEP 索引绑定。五个窗口的最终恢复求解都与中心 alpha 对应，误差最大小于2e−13；向量按每个法向节点五个 primitive variables 的 block5 顺序解释。小残差只约束代数方程，不证明物理模态身份。','',
           '|中心|alpha|Re reciprocal derivative|Im reciprocal derivative|峰 n_N|W99尺度 n_N|',
           '|---|---|---:|---:|---:|---:|']
    for r in compact:lines.append('|%s|%.8f%+.8fi|%.8g|%.8g|%.8g|%.8g|'%(r['label'],r['alpha_r'],r['alpha_i'],r['reciprocal_derivative_real'],r['reciprocal_derivative_imag'],r['peak_n_N'],r['W99_n_N']))
    lines += ['',
       'i119 两个原 bootstrap 根的模方峰在 j383/j385、n=3.82665/4.11925，域端 nmax=7.69917；它们距外边界仍3.87252/3.57992。因此“末40个索引点占比接近1”不能说明贴外边界。增长根峰与外层密度梯度最大位置 j383 相同，提示它与原基流外层陡变相关；尚未判定是物理激波/外层模态、连续谱离散或网格效应。两者仍保留，须做法向网格/外域/模型适用性验证。',
       'i1 大增长根的 reciprocal derivative≈−0.0396985−0.0357413i，其峰在 n=.760219，而域端5.85714。负实部是色散诊断，不能单凭该量或 alpha_r 符号判定可接受/不可接受。',
       'W99 指原基流展向速度首次达到外边界 W 的99%的位置，只是描述尺度；曲面三维流动及外层激波下它不是唯一的边界层厚度。',
       '',
       '补充搜索 i119、omega=0、beta=−7：alpha 实部目标10、20、30、40、50、60、80、100，虚部初猜−.3，每个要求 nev8/ncv24。8次实际单点求解均无 PETSc 错误和未使用选项，返回66个 raw pairs，按 |delta alpha|<1e−5 得到58个候选。稳定根和不稳定根全部记录；有限shift搜索不是全谱。',
       '新候选 alpha≈35.9340816817−4.0232078978i 的绝对多项式残差≈1.06e−9，峰 n=.03903696，W99尺度内无权模方≈98.41%。它更靠近壁面，但其物理身份与沿程连续性、网格独立性仍待验证，不能仅凭结构或增长直接认定为目标模态。',
       '',
       '理论依据：Mack (1984), Boundary-Layer Linear Stability Theory, AGARD R709 Part3，§2.3.3、p20、式(2.31)–(2.35)，说明非保守稳定性波的群速通常为复数，Gaster关系仅在小增长近似下使用。本文导出的固定实 beta 分支 d alpha/d omega 及其倒数只作为局部复色散诊断；不把 Re倒数当单独辐射/因果筛选器。原资料：[Mack全文](https://transition.arizona.edu/sites/default/files/2021-08/mack-bllst.pdf)。',
       'Gaster (1962), JFM 14(2),222–224，DOI10.1017/S0022112062001184；本次查到出版商摘要，未获取全文。摘要明确限定小增长率。引用限于这一适用性：[Cambridge出版商页面](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/note-on-the-relation-between-temporallyincreasing-and-spatiallyincreasing-disturbances-in-hydrodynamic-stability/0777D0DF4F53DB5542CCEBEA551E1F17)。',
       '检索记录：DOI直接open及一次猜测的Cambridge URL均失败；后通过确切标题检索获得上述出版商页面。WSL rg不存在，已使用Python定点读取源码核对block5顺序。',
       '',
       '交付：local_mode_audit.json 五窗口详细结果，local_mode_audit_targets.json 全部shift/原始候选/结构与去重结果；centers.csv、candidates.csv、profiles.csv及profiles.png；可复现脚本 work/local_mode_audit.py；全部 runs/admissibility_* 原始输入/日志/原生向量保留。',
       '结束审查：本小诊断完成；没有删除任何用户文件，没有生成需删除的失败尝试代码。保留一个复现脚本和必要输出。未完成事项为根的物理身份、因果/辐射条件及法向网格/外域验证，已交回主任务。']
    (ROOT/'outputs/local_mode_audit_说明.md').write_text('\n\n'.join(lines)+'\n',encoding='utf-8')
    log('最终审查完成：五个窗口50个采样行有效，8次有限shift单点求解、58个去重根及native向量结构已保存。实际normal距离核对纠正了按索引误判远场的风险。物理身份/法向独立性未证明，明确保留待核验。生成CSV、PNG及说明文档；仅保留复现代码和必要原始证据，无需清理的本任务尝试代码。')
    print(json.dumps(dict(centers=len(compact),target_runs=len(targets['targets']),raw_pairs=sum(len(x['pairs']) for x in targets['targets']),unique_pairs=len(candidates)),indent=2))

def main():
    assert driver.sha(BIN)==EXPECTED
    log('SHA-256 核验通过；处理五个中心的实际频率窗口，最多三个并发，实 beta 固定到 ±1e-8 采样误差；已存在且spec完全一致的结果只读复用。')
    report=dict(started=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),
                solver=str(BIN),solver_sha256=EXPECTED,results=[],
                qualification='Finite frequency/MAC diagnostic only. Neither alpha_r sign nor Re reciprocal derivative automatically removes roots. No physical causal mode classification is claimed.')
    with cf.ThreadPoolExecutor(max_workers=3) as ex:
        jobs={ex.submit(probe,c):c for c in CENTERS}
        for job in cf.as_completed(jobs):
            try:r=job.result()
            except Exception as e:r=dict(label=jobs[job][0],failed=True,error=repr(e))
            report['results'].append(r)
            OUTPUT.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
            log('完成 '+r['label']+'；结果 '+json.dumps({k:r.get(k) for k in ('failed','error','dispersion','native_center_structure')},ensure_ascii=False))
    report['solver_unchanged']=driver.sha(BIN)==EXPECTED
    report['completed']=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
    OUTPUT.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps([dict(label=r['label'],failed=r.get('failed'),dispersion=r.get('dispersion'),native_center_structure=r.get('native_center_structure')) for r in report['results']],indent=2),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--targets',action='store_true');parser.add_argument('--report',action='store_true');parser.add_argument('--nearwall',action='store_true');args=parser.parse_args()
    if args.nearwall:nearwall()
    elif args.report:report_files()
    elif args.targets:targets()
    else:main()
