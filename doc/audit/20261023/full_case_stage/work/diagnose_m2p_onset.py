"""Finite local (omega,beta) maximum of the continued M2p search branch."""
import argparse,csv,datetime,hashlib,json,math,re,struct
from pathlib import Path
import numpy as np
import h5py
from scipy.interpolate import CubicSpline
import full_driver as driver
import fd_envelope as fd

ROOT=driver.ROOT
BIN=ROOT/'outputs/bin/START_LST_wall_local_v7'
HASH='7aa6016570f406587156f229107c0e5e4f76c8a85ba89dda267039ff1d688d79'
PROGRESS=ROOT/'outputs/neutral_x_cont_v7a_M2p_upstream_progress.json'
GAP=ROOT/'outputs/root_gap_diagnostic_M2p_upstream.json'
OUT=ROOT/'outputs/local_mode_audit_M2p_onset_diagnostic.json'
LOG=ROOT/'outputs/local_mode_audit_工作日志_2026-10-23_M2p起始区诊断.md'
BOUNDS=np.array([[24.,34.],[30.,40.]])
FLOAT=r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?'

def stamp():return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
def log(message):
    with LOG.open('a',encoding='utf-8') as f:f.write('\n'+stamp()+' '+message+'\n')

def metrics(point):
    center=point['center'];h=point['h'];values={}
    for r in point['samples']:
        i=round((r['omega_r']-center['omega_r'])/h);j=round((r['beta_r']-center['beta_r'])/h)
        values[i,j]=r['growth']
    g=lambda i,j:values[i,j]
    grad=np.array([(g(1,0)-g(-1,0))/(2*h),(g(0,1)-g(0,-1))/(2*h)])
    cross=(g(1,1)-g(1,-1)-g(-1,1)+g(-1,-1))/(4*h*h)
    H=np.array([[(g(1,0)-2*g(0,0)+g(-1,0))/h**2,cross],
                [cross,(g(0,1)-2*g(0,0)+g(0,-1))/h**2]])
    return dict(gradient=grad.tolist(),Hessian=H.tolist(),Hessian_eigenvalues=np.linalg.eigvalsh(H).tolist())

def maximize(station,seed):
    result=dict(station=station,status='unresolved',bounds=BOUNDS.tolist(),trace=[],seed=seed,
                solver=str(BIN),solver_sha256=HASH,qualification='One locally continued M2p branch; finite window only, no global spectrum or all-mode stability claim.')
    current=dict(seed);w=seed['omega_r'];b=seed['beta_r'];h=.02
    try:
        for it in range(12):
            p=fd.stencil(station,current,w,b,BIN,'wall_local',h=h,label=f'M2ponset_max_i{station}')
            p.update(metrics(p));result['trace'].append(p)
            grad=np.array(p['gradient']);H=np.array(p['Hessian']);current=p['center']
            if max(p['Hessian_eigenvalues'])>=0:raise RuntimeError('Hessian not negative definite; local Newton maximum is unverified')
            if np.max(abs(grad))<1e-7:
                verify=fd.stencil(station,current,w,b,BIN,'wall_local',h=h/2,tight=True,label=f'M2ponset_verify_i{station}')
                verify.update(metrics(verify));result['verification']=verify
                correction=np.linalg.solve(np.array(verify['Hessian']),-np.array(verify['gradient']))
                if max(verify['Hessian_eigenvalues'])<0 and max(abs(np.array(verify['gradient'])))<1e-6 and max(abs(correction))<1e-4:
                    result.update(status='local_maximum_candidate',point=verify['center'],gradient=verify['gradient'],Hessian=verify['Hessian'],
                        Hessian_eigenvalues=verify['Hessian_eigenvalues'],verification_correction=correction.tolist(),
                        growth_difference_on_halving=verify['center']['growth']-p['center']['growth'])
                    break
                h/=2;p=verify;current=verify['center'];grad=np.array(verify['gradient']);H=np.array(verify['Hessian'])
            step=np.linalg.solve(H,-grad);step*=min(1.,1/max(abs(step)))
            proposed=np.array([w,b])+step
            if np.any(proposed<=BOUNDS[:,0]) or np.any(proposed>=BOUNDS[:,1]):
                result.update(status='boundary_truncated',reason='Maximum Newton step touches/leaves finite window',proposed=proposed.tolist());break
            w,b=map(float,proposed)
        else:result['reason']='maximum iterations reached; not evidence of no root'
    except Exception as e:result['reason']=repr(e)
    return result

def survey(station,seed,wl,wh,nw,bl,bh,nb,label,reverse_w=False,reverse_b=False,export=False):
    spec=dict(station=station,mode='M1',coordinates='wall_local',seed={k:seed[k] for k in ('omega_r','beta_r','alpha_r','alpha_i')},
        seed_shift=1e-4,solver=str(BIN),solver_sha256=HASH,task=33,sample_only=True,
        omega_min=wl,omega_max=wh,omega_points=nw,beta_min=bl,beta_max=bh,beta_points=nb,
        reverse_outer=reverse_w,reverse_inner=reverse_b,anchor_row_seed=True,adaptive=True,vector_tracking=True,
        minimum_overlap=.95,minimum_step=1e-6,retry_eigenpairs=8,step=.4,jump=2,
        options=['-pep_tol','1e-12','-pep_max_it','1000','-pep_nev','1','-pep_ncv','8','-st_pc_factor_mat_solver_type','superlu'],timeout=180)
    if export:spec['options']+=['-pep_view_vectors','binary:out/all_vectors.petsc','-pep_view_values','ascii:out/all_values.txt','-options_left']
    digest=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()[:12]
    spec['name']=f'M2ponset_{label}_i{station}_{digest}'
    s=driver.run(spec);path=ROOT/'runs'/spec['name']/'out/neutral_field.plt';rows=driver.table(path)
    valid=[r for r in rows if r['status']==0]
    return dict(name=spec['name'],summary=s,source_path=str(path),samples=rows,
        complete=bool(len(rows)==nw*nb and len(valid)==len(rows) and s['petsc_errors']==0 and s['unused_options']==0 and not s['coordinate_metadata_missing']),
        winner=max(valid,key=lambda r:r['growth']) if valid else None)

def save(state):
    state['updated']=stamp();temporary=OUT.with_suffix('.tmp')
    temporary.write_text(json.dumps(state,indent=2,allow_nan=False)+'\n');temporary.replace(OUT)

def native(survey_result,point):
    d=ROOT/'runs'/survey_result['name'];text=(d/'solver.log').read_text(errors='replace');blocks=[];block=[]
    for line in text.splitlines():
        m=re.search(r'^\s*(\d+)\s+\(\s*('+FLOAT+r')\s*,\s*('+FLOAT+r')\s*\)\s+The error is\s+('+FLOAT+r')',line)
        if not m:continue
        index=int(m[1])
        if index==0 and block:blocks.append(block);block=[]
        block.append(dict(index=index,alpha_r=float(m[2].replace('D','E')),alpha_i=float(m[3].replace('D','E')),residual=float(m[4].replace('D','E'))))
    if block:blocks.append(block)
    pairs=blocks[-1];data=(d/'out/all_vectors.petsc').read_bytes();off=0;vv=[]
    while off<len(data):
        cls,n=struct.unpack_from('>ii',data,off);off+=8;assert cls==1211214 and n==2000
        end=off+16*n;vv.append(np.frombuffer(data[off:end],dtype='>c16').astype(np.complex128).reshape(400,5));off=end
    assert len(pairs)==len(vv)
    z=complex(point['alpha_r'],point['alpha_i']);k=min(range(len(pairs)),key=lambda k:abs(complex(pairs[k]['alpha_r'],pairs[k]['alpha_i'])-z))
    error=abs(complex(pairs[k]['alpha_r'],pairs[k]['alpha_i'])-z);assert error<1e-7
    station=int(point['iloc']);n,q=normal_profile(station);a=np.sum(abs(vv[k])**2,axis=1);peak=int(np.argmax(a))
    j99=int(np.flatnonzero(q[:,3]>=.99*q[-1,3])[0])
    runlog=(d/'out/neutral_run.log').read_text();overlaps=[float(x.replace('D','E')) for x in re.findall(r'vector_overlap=\s*('+FLOAT+r')',runlog)]
    result=dict(station=station,source_run=survey_result['name'],raw_index=k,raw_pair=pairs[k],point_alpha_match_error=error,
                native_vector_file=str(d/'out/all_vectors.petsc'),native_vector_sha256=driver.sha(d/'out/all_vectors.petsc'),
                peak_j=peak+1,peak_n_N=float(n[peak]),nmax_N=float(n[-1]),W99_n_N=float(n[j99]),
                primitive_amplitude_fraction_below_W99=float(a[n<=n[j99]].sum()/a.sum()),
                minimum_step_MAC=min(overlaps) if overlaps else None,
                metric='Euclidean solver-vector MAC per continuation step; primitive amplitude localization is not energy or physical mode classification.')
    return result,vv[k],n

def normal_profile(station):
    case=ROOT.parent/'M1/in'
    with h5py.File(case/'grid_vertical2d_N.h5','r') as f:g=f['grid'][:,station-1,:,0]
    with h5py.File(case/'flow_vertical2d_N.h5','r') as f:q=f['baseflow'][:,station-1,:,0]
    normal=g[0]/np.linalg.norm(g[0]);n=(g-g[0])@normal
    return n,q

def shapes(state,stations):
    state.setdefault('native_peaks',[]);cache={}
    for maximum in state['maxima']:
        station=maximum['station']
        if station not in stations or maximum['status']!='local_maximum_candidate':continue
        p=maximum['point'];b=p['beta_r'];w=p['omega_r']
        s=survey(station,p,w,w,1,b-1e-8,b+1e-8,2,'native_peak',export=True)
        record,v,n=native(s,p);cache[station]=(record,v,n)
        if not any(x['station']==station for x in state['native_peaks']):state['native_peaks'].append(record)
        save(state)
    state['peak_shape_comparisons']=[]
    if 205 in cache:
        _,a,n0=cache[205]
        for station,(_,v,n1) in cache.items():
            if station==205:continue
            keep=n0<=n1[-1];n=n0[keep];aa=a[keep];bb=CubicSpline(n1,v,axis=0,extrapolate=False)(n)
            weights=np.empty(len(n));weights[1:-1]=(n[2:]-n[:-2])/2;weights[0]=(n[1]-n[0])/2;weights[-1]=(n[-1]-n[-2])/2;weights=weights[:,None]
            mac=float(abs(np.sum(weights*np.conj(aa)*bb))**2/(np.sum(weights*abs(aa)**2)*np.sum(weights*abs(bb)**2)))
            state['peak_shape_comparisons'].append(dict(reference=205,station=station,normal_integral_equal_primitive_MAC=mac,
                qualification='Physical-normal CubicSpline map, no extrapolation; equal primitive weights and local tangent/normal components, not energy.'))
        save(state)
    state.setdefault('cross_station_checks',[])
    bystation={r['station']:r for r in state['maxima']}
    for station in [s for s in stations if s in [206,207]]:
        if any(x['station']==station for x in state['cross_station_checks']):continue
        p=bystation[205]['point'];target=bystation[station]['point'];w=p['omega_r'];b=p['beta_r']
        spec=dict(station=station,seed_station=205,mode='M1',coordinates='wall_local',seed={k:p[k] for k in ('omega_r','beta_r','alpha_r','alpha_i')},
             seed_shift=1e-4,solver=str(BIN),solver_sha256=HASH,task=33,sample_only=True,
             outer_parameter='iloc',inner_parameter='beta',outer_min=station,outer_max=station,outer_points=1,
             inner_min=b-1e-8,inner_max=b+1e-8,inner_points=2,anchor_row_seed=True,adaptive=True,
             vector_tracking=True,minimum_overlap=.95,minimum_step=1e-6,retry_eigenpairs=8,step=.4,jump=2,
             options=['-pep_tol','1e-12','-pep_max_it','1000','-pep_nev','1','-pep_ncv','8','-st_pc_factor_mat_solver_type','superlu'],timeout=180)
        digest=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()[:12];spec['name']=f'M2ponset_transport205to{station}_{digest}'
        summary=driver.run(spec);field=driver.table(ROOT/'runs'/spec['name']/'out/neutral_field.plt')
        valid=[v for v in field if v['status']==0 and int(v['iloc'])==station]
        check=dict(station=station,from_station=205,transport_run=spec['name'],transport_summary=summary,transport_field=field,status='unknown')
        if len(valid)==2 and summary['petsc_errors']==0:
            transported=valid[0]
            follow=survey(station,transported,target['omega_r'],target['omega_r'],1,target['beta_r']-1e-8,target['beta_r']+1e-8,2,'transport_to_peak')
            check['follow']=follow
            if follow['complete']:
                mid=follow['winner'];error=abs(complex(mid['alpha_r'],mid['alpha_i'])-complex(target['alpha_r'],target['alpha_i']))
                check.update(status='endpoint_matched' if error<1e-6 else 'endpoint_mismatch',endpoint_alpha_difference=error,
                             note='Actual iloc/MAC transfer at fixed parameters followed by same-station parameter continuation to previously found maximum; beta differs by 1e-8 in sampled endpoint.')
        state['cross_station_checks'].append(check);save(state)

def finalize(state):
    assert len(state['maxima'])==3 and all(r['status']=='local_maximum_candidate' for r in state['maxima'])
    assert len(state['surveys'])==15 and all(r['complete'] for r in state['surveys'])
    assert len(state.get('native_peaks',[]))==3 and all(c['status']=='endpoint_matched' for c in state.get('cross_station_checks',[]))
    inventory=[]
    for d in sorted((ROOT/'runs').iterdir()):
        if not d.name.startswith('M2ponset_'):continue
        path=d/'summary.json';assert path.exists(),str(d)
        s=json.loads(path.read_text());assert s['solver_sha256']==HASH
        text=(d/'out/neutral_run.log').read_text()
        accepted=[float(v.replace('D','E')) for v in re.findall(r'vector_overlap=\s*('+FLOAT+r')\s+accepted=\s*T',text)]
        inventory.append(dict(name=d.name,summary_path=str(path),summary_sha256=driver.sha(path),
            started=s['started'],seconds=s['seconds'],terminal=s['terminal'],exit=s['exit'],
            solves=s['solves'],evaluations=s['evaluations'],rejected=s['rejected'],field_records=s['field_records'],
            failed_field_records=s['failed_field_records'],petsc_errors=s['petsc_errors'],unused_options=s['unused_options'],
            accepted_step_MAC_min=min(accepted) if accepted else None,
            adaptive_rejected_step_overlap_records=len(re.findall(r'vector_overlap=.*accepted=\s*F',text))))
    assert len(inventory)==42
    assert all(r['failed_field_records']==0 and r['rejected']==0 and r['petsc_errors']==0 and r['unused_options']==0 for r in inventory)
    state['run_inventory']=inventory;state['finite_window_summaries']=[]
    compact=[]
    for m in sorted(state['maxima'],key=lambda r:r['station']):
        i=m['station'];p=m['point'];surveys=[s for s in state['surveys'] if s['station']==i and s['kind'].startswith('finitewindow')]
        allrows=[v for s in surveys for v in s['samples']];unique={(round(v['omega_r'],9),round(v['beta_r'],9)):v for v in allrows}
        edge=[v for v in unique.values() if min(abs(v['omega_r']-24),abs(v['omega_r']-34),abs(v['beta_r']-30),abs(v['beta_r']-40))<1e-8]
        winner=max(unique.values(),key=lambda v:v['growth']);edgewinner=max(edge,key=lambda v:v['growth'])
        small=next(s for s in state['surveys'] if s['station']==i and s['kind']=='smallwindow')
        near_edge=abs(small['winner']['omega_r']-p['omega_r'])>.24 or abs(small['winner']['beta_r']-p['beta_r'])>.24
        boundary_winner=min(abs(winner['omega_r']-24),abs(winner['omega_r']-34),abs(winner['beta_r']-30),abs(winner['beta_r']-40))<1e-8
        finite=dict(station=i,raw_records=len(allrows),unique_samples=len(unique),boundary_samples=len(edge),
            sampled_winner=winner,boundary_winner=edgewinner,winner_touches_window_boundary=boundary_winner,
            smallwindow_winner_touches_boundary=near_edge,
            interpretation='Finite sampled window plus local maximum evidence; no mathematical proof that all unsampled interior modes/points have lower growth.')
        state['finite_window_summaries'].append(finite)
        assert len(unique)==25 and len(edge)==16 and not boundary_winner and not near_edge
        H=np.array(m['Hessian']);effective=float(H[0,0]-H[0,1]**2/H[1,1])
        compact.append(dict(station=i,x_mm=35*p['x'],omega=p['omega_r'],beta=p['beta_r'],alpha_r=p['alpha_r'],alpha_i=p['alpha_i'],
            max_growth=p['growth'],gradient_inf_norm=float(max(abs(np.array(m['gradient'])))),
            Hessian_eigenvalue_1=m['Hessian_eigenvalues'][0],Hessian_eigenvalue_2=m['Hessian_eigenvalues'][1],
            finite_boundary_sampled_max_growth=edgewinner['growth'],beta_optimized_effective_omega_curvature=effective,
            local_quadratic_half_width_estimate=float(math.sqrt(-2*p['growth']/effective)) if p['growth']>0 else None))
    bystation={r['station']:r for r in compact}
    state['onset_interpretation']=dict(status='local_branch_onset_bracket_evidence',
        negative_peak_station=205,positive_peak_station=206,
        x_mm_bracket=[bystation[205]['x_mm'],bystation[206]['x_mm']],
        statement='The continued M2p local growth maximum changes from negative at i205 to positive at i206; consistent with onset/disappearance of this local instability island. This does not establish all-branch stability, prove complete absence of roots in the finite window, or solve an onset at fractional i/x.',
        exact_critical_x_solved=False)
    state['input_integrity_against_original_manifest']=driver.snapshot()
    state['solver_unchanged']=driver.sha(BIN)==HASH
    state['source_documents_current_hashes']={str(p):driver.sha(p) for p in (PROGRESS,GAP)}
    state['station_summary']=compact;state['completed']=stamp();save(state)
    with (ROOT/'outputs/local_mode_audit_M2p_onset_summary.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(compact[0]));w.writeheader();w.writerows(compact)
    lines=['# M2p 上游停止点：局部增长峰与沿程起始区诊断','',stamp(),'',
        '三个站的实际局部峰由负转正：i205为−0.01395243510，i206为+0.00671119561，i207为+0.02700309694。结合负定Hessian、小窗及有限域边界样本、实际iloc/MAC转移和原生向量比较，这支持原中性候选续算在i205附近停止是该连续分支的局部不稳定岛起始/消失，而非仅因Newton不佳。该结论限定于原400点离散和已跟踪M2p支；不代表整个联合谱稳定，不把原迭代失败写成无根。','',
        '|i|x (mm)|omega_peak|beta_peak|max g|alpha_r|Hessian特征值|梯度无穷范数|',
        '|---:|---:|---:|---:|---:|---:|---|---:|']
    for r in compact:lines.append('|%d|%.10f|%.11f|%.11f|%+.11g|%.11f|(%.7g, %.7g)|%.3g|'%(r['station'],r['x_mm'],r['omega'],r['beta'],r['max_growth'],r['alpha_r'],r['Hessian_eigenvalue_1'],r['Hessian_eigenvalue_2'],r['gradient_inf_norm']))
    lines += ['',
        '方法：从root已有实际同站seed出发解grad_(omega,beta) g=0；梯度/Hessian由同支9点中心差分场计算，Newton单次步长限制为1。确认Hessian负定后，将h逐级从.02减到.005并用更严PEP tol/NCV复核，梯度≤5.10e−7。各峰最后h减半中心g变化≤7.4e−14，修正量每轴约5e−6；这仅量化本次局部差分/求解器核验，不是法向网格误差界。',
        '判别依据：对充分光滑的实标量函数，驻点的负定Hessian是局部极大的充分条件；本次只是数值近似检验，不保证全域全局极大。已查MIT教材§7.3.1的Taylor展开/驻点分类：[Classification of stationary points](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node138.html)。另一次Stanford凸优化教材PDF读取超时，未将未读内容用于结论。','',
        '有限窗实算：每站原峰±.25的5×5小窗25点均接受且最高点在内部；随后从峰向四个象限外展，覆盖omega[24,34]、beta[30,40]，每站36条原始记录/25个去重点，其中16个边界点。所有采样点接受，未把缺失点当稳定。','',
        '|i|有限窗边界采样最高g|有限窗采样winner触边|小窗winner触边|','|---:|---:|---|---|']
    for r in state['finite_window_summaries']:lines.append('|%d|%.11g|%s|%s|'%(r['station'],r['boundary_winner']['growth'],r['winner_touches_window_boundary'],r['smallwindow_winner_touches_boundary']))
    lines += ['',
        '有限窗粗样本中未见比已解局部峰更高的值，边界均低于峰；这不是对未采样内部其他峰的严格排除，也未覆盖该窗外omega/beta或其他谱支。没有因非触边就宣称“整个联合中性边界闭合”。',
        '同支连续性：i205/206/207各峰原生向量与最后PEP raw index/alpha匹配误差<3.7e−13，原多项式绝对残差≤7.22e−13。三站峰均在j282、n≈.06442/.06490/.06538，描述性W99尺度内primitive模方约99.745%–99.747%。这些只是原始变量结构指标，不是能量或物理模态分类。',
        '将i205峰的真实状态实际沿iloc转移到206/207，再在同站沿omega/beta回到各自峰，全部MAC/残差门限接受，终点alpha与独立峰相差约3.1e−8（采样beta偏移1e−8）。跨站物理n映射、等primitive梯形积分MAC：205→206 %.9f、205→207 %.9f；坐标采用各站局部切/法向分量，非物理能量范数。'%(state['peak_shape_comparisons'][0]['normal_integral_equal_primitive_MAC'],state['peak_shape_comparisons'][1]['normal_integral_equal_primitive_MAC']),
        '解释限度：起始区由原离散壁点夹在x=%.10f至%.10f mm。未求非整数x的临界点，也没有实际在该点求解两条频率根的合并。原206/207接受的上频率中性点与本次正增长峰相容；负定局部二次模型和峰符号翻转可解释原g=0 Newton在205反复摆动，但不能推出此处其他起点/模式没有根。'%(bystation[205]['x_mm'],bystation[206]['x_mm']),
        '来源：`'+str(PROGRESS)+'`和`'+str(GAP)+'`的实际场/trace；两文件不修改。所有本次独立42次调用的summary路径/hash、原输入和冻结binary完整性核对都在JSON中。实际总求解调用耗时求和%.3fs，field总记录%d；0失败field、0最终拒绝、0 PETSc错误、0未使用参数。sample_only task33的终态no bracketed root仅为不执行根搜索的控制结果，不能解读成没有中性根。'%(sum(r['seconds'] for r in inventory),sum(r['field_records'] for r in inventory)),
        '本轮未开展该M2p峰法向加密、基流独立网格或物理因果/辐射分类，故峰符号首先是现有400点模型的观察；没有套用此前M3near的法向误差结果。',
        '结束审查：局部峰/小窗/有限边界、形态与实际跨站连续性已完成，原求解/继续失败证据全部保留，未把未知改为稳定。一个复现脚本及全部配置/原生向量/日志、JSON/CSV与工作日志均保留；没有改原H5或源码，没有本任务生成且不再需要的尝试性代码，无需删除用户/其他任务文件。']
    paragraphs=[];table=[]
    for line in lines:
        if line.startswith('|'):table.append(line);continue
        if table:paragraphs.append('\n'.join(table));table=[]
        if line:paragraphs.append(line)
    if table:paragraphs.append('\n'.join(table))
    (ROOT/'outputs/local_mode_audit_M2p_onset_说明.md').write_text('\n\n'.join(paragraphs)+'\n',encoding='utf-8')
    log('结束审查：42个实际短调用、3负定局部峰、75小窗+108有限窗记录与实际跨站/MAC终点核对完成。源H5/binary对原manifest核验通过；原root进度/失败文档未写入。205负峰/206正峰仅作此局部支起始区证据，保留未采样/其他支/法向未核验限制。无需清理的尝试代码。')

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stations',nargs='+',type=int,default=[205,206,207]);parser.add_argument('--boundary',action='store_true');parser.add_argument('--shapes',action='store_true');parser.add_argument('--report',action='store_true');args=parser.parse_args()
    assert driver.sha(BIN)==HASH
    if args.report:finalize(json.loads(OUT.read_text()));return
    if not LOG.exists():LOG.write_text('# M2p 上游中性候选停止点局部峰诊断\n\n日期2026-10-23，Asia/Shanghai。\n\n目标/评估：root FD在i205未收敛不能解释成无根；独立寻找该实际连续分支g(omega,beta)局部最大。固定有限窗omega[24,34]、beta[30,40]，优先205/206，便宜时207。严格保持真实同站seed及来源，负alpha_r根不删除。\n\n策略：从原9点实际场解grad g=0，用负定Hessian及减半差分/tight核验；另以围绕真峰的小窗和有限域边界样本检查形态、触边/未知。单线程、有限短调用，不重新搜全谱。完成标准：实际峰growth符号、梯度/Hessian/来源、窗口样本及失败、限定的起始区解释和结束审查。\n',encoding='utf-8')
    progress=json.loads(PROGRESS.read_text());gap=json.loads(GAP.read_text())
    same={r['station']:dict(r['point'],source_path=str(PROGRESS)) for r in progress['records'] if r['status']=='candidate_converged'}
    same[205]=dict(gap['result']['trace'][-1]['center'],source_path=str(GAP),source_run=gap['result']['trace'][-1]['name'])
    state=json.loads(OUT.read_text()) if OUT.exists() else dict(created=stamp(),bounds=BOUNDS.tolist(),solver=str(BIN),solver_sha256=HASH,
        provenance={str(p):driver.sha(p) for p in [PROGRESS,GAP]},maxima=[],surveys=[],qualification='Finite local single-branch diagnosis, not entire joint-spectrum stability.')
    for station in args.stations:
        existing=next((r for r in state['maxima'] if r['station']==station),None)
        if existing:continue
        log('从实际同站seed开始局部峰：'+str(station)+' '+json.dumps(same[station],ensure_ascii=False))
        r=maximize(station,same[station]);state['maxima'].append(r);save(state)
        print('M2P MAXIMUM '+json.dumps({k:v for k,v in r.items() if k not in ['trace','verification','seed']},indent=2),flush=True)
        log('局部峰结果：'+json.dumps({k:v for k,v in r.items() if k not in ['trace','verification','seed']},ensure_ascii=False))
        if r['status']!='local_maximum_candidate':continue
        p=r['point'];small=survey(station,p,p['omega_r']-.25,p['omega_r']+.25,5,p['beta_r']-.25,p['beta_r']+.25,5,'smallwindow')
        state['surveys'].append(dict(station=station,kind='smallwindow',**small));save(state)
    if args.boundary:
        for maximum in state['maxima']:
            station=maximum['station']
            if station not in args.stations or maximum['status']!='local_maximum_candidate':continue
            p=maximum['point'];w=p['omega_r'];b=p['beta_r']
            # Four outward arms visit both omega boundaries and beta boundaries;
            # include the exact seed at each arm's first beta/omega point.
            for wl,wh,nw,rw,wlabel in [(24.,w,3,True,'wdown'),(w,34.,3,False,'wup')]:
                for bl,bh,nb,rb,blabel in [(30.,b,3,True,'bleft'),(b,40.,3,False,'bright')]:
                    label='finitewindow_'+wlabel+'_'+blabel
                    if any(s['station']==station and s['kind']==label for s in state['surveys']):continue
                    s=survey(station,p,wl,wh,nw,bl,bh,nb,label,rw,rb)
                    state['surveys'].append(dict(station=station,kind=label,**s));save(state)
                    print('M2P WINDOW '+json.dumps(dict(station=station,kind=label,complete=s['complete'],winner=s['winner'])),flush=True)
    if args.shapes:shapes(state,args.stations)
    assert driver.sha(BIN)==HASH

if __name__=='__main__':main()
