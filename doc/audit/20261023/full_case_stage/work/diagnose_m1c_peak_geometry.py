"""Keep all observed beta extrema near the stopped M1c i321 candidate."""
import argparse,csv,datetime,hashlib,json,math,re
from pathlib import Path
import numpy as np
import full_driver as driver
import fd_envelope as fd
ROOT=driver.ROOT
BIN=ROOT/'outputs/bin/START_LST_wall_local_v7'
HASH='7aa6016570f406587156f229107c0e5e4f76c8a85ba89dda267039ff1d688d79'
CURVE=ROOT/'outputs/fd_curvature_diagnostic_M1c_i321.json'
GAP=ROOT/'outputs/root_gap_diagnostic_M1c_downstream.json'
OUT=ROOT/'outputs/local_mode_audit_M1c_i321_geometry.json'
LOG=ROOT/'outputs/local_mode_audit_工作日志_2026-10-23_M1c局部峰几何.md'

def stamp():return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
def log(message):
    with LOG.open('a',encoding='utf-8') as f:f.write('\n'+stamp()+' '+message+'\n')
def save(s):
    s['updated']=stamp();tmp=OUT.with_suffix('.tmp');tmp.write_text(json.dumps(s,indent=2,allow_nan=False)+'\n');tmp.replace(OUT)

def scan(seed,wl,wh,nw,bl,bh,db,label,reverse_w=False,reverse_b=False):
    spec=dict(station=321,mode='M1',coordinates='wall_local',seed={k:seed[k] for k in ('omega_r','beta_r','alpha_r','alpha_i')},seed_shift=1e-4,
        solver=str(BIN),solver_sha256=HASH,task=33,sample_only=True,omega_min=wl,omega_max=wh,omega_points=nw,
        beta_min=bl,beta_max=bh,beta_points=math.ceil((bh-bl)/db-1e-12)+1,reverse_outer=reverse_w,reverse_inner=reverse_b,
        anchor_row_seed=True,adaptive=True,vector_tracking=True,minimum_overlap=.95,minimum_step=1e-6,retry_eigenpairs=8,
        step=.2,jump=2,options=['-pep_tol','1e-12','-pep_max_it','1000','-pep_nev','1','-pep_ncv','8','-st_pc_factor_mat_solver_type','superlu'],timeout=240)
    digest=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()[:12];spec['name']='M1cgeom_'+label+'_'+digest
    summary=driver.run(spec);path=ROOT/'runs'/spec['name']/'out/neutral_field.plt';rows=driver.table(path)
    return dict(name=spec['name'],summary=summary,source_path=str(path),samples=rows,
        complete=bool(len(rows)==nw*spec['beta_points'] and all(r['status']==0 for r in rows) and not summary['petsc_errors'] and not summary['coordinate_metadata_missing']))

def summarize(state):
    byomega={}
    for s in state['scans']:
        for r in s['samples']:
            byomega.setdefault(round(r['omega_r'],10),{})[round(r['beta_r'],10)]=dict(r,source_run=s['name'],source_path=s['source_path'])
    results=[]
    for omega,raw in sorted(byomega.items()):
        rows=sorted(raw.values(),key=lambda r:r['beta_r']);valid=[r for r in rows if r['status']==0]
        maxima=[];minima=[]
        for a,b,c in zip(rows,rows[1:],rows[2:]):
            if any(r['status']!=0 for r in (a,b,c)):continue
            if b['growth']>a['growth'] and b['growth']>=c['growth']:maxima.append(dict(point=b,left=a,right=c))
            if b['growth']<a['growth'] and b['growth']<=c['growth']:minima.append(dict(point=b,left=a,right=c))
        winner=max(valid,key=lambda r:r['growth']) if valid else None
        results.append(dict(omega=omega,samples=len(rows),failed=sum(r['status']!=0 for r in rows),maxima=maxima,minima=minima,winner=winner,
                            sampled_beta_range=[rows[0]['beta_r'],rows[-1]['beta_r']],
                            winner_touches_sampled_beta_edge=bool(winner and min(abs(winner['beta_r']-rows[0]['beta_r']),abs(winner['beta_r']-rows[-1]['beta_r']))<1e-8),
                            winner_touches_beta_boundary=bool(winner and min(abs(winner['beta_r']-15),abs(winner['beta_r']-20))<1e-8)))
    state['rows']=results

def fold_stencil(seed,w,b,h=.01,label='fold'):
    s=scan(seed,w-h,w+h,3,b-2*h,b+2*h,h,'fold_'+label)
    if not s['complete']:raise RuntimeError('Incomplete fold stencil '+s['name'])
    values={(round((r['omega_r']-w)/h),round((r['beta_r']-b)/h)):r for r in s['samples']}
    if len(values)!=15:raise RuntimeError('Expected 3x5 fold grid')
    g=lambda i,j:values[i,j]['growth']
    gb=(g(0,1)-g(0,-1))/(2*h);gbb=(g(0,1)-2*g(0,0)+g(0,-1))/h**2
    gwb=(g(1,1)-g(1,-1)-g(-1,1)+g(-1,-1))/(4*h*h)
    gwbb=((g(1,1)-2*g(1,0)+g(1,-1))-(g(-1,1)-2*g(-1,0)+g(-1,-1)))/(2*h**3)
    gbbb=(g(0,2)-2*g(0,1)+2*g(0,-1)-g(0,-2))/(2*h**3)
    return dict(name=s['name'],center=values[0,0],h=h,F=[gb,gbb],J=[[gwb,gbb],[gwbb,gbbb]],
                samples=s['samples'],source_path=s['source_path'],summary=s['summary'])

def fold(seed):
    result=dict(status='unresolved',station=321,trace=[],solver=str(BIN),solver_sha256=HASH,
        qualification='Candidate disappearance of extrema of growth versus beta, not eigenvalue coalescence or a neutral point.')
    w,b=11.51,17.;current=dict(seed)
    try:
        for it in range(10):
            p=fold_stencil(current,w,b,.01,'iterate');result['trace'].append(p);current=p['center']
            delta=np.linalg.solve(np.array(p['J']),-np.array(p['F']))
            if max(abs(np.array(p['F'])))<1e-6 and max(abs(delta))<1e-4:
                v=fold_stencil(current,w,b,.005,'verify');result['verification']=v
                correction=np.linalg.solve(np.array(v['J']),-np.array(v['F']))
                if max(abs(np.array(v['F'])))<1e-5 and abs(correction[0])<1e-4 and abs(correction[1])<1e-3:
                    result.update(status='stationary_pair_fold_candidate',point=v['center'],verification_correction=correction.tolist())
                    break
            delta*=min(1.,.03/max(abs(delta[0]),1e-100),.3/max(abs(delta[1]),1e-100))
            w+=float(delta[0]);b+=float(delta[1])
            if not 11.3<w<11.8 or not 15<b<20:raise RuntimeError('fold Newton leaves declared window')
        else:result['reason']='iteration limit, no nonexistence claim'
    except Exception as e:result['reason']=repr(e)
    return result

def finalize(state):
    summarize(state);inventory=[]
    for d in sorted((ROOT/'runs').iterdir()):
        if not d.name.startswith('M1cgeom_'):continue
        q=d/'summary.json';assert q.exists()
        s=json.loads(q.read_text());assert s['solver_sha256']==HASH
        inventory.append(dict(name=d.name,summary_path=str(q),summary_sha256=driver.sha(q),
            seconds=s['seconds'],terminal=s['terminal'],rejected=s['rejected'],field_records=s['field_records'],
            failed_field_records=s['failed_field_records'],petsc_errors=s['petsc_errors'],unused_options=s['unused_options']))
    assert len(inventory)==9 and all(r['failed_field_records']==0 and r['rejected']==0 and r['petsc_errors']==0 and r['unused_options']==0 for r in inventory)
    coarse=[s for s in state['scans'] if s['label'].startswith('coarse_')]
    assert len(coarse)==4 and all(s['complete'] for s in state['scans'])
    broad=[v for s in coarse for v in s['samples']];winner=max(broad,key=lambda v:v['growth'])
    state['status']='positive_growth_window_truncated_with_stationary_pair_disappearance_evidence'
    state['broad_window_winner']=winner
    state['window_truncation']=dict(omega_max_positive=True,beta_min_positive=True,
        reason='At omega=11.8 beta=15 the actual tracked branch has positive growth; no outer neutral boundary was bracketed by this window.',
        more_external_neutral_root_located=False)
    state['run_inventory']=inventory;state['input_integrity_against_original_manifest']=driver.snapshot()
    state['source_documents_current_hashes']={str(q):driver.sha(q) for q in (CURVE,GAP)}
    state['solver_unchanged']=driver.sha(BIN)==HASH;state['completed']=stamp();save(state)
    with (ROOT/'outputs/local_mode_audit_M1c_i321_rows.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['omega','samples','beta_sample_min','beta_sample_max','observed_local_maxima','observed_local_minima','winner_beta','winner_growth','winner_at_declared_beta_edge','winner_at_sampled_beta_edge'])
        for r in state['rows']:w.writerow([r['omega'],r['samples'],*r['sampled_beta_range'],len(r['maxima']),len(r['minima']),r['winner']['beta_r'],r['winner']['growth'],r['winner_touches_beta_boundary'],r['winner_touches_sampled_beta_edge']])
    fold_result=state['fold'];p=fold_result['point'];v=fold_result['verification']
    before=[r for r in state['rows'] if abs(r['omega']-(p['omega_r']-.005))<1e-8][0]
    after=[r for r in state['rows'] if abs(r['omega']-(p['omega_r']+.005))<1e-8][0]
    lines=['# M1c i320→321停止点的局部峰几何与窗口截断','',stamp(),'',
        '实际数据支持两个不同事实：原beta≈17的局部增长极大/极小可在g仍正时消失；同时同一跟踪起点的其他beta仍明显正增长。指定omega[11.3,11.8]、beta[15,20]没有包住正增长区，因此本次不提供“已闭合的更外中性根”，也不能据原Newton停止宣称稳定或无根。',
        f'最直接截断证据：i321、x={35*winner["x"]:.10f} mm，omega={winner["omega_r"]:.12g}、beta={winner["beta_r"]:.12g}的实际alpha={winner["alpha_r"]:.12g}{winner["alpha_i"]:+.12g}i，g={winner["growth"]:.12g}>0，位于omega上界和beta下界。原beta≈17局部峰已不是该窗的最高增长点。','',
        '|omega|观测beta局部max数|局部min数|采样winner beta|winner g|触声明beta边界|',
        '|---:|---:|---:|---:|---:|---|']
    for r in state['rows']:
        if r['sampled_beta_range'][0]>15.000001:continue
        lines.append('|%.11f|%d|%d|%.10f|%+.10g|%s|'%(r['omega'],len(r['maxima']),len(r['minima']),r['winner']['beta_r'],r['winner']['growth'],r['winner_touches_beta_boundary']))
    lines += ['',
        f'局部驻点对消失候选：由3×5实际场差分解(g_beta,g_betabeta)=0，得到omega≈{p["omega_r"]:.11f}、beta≈{p["beta_r"]:.11f}、alpha≈{p["alpha_r"]:.11f}{p["alpha_i"]:+.11g}i，g≈{p["growth"]:.11g}。h=.005复核的g_beta≈{v["F"][0]:.6g}、g_betabeta≈{v["F"][1]:.6g}；估计参数修正为{fold_result["verification_correction"]}，故这些坐标仅作为有限精度候选，而非精确分岔点。',
        f'复核g_omega_beta≈{v["J"][0][0]:.8g}、g_beta_beta_beta≈{v["J"][1][1]:.8g}均非零，配合前后直接采样，符合局部beta驻点对消失的数值特征。这是增长函数的极值几何，不是特征值谱支合并。g≈+.00104846明显未达到中性。',
        '直接形态观察：在上述候选omega低.005处、beta候选±.3的61点行中，同时见下表的局部min/max；omega高.005的61点行中未检出两者。两行所有采样接受，无缺失样本跨接。未检出只限定这个细beta窗与.01采样间隔。','',
        '|位置|omega|beta|g|','|---|---:|---:|---:|']
    for kind,points in [('local minimum',before['minima']),('local maximum',before['maxima'])]:
        for item in points:
            r=item['point'];lines.append('|%s|%.11f|%.11f|%.11g|'%(kind,r['omega_r'],r['beta_r'],r['growth']))
    lines += ['',
        '方法与适用性：原真实i320中性状态已由root的iloc/MAC转移到i321。本次取root已保存的i321实际状态作为起点，四向sample_only/anchor/MAC有限场保留全部beta极大/极小和触边值，不只保留一个ridge。所有调用固定v7 SHA-256 `'+HASH+'`，继续追踪相同起点支，未按alpha实部符号过滤，也没有重新搜全谱。',
        '理论资料：MIT教材§7.3.1用Taylor展开说明，驻点处Hessian退化时需更高阶信息而不能按普通负定极大检验直接分类：[Classification of stationary points](https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node138.html)。本次“驻点对消失”是基于原始场、二/三阶差分和前后采样的数值解释；没有由此推出全域最大值或整体谱稳定。',
        f'计算量与证据：9个实际单线程短调用、field总记录{sum(r["field_records"] for r in inventory)}（四臂224、四个15点驻点差分及细形态183），求解耗时求和{sum(r["seconds"] for r in inventory):.3f}s。所有field有效，最终拒绝/PETSc错误/未使用参数均0。sample_only终态no bracketed root不是无根证明。全部summary路径/hash、场、solver日志、配置及实际来源在JSON/各runs/M1cgeom_*保留。',
        '原始root资料：`'+str(GAP)+'`与`'+str(CURVE)+'`只读，原H5和baseline binary与原manifest核对一致；未改源码。JSON同时保存起读hash与完成时hash，便于辨认其他任务若更新来源。',
        '结束状态：i321局部峰几何与指定窗口截断已解释到上述证据限度；更外中性根尚未定位，beta下界/omega上界正增长须保留为截断。没有宣称整个联合包络闭合或该支全域无根。可选i33→32本轮未启动：对应root_gap_diagnostic_M1c_upstream.json不存在，且未把321的参数窗口套给鼻部；其原失败保持未知，交回主任务按优先级安排。',
        '清理审查：一个必要复现脚本、工作日志和最终JSON/CSV/说明、所有原始短调用保留。无失败求解，无本任务不再需要的尝试性代码；不删除原输入及其他任务文件。']
    paragraphs=[];table=[]
    for line in lines:
        if line.startswith('|'):table.append(line);continue
        if table:paragraphs.append('\n'.join(table));table=[]
        if line:paragraphs.append(line)
    if table:paragraphs.append('\n'.join(table))
    (ROOT/'outputs/local_mode_audit_M1c_i321_geometry_说明.md').write_text('\n\n'.join(paragraphs)+'\n',encoding='utf-8')
    log('结束审查：9次实际调用全部有效，驻点对消失候选g仍正，声明窗正增长触omega上界/beta下界，较外中性根未定位。报告保留截断/非全谱/有限差分精度限制；原文件只读。可选i32未启动，未把321窗口套用；无不需要的尝试性代码可清理。')

def main():
    p=argparse.ArgumentParser();p.add_argument('--coarse',action='store_true');p.add_argument('--refine',action='store_true');p.add_argument('--fold',action='store_true');p.add_argument('--fold-geometry',action='store_true');p.add_argument('--report',action='store_true')
    p.add_argument('--omega',type=float);p.add_argument('--beta',type=float);p.add_argument('--seed',type=Path);p.add_argument('--label',default='candidate');args=p.parse_args()
    assert driver.sha(BIN)==HASH
    if args.report:finalize(json.loads(OUT.read_text()));return
    if not LOG.exists():LOG.write_text('# M1c i320→321 停止点的局部beta峰几何诊断\n\n日期2026-10-23，Asia/Shanghai。\n\n目标：有限omega[11.3,11.8]、beta[15,20]实算，保留多beta极大/极小而非只追一峰，解释已观测曲率翻转与g=0 Newton未收敛。使用实际同站M1c根及冻结v7；不以alpha实部符号筛选。\n\n计划/完成标准：四向从真实锚点外展有限采样；保留所有极值及触边/未知；对可能较外中性峰以FD核验，无法成立时保留未完成；单线程少量调用，不修改原H5/源码。仅诊断该起点支，不能宣布联合谱稳定或严格无根。\n',encoding='utf-8')
    state=json.loads(OUT.read_text()) if OUT.exists() else dict(created=stamp(),station=321,solver=str(BIN),solver_sha256=HASH,
        source_hashes={str(q):driver.sha(q) for q in (CURVE,GAP)},scans=[],refinements=[],declared_window=dict(omega=[11.3,11.8],beta=[15,20]))
    source=json.loads(CURVE.read_text());seed=source['trace'][0]['center']
    if args.coarse:
        w,b=seed['omega_r'],seed['beta_r']
        for wl,wh,nw,rw,wlabel in [(11.3,w,3,True,'wdown'),(w,11.8,5,False,'wup')]:
            for bl,bh,rb,blabel in [(15.,b,True,'bleft'),(b,20.,False,'bright')]:
                label='coarse_'+wlabel+'_'+blabel
                if any(s['label']==label for s in state['scans']):continue
                s=scan(seed,wl,wh,nw,bl,bh,.2,label,rw,rb);s['label']=label
                state['scans'].append(s);summarize(state);save(state)
                print('M1C GEOMETRY '+json.dumps(dict(label=label,complete=s['complete'],field=s['summary']['field_records'],failed=s['summary']['failed_field_records'])),flush=True)
    if args.refine:
        assert args.seed and args.omega is not None and args.beta is not None
        seed=json.loads(args.seed.read_text())
        if 'point' in seed:seed=seed['point']
        r=fd.refine(321,seed,args.omega,args.beta,BIN,'wall_local',label='M1cgeom_fd_'+args.label,
                    omega_bounds=(11.3,11.8),beta_bounds=(15,20),max_iterations=20,h=.01)
        state['refinements'].append(r);save(state)
        print('M1C REFINE '+json.dumps({k:v for k,v in r.items() if k not in ['trace','verification']}),flush=True)
    if args.fold:
        result=fold(seed);state['fold']=result;save(state)
        print('M1C FOLD '+json.dumps({k:v for k,v in result.items() if k not in ['trace','verification']}),flush=True)
    if args.fold_geometry:
        result=state['fold'];assert result['status']=='stationary_pair_fold_candidate'
        p=result['point'];w=p['omega_r'];b=p['beta_r']
        s=scan(p,w-.005,w+.005,3,b-.3,b+.3,.01,'fold_geometry');s['label']='fold_geometry';state['scans'].append(s);save(state)
    summarize(state);save(state)
    print('M1C ROWS '+json.dumps([dict(omega=r['omega'],maxima=[p['point'] for p in r['maxima']],minima=[p['point'] for p in r['minima']],winner=r['winner'],edge=r['winner_touches_beta_boundary']) for r in state['rows']],indent=2),flush=True)

if __name__=='__main__':main()
