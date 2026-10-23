"""Independent M3near normal-grid/domain tests; never writes original HDF5."""
import concurrent.futures as cf
import configparser,csv,datetime,hashlib,json,math,os,re,struct,subprocess,sys,time
from pathlib import Path
import h5py,numpy as np,scipy
from scipy.interpolate import PchipInterpolator,CubicSpline

HERE=Path(__file__).resolve().parent
RUN=HERE.parent.parent
sys.path.insert(0,str(RUN/'work'))
import full_driver as driver
SOURCE=RUN.parent/'M1/in'
BIN=RUN/'outputs/bin/START_LST_wall_local_v5'
HASH='528d6d6a1b7aab34c73a70b09721295089ba293b2ef2c9bec173aba8e589c266'
TARGET=complex(35.93408168168585,-4.023207897824065)
FLOAT=r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?'
OUT=HERE/'outputs'
OUT.mkdir(exist_ok=True)
LOG=OUT/'工作日志_2026-10-22_M3near法向敏感性.md'

def stamp():return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
def log(message):
    with LOG.open('a',encoding='utf-8') as f:f.write('\n'+stamp()+' '+message+'\n')

def vectors(path,jn):
    data=path.read_bytes();offset=0;vv=[]
    while offset<len(data):
        cls,n=struct.unpack_from('>ii',data,offset);offset+=8
        assert cls==1211214 and n==jn*5,(cls,n,jn)
        end=offset+n*16;assert end<=len(data)
        vv.append(np.frombuffer(data[offset:end],dtype='>c16').astype(np.complex128).reshape(jn,5));offset=end
    return vv

def geometry(grid):
    g=grid[:,118,:,0];normal=g[0]/np.linalg.norm(g[0]);n=(g-g[0])@normal
    assert np.all(np.diff(n)>0)
    return n

def generate_variant(name,method,jn,source_hashes):
    d=HERE/'inputs'/name;d.mkdir(parents=True,exist_ok=True)
    manifest=d/'generation.json'
    if manifest.exists():
        result=json.loads(manifest.read_text())
        assert result['source_hashes']==source_hashes
        assert all(driver.sha(Path(x['path']))==x['sha256'] for x in result['files'])
        return result
    eta=np.arange(400,dtype=float)
    positions=np.linspace(0,399,jn) if method!='crop' else np.arange(jn,dtype=float)
    info=dict(name=name,method=method,jn=jn,in_=841,eta_min=float(positions[0]),eta_max=float(positions[-1]),
              created=stamp(),source_hashes=source_hashes,scipy_version=scipy.__version__,numpy_version=np.__version__,files=[])
    for filename,dataset in [('grid_vertical2d_N.h5','grid'),('flow_vertical2d_N.h5','baseflow')]:
        with h5py.File(SOURCE/filename,'r') as f:
            original=f[dataset][:];attrs=dict(f[dataset].attrs)
        if method=='crop':data=original[:jn].copy()
        elif method=='PCHIP':data=PchipInterpolator(eta,original,axis=0,extrapolate=False)(positions)
        elif method=='CubicSpline':data=CubicSpline(eta,original,axis=0,bc_type='not-a-knot',extrapolate=False)(positions)
        else:raise ValueError(method)
        assert data.shape==(jn,841,original.shape[2],2) and np.isfinite(data).all()
        assert np.all(data[:,:,:,1]==0)
        if dataset=='grid':n=geometry(data);info['nmax_N']=float(n[-1]);info['nmax_mm']=35*float(n[-1]);info['min_delta_n_N']=float(np.diff(n).min())
        if dataset=='baseflow':
            assert np.all(data[:,:,0,0]>0) and np.all(data[:,:,4,0]>0)
            info['minimum_density']=float(data[:,:,0,0].min());info['minimum_temperature']=float(data[:,:,4,0].min())
        # Every original source node is exactly retained in nested refinements.
        if method!='crop':
            stride=(jn-1)//399
            assert (jn-1)%399==0
            info[dataset+'_source_node_max_difference']=float(np.max(abs(data[::stride]-original)))
        output=d/filename
        assert not output.exists(),'refuse overwriting a partial derived input'
        with h5py.File(output,'w') as f:
            ds=f.create_dataset(dataset,data=data,compression='gzip',compression_opts=1)
            for k,v in attrs.items():ds.attrs[k]=v
            ds.attrs['plot3d_shape']=np.array([841,jn],dtype=np.int64)
            ds.attrs['fortran_logical_shape']=np.array([2,data.shape[2],841,jn],dtype=np.int64)
            ds.attrs['codex_derivation']=method+' on the original uniform computational j coordinate; no extrapolation'
            ds.attrs['codex_original_sha256']=source_hashes[filename]
        info['files'].append(dict(path=str(output),dataset=dataset,sha256=driver.sha(output),shape=list(data.shape),bytes=output.stat().st_size))
        del original,data
    manifest.write_text(json.dumps(info,indent=2)+'\n')
    log('生成独立输入 '+name+' '+json.dumps(info,ensure_ascii=False))
    return info

def solve(name,input_dir,jn,method):
    d=HERE/'runs'/name;d.mkdir(parents=True,exist_ok=True);result_path=d/'result.json'
    if result_path.exists():return json.loads(result_path.read_text())
    assert not (d/'solver.log').exists(),'partial existing run retained, no overwrite'
    cfg=configparser.ConfigParser();cfg.read(RUN.parent/'M1/config.ini')
    cfg['filenames']['dir']=str(input_dir)+'/'
    cfg['domain']['jn']=str(jn);cfg['domain']['iloc']='119'
    cfg['lst']['profile_coordinates']='wall_local';cfg['lst']['alpha']=f'{TARGET.real+1e-4:.16g}, {TARGET.imag:.16g}'
    cfg['lst']['omega']='0, 0';cfg['lst']['beta']='-7, 0'
    cfg['neutral']={k:str(v) for k,v in dict(enabled='T',track_initial_target='T',task=33,
        outer_parameter='omega',inner_parameter='beta',outer_min=0.,outer_max=0.,outer_points=1,
        inner_min=-7.-1e-8,inner_max=-7.+1e-8,scan_points=2,growth_tol=1e-8,parameter_tol=1e-7,
        max_iterations=70,eigen_error_tol=1e-7,max_eigen_jump=2,continuation_step=.005,
        output_prefix='out/neutral',sample_only='T',adaptive_steps='T',minimum_step=1e-6,
        vector_tracking='T',minimum_overlap=.95,retry_eigenpairs=0).items()}
    with (d/'config.ini').open('w') as f:cfg.write(f)
    cmd=[str(BIN),'-f','config.ini','-pep_tol','1e-13','-pep_max_it','1000','-pep_nev','1','-pep_ncv','8',
         '-st_pc_factor_mat_solver_type','superlu','-pep_view_vectors','binary:out/all_vectors.petsc',
         '-pep_view_values','ascii:out/all_values.txt','-options_left']
    started=stamp();t=time.time()
    with (d/'solver.log').open('w') as f:
        f.write('COMMAND '+json.dumps(cmd)+'\n');f.flush()
        try:p=subprocess.run(cmd,cwd=d,env=driver.ENV,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,timeout=180);code=p.returncode
        except subprocess.TimeoutExpired:code=124
    text=(d/'solver.log').read_text(errors='replace')
    blocks=[];block=[]
    for line in text.splitlines():
        m=re.search(r'^\s*(\d+)\s+\(\s*('+FLOAT+r')\s*,\s*('+FLOAT+r')\s*\)\s+The error is\s+('+FLOAT+r')',line)
        if not m:continue
        index=int(m[1])
        if index==0 and block:blocks.append(block);block=[]
        block.append(dict(index=index,alpha_r=float(m[2].replace('D','E')),alpha_i=float(m[3].replace('D','E')),absolute_polynomial_residual=float(m[4].replace('D','E'))))
    if block:blocks.append(block)
    pairs=blocks[-1] if blocks else []
    fields=driver.table(d/'out/neutral_field.plt')
    runlog=(d/'out/neutral_run.log').read_text() if (d/'out/neutral_run.log').exists() else ''
    m=re.search(r'evaluations=\s*(\d+)\s+solves=\s*(\d+)\s+rejected=\s*(\d+)\s+result=(.+)',runlog)
    overlaps=[float(v.replace('D','E')) for v in re.findall(r'vector_overlap=\s*('+FLOAT+r')',runlog)]
    result=dict(name=name,method=method,jn=jn,input_dir=str(input_dir),started=started,seconds=round(time.time()-t,3),
                solver=str(BIN),solver_sha256=driver.sha(BIN),command=cmd,exit=code,
                petsc_errors=text.count('PETSC ERROR'),unused_options=text.count('Option left:'),
                terminal=m[4].strip() if m else 'MISSING',rejected=int(m[3]) if m else None,
                field_records=len(fields),failed_field_records=sum(r['status']!=0 for r in fields),
                minimum_tiny_beta_MAC=min(overlaps) if overlaps else None,pairs=pairs,passed=False)
    path=d/'out/all_vectors.petsc'
    if pairs and path.exists():
        vv=vectors(path,jn)
        result['raw_vector_count']=len(vv);result['raw_vector_index_association_verified']=len(pairs)==len(vv)
        k=min(range(len(pairs)),key=lambda k:abs(complex(pairs[k]['alpha_r'],pairs[k]['alpha_i'])-TARGET))
        a=complex(pairs[k]['alpha_r'],pairs[k]['alpha_i'])
        with h5py.File(input_dir/'grid_vertical2d_N.h5','r') as f:g=f['grid'][:]
        with h5py.File(input_dir/'flow_vertical2d_N.h5','r') as f:q=f['baseflow'][:,118,:,0]
        n=geometry(g);vec=vv[k];amplitude=np.sum(abs(vec)**2,axis=1);peak=int(np.argmax(amplitude))
        j99=int(np.flatnonzero(q[:,3]>=.99*q[-1,3])[0])
        result.update(alpha_r=a.real,alpha_i=a.imag,growth=-a.imag,delta_alpha_to_original=abs(a-TARGET),
                      absolute_polynomial_residual=pairs[k]['absolute_polynomial_residual'],raw_index=k,
                      native_vector_file=str(path),vector_sha256=driver.sha(path),
                      nmax_N=float(n[-1]),peak_j=peak+1,peak_n_N=float(n[peak]),W99_n_N=float(n[j99]),
                      interpretation='Unweighted primitive vector and descriptive W99 only; no physical identity certificate.')
        result['passed']=bool(code in (0,1) and result['petsc_errors']==0 and result['unused_options']==0 and
             result['rejected']==0 and len(fields)==2 and result['failed_field_records']==0 and
             result['raw_vector_index_association_verified'] and result['absolute_polynomial_residual']<=1e-7 and
             result['delta_alpha_to_original']<2 and result['minimum_tiny_beta_MAC']>=.95)
    result_path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('NORMAL '+json.dumps({k:result.get(k) for k in ('name','jn','seconds','passed','alpha_r','alpha_i','delta_alpha_to_original','absolute_polynomial_residual')}),flush=True)
    return result

def compare(reference,other):
    def data(r):
        with h5py.File(Path(r['input_dir'])/'grid_vertical2d_N.h5','r') as f:g=f['grid'][:]
        return geometry(g),vectors(Path(r['native_vector_file']),r['jn'])[r['raw_index']]
    n0,v0=data(reference);n1,v1=data(other)
    mask=n0<=n1[-1]+1e-12;x=n0[mask];a=v0[mask];b=CubicSpline(n1,v1,axis=0,extrapolate=False)(x)
    weights=np.empty(len(x));weights[1:-1]=(x[2:]-x[:-2])/2;weights[0]=(x[1]-x[0])/2;weights[-1]=(x[-1]-x[-2])/2
    weights=weights[:,None]
    mac=float(abs(np.sum(weights*np.conj(a)*b))**2/(np.sum(weights*abs(a)**2)*np.sum(weights*abs(b)**2)))
    euclidean=float(abs(np.vdot(a,b))**2/(np.vdot(a,a).real*np.vdot(b,b).real))
    z0=complex(reference['alpha_r'],reference['alpha_i']);z1=complex(other['alpha_r'],other['alpha_i'])
    return dict(reference=reference['name'],other=other['name'],delta_alpha=abs(z1-z0),
                relative_complex_alpha_difference=abs(z1-z0)/abs(z0),
                growth_difference=other['growth']-reference['growth'],relative_growth_difference=(other['growth']-reference['growth'])/abs(reference['growth']),
                mapped_to_original_grid_points=len(x),normal_integral_primitive_MAC=mac,mapped_Euclidean_MAC=euclidean,
                mapping='CubicSpline of complex vector on physical normal distance, no extrapolation; integral uses trapezoid weights and equal primitive-component weights, not physical energy.')

def finalize(report):
    order=['original400','crop340','crop360','pchip799','pchip1597','cubic799']
    rows=sorted(report['results'],key=lambda r:order.index(r['name']))
    assert report['source_unchanged'] and report['solver_unchanged']
    assert len(rows)==6 and all(r['passed'] for r in rows)
    byname={r['name']:r for r in rows}
    comparisons={(r['reference'],r['other']):r for r in report['comparisons']}
    fine=comparisons['pchip799','pchip1597'];allfine=comparisons['original400','pchip1597']
    crop=max(comparisons['original400',name]['delta_alpha'] for name in ['crop340','crop360'])
    method=comparisons['pchip799','cubic799']
    notes=['# M3near 代表点法向离散与外边界敏感性验证','',stamp(),'',
        '本次实际检查支持：在i119、omega=0、beta=−7这个代表点，M3near近壁根在两次裁域后保持同一向量与特征值；法向400→799→1597插值加密后增长率变化缩小。此结论限定于现有CFD输入重构及此根、此参数点，不能替代CFD网格独立、整条中性边界收敛或物理身份/因果条件验证。','',
        '原400点基准只读M1/in；五组派生H5均存于本目录inputs/。保留原841个i站，避免把流向导数/几何邻域改成另一问题。jn由实际config明确设置，native PEP向量长度逐一核对为5*jn，从而确认不同离散实际进入求解。冻结v5 SHA-256为 `'+HASH+'`，没有修改源码或v7。','',
        '|案例|jn|nmax（无量纲）|alpha_r|alpha_i|sigma=−alpha_i|原多项式绝对残差|',
        '|---|---:|---:|---:|---:|---:|---:|']
    for r in rows:notes.append('|%s|%d|%.9g|%.12f|%.12f|%.12f|%.3g|'%(r['name'],r['jn'],r['nmax_N'],r['alpha_r'],r['alpha_i'],r['growth'],r['absolute_polynomial_residual']))
    notes += ['',
        '实际比较：400→1597的复alpha差 %.9g，相对复alpha差 %.6g%%，增长率差 %.6g%%；799→1597的复alpha差 %.9g，增长率差 %.6g%%。两个裁域（nmax=.4410353及1.2049509，相比原域7.6991726）的复alpha差最大 %.3g，映射向量MAC约1。799点CubicSpline相对同点数PCHIP的复alpha差 %.9g，增长率差 %.6g%%。'%(allfine['delta_alpha'],100*allfine['relative_complex_alpha_difference'],100*allfine['relative_growth_difference'],fine['delta_alpha'],100*fine['relative_growth_difference'],crop,method['delta_alpha'],100*method['relative_growth_difference']),
        '原400点与各派生离散的物理n梯形积分primitive向量MAC最小 %.12f；两最细PCHIP之间MAC=%.12f。向量插值使用物理n上的复CubicSpline且不外推，积分对五个原始变量等权，不是物理能量。'%(min(r['normal_integral_primitive_MAC'] for r in report['comparisons']),fine['normal_integral_primitive_MAC']),
        '所有6个短调用各有2条有效微扰beta采样行、0拒绝，PETSc错误及未使用选项为0，实际初值附近只求nev1/ncv8；没有重新搜谱。最终报告alpha取精确omega=0、beta=−7恢复求解的raw PEP pair，其native索引逐一对齐。裁域法向边界明显移近且删除外层陡变区，而本候选近壁根不变，是此参数点对这些裁域不敏感的数值证据；不能据此判断所有外层根均无物理意义。','',
        '加密生成方法：原计算坐标eta_j=j−1，j=1..400；对grid及全部六个baseflow分量沿这个统一eta轴插值，eta范围保持0..399，没有外推。799点为间隔1/2、1597点为间隔1/4，保留全部原节点。原节点最大坐标差≤1.43e−14、baseflow差≤2.23e−16。虚部维持精确零，全部值有限，密度和温度保持正值；各H5形状/属性中的jn已更新。CubicSpline使用not-a-knot端条件，裁域直接取原j的前340/360点而不插值。','',
        '资料及局限：运行环境SciPy '+scipy.__version__+'、NumPy '+np.__version__+'。PCHIP保持数据单调性而不在陡变处过冲，但只有一阶导数保证连续、二阶导数可在节点跳变；CubicSpline为二阶连续。已核对与实际版本一致的官方资料：[PCHIP 1.10.1](https://docs.scipy.org/doc/scipy-1.10.1/reference/generated/scipy.interpolate.PchipInterpolator.html)、[CubicSpline 1.10.1](https://docs.scipy.org/doc/scipy-1.10.1/reference/generated/scipy.interpolate.CubicSpline.html)。因此本次将插值方法差异与离散差异都列出，不将“新增插值节点”写成独立CFD加密解。',
        '代表点的非零增长保持正号、向量归属高度一致，但中性根附近可能比本点敏感，后续实际包络关键中性点仍应做同类法向检查。没有据此宣布整个流动或整个谱网格无关。','',
        '实际计算量：6个求解器短调用累计墙钟（并发任务各自耗时求和）%.3f秒，最大单调用%.3f秒，最多3并发；派生H5合计%d字节。整个首次命令约9.7秒。原H5与二进制运行前后SHA-256核对一致。'%(sum(r['seconds'] for r in rows),max(r['seconds'] for r in rows),sum(x['bytes'] for v in report['variants'] for x in v['files'])),
        '原始文件与派生文件SHA-256如下。生成方法、版本、原节点误差等完整存于inputs/*/generation.json及outputs/normal_validation.json。','',
        '|文件|SHA-256|','|---|---|']
    for name,h in report['source_hashes'].items():notes.append('|原M1/in/'+name+'|`'+h+'`|')
    for v in report['variants']:
        for f in v['files']:notes.append('|inputs/'+v['name']+'/'+Path(f['path']).name+'|`'+f['sha256']+'`|')
    notes += ['',
        '复现：在WSL运行本目录validate_normal.py；派生输入及已完成结果有精确hash缓存检查，不覆盖部分失败产物。各runs/*保存实际config、命令、全部stdout、neutral field/run log及native vectors。主输出副本为RUN/outputs/local_mode_audit_normal_validation.json。',
        '结束审查：六个实际案例及七组向量/特征值比较已完成，必要原始证据和生成脚本/H5全部保留。没有失败求解，没有修改或删除原文件，没有需要删除的本任务尝试性代码。派生H5约160MiB，属于复现输入而保留。未完成事项：独立CFD网格解、其他站/中性点的法向收敛及物理模态身份。']
    text='\n\n'.join(notes)+'\n'
    (OUT/'验证说明.md').write_text(text,encoding='utf-8')
    (RUN/'outputs/local_mode_audit_normal_validation_说明.md').write_text(text,encoding='utf-8')
    with (OUT/'cases.csv').open('w',newline='') as f:
        keys=['name','jn','nmax_N','alpha_r','alpha_i','growth','absolute_polynomial_residual','seconds','passed']
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows({k:r[k] for k in keys} for r in rows)
    log('最终审查：6case/12有效field、7比较、10派生H5 hash与冻结binary/source核对通过；中文方法/误差/限制和CSV已保存。没有失败调用或需删除的尝试性代码，必要复现H5及原生向量保留。')

def main():
    assert driver.sha(BIN)==HASH
    if not LOG.exists():
        LOG.write_text('# M3near 法向分辨率与外边界敏感性工作日志\n\n日期2026-10-22，Asia/Shanghai。\n\n目标：i119、omega=0、beta=−7的已验近壁根进行独立法向离散和外边界敏感性检查。前提：只读原H5，冻结v5不改源。jn从config读取，H5/PETSc随jn加载，法向离散由派生坐标度量重组；源码已定点查阅cfgio_adapter.f90、petsc_loader.f90、lst_eqn.f90。\n\n策略：原400节点基准，PCHIP在原计算j坐标细分为799/1597节点并保持原节点值；原j前340/360节点裁域；799节点CubicSpline作为插值方法敏感性补充。只用已知根附近单target/nev1，不重新搜谱。比较实际alpha、growth、残差及映射到原物理n坐标的向量归属。\n\n完成标准：6个短调用的原始结果/失败、派生H5生成方法及hash、两级细化与两裁域定量误差、资料与限制、最终审查。预计输入生成几十秒、求解几秒至数十秒，3并发上限。\n\n重要限制：插值加密不产生新的CFD信息，最多检验同一输入剖面重构下的LST离散敏感性；不能证明原CFD解已网格独立。PCHIP只有C1，其二阶导数可在结点不连续，故另以CubicSpline比较。\n',encoding='utf-8')
    hashes={name:driver.sha(SOURCE/name) for name in ('grid_vertical2d_N.h5','flow_vertical2d_N.h5')}
    log('原始SHA核验记录 '+json.dumps(hashes)+'；准备独立派生数据（已有manifest与hash完全匹配时只读复用）。')
    variants=[]
    for name,method,jn in [('crop340','crop',340),('crop360','crop',360),('pchip799','PCHIP',799),('pchip1597','PCHIP',1597),('cubic799','CubicSpline',799)]:
        variants.append(generate_variant(name,method,jn,hashes))
    report=dict(created=stamp(),scope='i119 omega=0 beta=-7 near-wall seed only; fixed original CFD data',solver=str(BIN),solver_sha256=HASH,
                source_hashes=hashes,variants=variants,results=[],comparisons=[])
    specs=[('original400',SOURCE,400,'original')]+[(v['name'],HERE/'inputs'/v['name'],v['jn'],v['method']) for v in variants]
    with cf.ThreadPoolExecutor(max_workers=3) as ex:
        jobs={ex.submit(solve,*spec):spec for spec in specs}
        for job in cf.as_completed(jobs):
            try:r=job.result()
            except Exception as e:r=dict(name=jobs[job][0],passed=False,error=repr(e))
            report['results'].append(r);(OUT/'normal_validation.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
            log('实际结果就绪（精确缓存复用或本次求解，原started时间保留） '+json.dumps({k:r.get(k) for k in ('name','passed','seconds','error','alpha_r','alpha_i','delta_alpha_to_original','absolute_polynomial_residual')},ensure_ascii=False))
    reference=next(r for r in report['results'] if r['name']=='original400')
    if reference['passed']:
        for other in report['results']:
            if other['passed'] and other['name']!='original400':report['comparisons'].append(compare(reference,other))
        fine={r['name']:r for r in report['results']}
        if fine['pchip799']['passed'] and fine['pchip1597']['passed']:report['comparisons'].append(compare(fine['pchip799'],fine['pchip1597']))
        if fine['pchip799']['passed'] and fine['cubic799']['passed']:report['comparisons'].append(compare(fine['pchip799'],fine['cubic799']))
    report['source_unchanged']=all(driver.sha(SOURCE/name)==h for name,h in hashes.items());report['solver_unchanged']=driver.sha(BIN)==HASH
    report['completed']=stamp();(OUT/'normal_validation.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    (RUN/'outputs/local_mode_audit_normal_validation.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    if report['comparisons']:
        with (OUT/'comparisons.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(report['comparisons'][0]));w.writeheader();w.writerows(report['comparisons'])
    log('阶段审查：实际结果和全部失败已保存，源H5未变='+str(report['source_unchanged'])+'；冻结二进制未变='+str(report['solver_unchanged'])+'。尚须据定量误差写结论，不能以solver容差替代网格敏感性。')
    finalize(report)
    print(json.dumps(dict(results=[{k:r.get(k) for k in ('name','passed','error','alpha_r','alpha_i')} for r in report['results']],comparisons=report['comparisons']),indent=2),flush=True)

if __name__=='__main__':main()
