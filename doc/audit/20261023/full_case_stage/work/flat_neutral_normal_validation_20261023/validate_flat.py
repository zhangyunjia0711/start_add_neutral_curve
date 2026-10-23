"""Fixed-parameter normal-discretization checks at two actual flat neutral starts."""
import configparser,csv,datetime,json,re,struct,subprocess,sys,time
from pathlib import Path
import h5py,numpy as np
from scipy.interpolate import CubicSpline

HERE=Path(__file__).resolve().parent
RUN=HERE.parent.parent
sys.path.insert(0,str(RUN/'work'))
import full_driver as driver
OLD=RUN/'work/local_mode_normal_validation_20261022/inputs'
SOURCE=RUN.parent/'M1/in'
BIN=RUN/'outputs/bin/START_LST_wall_local_v7'
SHA='7aa6016570f406587156f229107c0e5e4f76c8a85ba89dda267039ff1d688d79'
FLOAT=r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?'
LOG=RUN/'outputs/工作日志_2026-10-23_平板中性代表点法向核验.md'

def stamp():return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
def log(t):
    with LOG.open('a') as f:f.write('\n'+stamp()+' '+t+'\n')

def vectors(path,jn):
    b=path.read_bytes();off=0;result=[]
    while off<len(b):
        cls,n=struct.unpack_from('>ii',b,off);off+=8
        assert cls==1211214 and n==5*jn
        end=off+16*n;assert end<=len(b)
        result.append(np.frombuffer(b[off:end],dtype='>c16').astype(np.complex128).reshape(jn,5));off=end
    return result

def metadata(directory,station):
    with (directory/'wall_local_profile_metadata.csv').open() as f:
        rows=list(csv.DictReader(line for line in f if not line.startswith('#')))
    return next(r for r in rows if int(r['iloc'])==station)

def normal(input_dir,station,meta):
    with h5py.File(input_dir/'grid_vertical2d_N.h5','r') as f:g=f['grid'][:,station-1,:,0]
    n=(g-g[0])@np.array([float(meta['n_x']),float(meta['n_y'])]);n[0]=0
    assert np.all(np.diff(n)>0)
    return n

def solve(name,input_dir,jn,station,point,here=HERE,binary=BIN,sha=SHA,domain_in=None,pep_nev=1,pep_ncv=8,pep_tol='1e-13'):
    d=here/'runs'/name;d.mkdir(parents=True,exist_ok=True);saved=d/'result.json'
    if saved.exists():
        old=json.loads(saved.read_text());assert old['fixed_parameters']=={k:point[k] for k in ('omega_r','beta_r')}
        assert old['station']==station and old['solver_sha256']==sha
        return old
    assert not (d/'solver.log').exists(),'partial run retained; choose a new name'
    assert driver.sha(binary)==sha
    cfg=configparser.ConfigParser();cfg.read(RUN.parent/'M1/config.ini')
    cfg['filenames']['dir']=str(input_dir)+'/'
    cfg['domain']['jn']=str(jn);cfg['domain']['iloc']=str(station)
    if domain_in is not None:
        cfg['domain']['in']=str(domain_in)
        cfg['domain']['istart']='1';cfg['domain']['iend']=str(domain_in)
    cfg['lst']['profile_coordinates']='wall_local'
    cfg['lst']['alpha']='%.16g, %.16g'%(point['alpha_r']+1e-4,point['alpha_i'])
    cfg['lst']['omega']='%.16g, 0'%point['omega_r'];cfg['lst']['beta']='%.16g, 0'%point['beta_r']
    cfg['neutral']={k:str(v) for k,v in dict(enabled='T',track_initial_target='T',task=33,
        outer_parameter='omega',inner_parameter='beta',outer_min=point['omega_r'],outer_max=point['omega_r'],outer_points=1,
        inner_min=point['beta_r']-1e-8,inner_max=point['beta_r']+1e-8,scan_points=2,
        growth_tol=1e-8,parameter_tol=1e-7,max_iterations=70,eigen_error_tol=1e-7,max_eigen_jump=64,
        continuation_step=.005,output_prefix='out/neutral',sample_only='T',adaptive_steps='T',
        minimum_step=1e-6,vector_tracking='T',minimum_overlap=.95,retry_eigenpairs=8,anchor_row_seed='T').items()}
    with (d/'config.ini').open('w') as f:cfg.write(f)
    cmd=[str(binary),'-f','config.ini','-pep_tol',str(pep_tol),'-pep_max_it','1000','-pep_nev',str(pep_nev),'-pep_ncv',str(pep_ncv),
         '-st_pc_factor_mat_solver_type','superlu','-pep_view_vectors','binary:out/all_vectors.petsc',
         '-pep_view_values','ascii:out/all_values.txt','-options_left']
    t=time.time();started=stamp()
    with (d/'solver.log').open('w') as f:
        f.write('COMMAND '+json.dumps(cmd)+'\n');f.flush()
        try:p=subprocess.run(cmd,cwd=d,env=driver.ENV,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,timeout=180);code=p.returncode
        except subprocess.TimeoutExpired:code=124
    text=(d/'solver.log').read_text(errors='replace');blocks=[];block=[]
    for m in re.finditer(r'^\s*(\d+)\s+\(\s*('+FLOAT+r')\s*,\s*('+FLOAT+r')\s*\)\s+The error is\s+('+FLOAT+r')',text,re.M):
        index=int(m[1])
        if index==0 and block:blocks.append(block);block=[]
        block.append(dict(index=index,alpha_r=float(m[2].replace('D','E')),alpha_i=float(m[3].replace('D','E')),absolute_polynomial_residual=float(m[4].replace('D','E'))))
    if block:blocks.append(block)
    pairs=blocks[-1] if blocks else []
    field=driver.table(d/'out/neutral_field.plt')
    runlog=(d/'out/neutral_run.log').read_text() if (d/'out/neutral_run.log').exists() else ''
    m=re.search(r'evaluations=\s*(\d+)\s+solves=\s*(\d+)\s+rejected=\s*(\d+)\s+result=(.+)',runlog)
    overlaps=[float(x.replace('D','E')) for x in re.findall(r'vector_overlap=\s*('+FLOAT+r')',runlog)]
    result=dict(name=name,station=station,jn=jn,input_dir=str(input_dir),started=started,seconds=round(time.time()-t,3),
        solver=str(binary),solver_sha256=driver.sha(binary),command=cmd,config_path=str(d/'config.ini'),exit=code,
        fixed_parameters={k:point[k] for k in ('omega_r','beta_r')},petsc_errors=text.count('PETSC ERROR'),
        unused_options=text.count('Option left:'),terminal=m[4].strip() if m else 'MISSING',
        solves=int(m[2]) if m else None,rejected=int(m[3]) if m else None,
        field_records=len(field),failed_field_records=sum(r['status']!=0 for r in field),
        minimum_tiny_beta_MAC=min(overlaps) if overlaps else None,pairs=pairs,passed=False)
    path=d/'out/all_vectors.petsc'
    if pairs and path.exists():
        vv=vectors(path,jn);assert len(vv)==len(pairs)
        target=complex(point['alpha_r'],point['alpha_i']);k=min(range(len(pairs)),key=lambda k:abs(complex(pairs[k]['alpha_r'],pairs[k]['alpha_i'])-target))
        a=complex(pairs[k]['alpha_r'],pairs[k]['alpha_i']);meta=metadata(d,station);n=normal(input_dir,station,meta)
        q=vv[k];amp=np.sum(abs(q)**2,axis=1);peak=int(np.argmax(amp))
        result.update(alpha_r=a.real,alpha_i=a.imag,growth=-a.imag,delta_alpha_to_original_candidate=abs(a-target),
            absolute_polynomial_residual=pairs[k]['absolute_polynomial_residual'],raw_index=k,
            native_vector_file=str(path),vector_sha256=driver.sha(path),nmax_N=float(n[-1]),metadata=meta,
            peak_j=peak+1,peak_n_N=float(n[peak]),raw_vector_index_association_verified=True)
        result['passed']=bool(code in (0,1) and result['petsc_errors']==0 and result['unused_options']==0 and result['rejected']==0 and
            len(field)==2 and result['failed_field_records']==0 and result['absolute_polynomial_residual']<=1e-7 and
            result['delta_alpha_to_original_candidate']<64 and result['minimum_tiny_beta_MAC']>=.95)
    saved.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('FIXED NORMAL '+json.dumps({k:result.get(k) for k in ('name','seconds','passed','alpha_r','alpha_i','growth')}),flush=True)
    return result

def compare(base,r,slope):
    n0=normal(Path(base['input_dir']),base['station'],base['metadata']);v0=vectors(Path(base['native_vector_file']),base['jn'])[base['raw_index']]
    n1=normal(Path(r['input_dir']),r['station'],r['metadata']);v1=vectors(Path(r['native_vector_file']),r['jn'])[r['raw_index']]
    keep=n0<=n1[-1]+1e-12;n=n0[keep];a=v0[keep];b=CubicSpline(n1,v1,axis=0,extrapolate=False)(n)
    weights=np.empty(len(n));weights[1:-1]=(n[2:]-n[:-2])/2;weights[0]=(n[1]-n[0])/2;weights[-1]=(n[-1]-n[-2])/2;weights=weights[:,None]
    mac=float(abs(np.sum(weights*np.conj(a)*b))**2/(np.sum(weights*abs(a)**2)*np.sum(weights*abs(b)**2)))
    dg=r['growth']-base['growth'];da=complex(r['alpha_r']-base['alpha_r'],r['alpha_i']-base['alpha_i'])
    return dict(reference=base['name'],other=r['name'],delta_alpha=abs(da),growth_difference=dg,
        normal_integral_equal_primitive_MAC=mac,mapped_points=len(n),original_g_omega=slope,
        linear_frequency_shift_estimate=-dg/slope,neutral_point_reconverged=False,
        qualification='Fixed original beta, original400 g_omega linearization only; no re-converged neutral/optimized-beta result')

def main():
    assert not LOG.exists(),'use exact saved results to resume; do not duplicate an active log'
    LOG.write_text('# 平板两个代表中性点法向敏感性\n\n目标：i321 M1p高频中性点与i841 M2p上根，原固定omega/beta分别用400/799/1597实际求解。采用既有PCHIP派生H5，不修改原始文件；不能把i119结果外推本两点。先核对实际来源/版本/站号及hash，单线程各3短调用，比较alpha/g/实际向量映射MAC，以原gomega仅估计频率位移，最后审查。完成标准是6个真实solve和误差/局限，不是每网格重新中性收敛或独立CFD网格验证。日期Asia/Shanghai。\n')
    provenance=[json.loads((OLD/name/'generation.json').read_text()) for name in ['pchip799','pchip1597']]
    for m in provenance:
        assert all(driver.sha(Path(f['path']))==f['sha256'] for f in m['files'])
        assert all(driver.sha(SOURCE/f)==h for f,h in m['source_hashes'].items())
    driver.snapshot();assert driver.sha(BIN)==SHA
    report=dict(created=stamp(),status='running',solver=str(BIN),solver_sha256=SHA,reused_input_provenance=provenance,cases=[])
    target=RUN/'outputs/flat_neutral_normal_validation_20261023.json'
    for family,station in [('M1p',321),('M2p',841)]:
        path=RUN/f'outputs/neutral_x_cont_v7a_{family}_downstream_progress.json'
        source=next(v for v in json.loads(path.read_text())['records'] if v['station']==station and v['status']=='candidate_converged')
        assert source['result']['coordinates']=='wall_local' and source['result']['solver_sha256']==SHA
        point=source['point'];slope=source['result']['verification']['J'][0][0]
        case=dict(family=family,station=station,point=point,source_path=str(path),source_sha256=driver.sha(path),source_record=source,original_g_omega=slope,results=[],comparisons=[])
        for name,inp,jn in [('original400',SOURCE,400),('pchip799',OLD/'pchip799',799),('pchip1597',OLD/'pchip1597',1597)]:
            r=solve(f'{family}_i{station}_{name}',inp,jn,station,point)
            case['results'].append(r)
            assert r['passed'],'keep actual failure; do not force association'
            assert float(r['metadata']['t_x'])==1 and float(r['metadata']['t_y'])==0
            assert float(r['metadata']['n_x'])==0 and float(r['metadata']['n_y'])==1
            log(json.dumps({k:r.get(k) for k in ['name','passed','growth','alpha_r','alpha_i']},ensure_ascii=False))
        base,medium,fine=case['results']
        case['comparisons']=[compare(base,medium,slope),compare(base,fine,slope),compare(medium,fine,slope)]
        for c in case['comparisons']:
            if c['reference']==base['name']:c['estimated_omega_at_fixed_beta']=point['omega_r']+c['linear_frequency_shift_estimate']
        report['cases'].append(case);target.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    assert all(driver.sha(Path(f['path']))==f['sha256'] for m in provenance for f in m['files'])
    assert all(driver.sha(SOURCE/f)==h for m in provenance for f,h in m['source_hashes'].items())
    assert all(driver.sha(Path(c['source_path']))==c['source_sha256'] for c in report['cases'])
    driver.snapshot();assert driver.sha(BIN)==SHA
    report.update(status='fixed_parameter_sensitivity_complete',completed=stamp(),source_and_derived_hashes_unchanged=True,
        total_calls=6,total_seconds=sum(r['seconds'] for c in report['cases'] for r in c['results']),
        qualification='Reconstruction of fixed CFD data; no independent CFD refinement or new-grid neutral/beta optimum convergence')
    target.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    lines=['# 平板代表中性点法向敏感性',stamp(),
        '从实际wall_local/v7沿程结果取i321 M1p高频根及i841 M2p上根；未用global/v4旧fd_validation2混配。固定各自原omega/beta，复用既有PCHIP799/1597输入，原CFD数据不变。每站用实际solver metadata的法向投影n；本两点t=(1,0)、n=(0,1)，未沿用鼻部径向近似。',
        '|来源/站号|jn|alpha_r|g=−alpha_i|原PEP残差|','|---|---:|---:|---:|---:|']
    for c in report['cases']:
        for r in c['results']:lines.append('|%s/i%d|%d|%.12g|%.12g|%.4g|'%(c['family'],c['station'],r['jn'],r['alpha_r'],r['growth'],r['absolute_polynomial_residual']))
    lines+=['','|来源/比较|复alpha差|g差|物理n积分primitive MAC|线性Delta omega|','|---|---:|---:|---:|---:|']
    for c in report['cases']:
        for x in c['comparisons']:lines.append('|%s %s→%s|%.9g|%.9g|%.12g|%.9g|'%(c['family'],x['reference'].split('_')[-1],x['other'].split('_')[-1],x['delta_alpha'],x['growth_difference'],x['normal_integral_equal_primitive_MAC'],x['linear_frequency_shift_estimate']))
    lines+=['','频率估计仅用Delta omega=−Delta g/(g_omega)_400，固定原beta，未在新网格重收敛中性点，也未重新优化beta。g接近0时不给无意义的相对g百分比。映射MAC以实际物理n坐标、复CubicSpline无外推、五primitive变量等权梯形积分；不代表物理能量范数或模态身份验证。',
        '派生H5在原计算j坐标PCHIP细分，原节点保留，仅检查固定CFD重构下的LST离散敏感性，不能代替独立CFD网格收敛。[SciPy1.10.1 PCHIP](https://docs.scipy.org/doc/scipy-1.10.1/reference/generated/scipy.interpolate.PchipInterpolator.html)与[CubicSpline](https://docs.scipy.org/doc/scipy-1.10.1/reference/generated/scipy.interpolate.CubicSpline.html)资料于前任务实际核验，本轮复用同版本同方法。',
        '六次单线程调用合计%.3f秒；全部真实样本、原PEP残差、native向量索引/长度和MAC通过。原H5、派生H5、binary与来源trace hash均保持。保留必要复现代码、配置、原生向量、日志与结果；无不再需要的尝试文件。结论范围只覆盖这两个参数点。'%report['total_seconds']]
    md=[];table=[]
    for line in lines:
        if line.startswith('|'):table.append(line);continue
        if table:md.append('\n'.join(table));table=[]
        if line:md.append(line)
    if table:md.append('\n'.join(table))
    target.with_name(target.stem+'_说明.md').write_text('\n\n'.join(md)+'\n')
    log('结束审查：6次固定参数调用完成，全部hash保持；误差和线性估计见JSON/说明，无新网格中性收敛声明，无不必要文件可清理。')
    print('FLAT NORMAL RESULT '+str(target),flush=True)

if __name__=='__main__':main()
