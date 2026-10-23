"""No-core exact spatial-PEP coefficient export and finite temporal probe."""
import configparser,datetime,json,re,struct,subprocess,sys,time,traceback
from pathlib import Path
import numpy as np
import scipy
from scipy import sparse
from scipy.sparse.linalg import LinearOperator,eigs,splu

HERE=Path(__file__).resolve().parent
RUN=HERE.parent.parent
sys.path.insert(0,str(RUN/'work'))
sys.path.insert(0,str(RUN/'work/flat_neutral_normal_validation_20261023'))
import full_driver as driver
import validate_flat as vf
BIN=RUN/'outputs/bin/START_LST_wall_local_v10'
SHA='634b45c5d30e076f5b86ccb136e0734c1905beedde108477dceb8a75adc3853c'
OUT=RUN/'outputs/pep_matrix_temporal_probe_20261023.json'
LOG=RUN/'outputs/工作日志_2026-10-23_PEP矩阵重建temporal小试.md'

def stamp():return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
def log(t):
    with LOG.open('a') as f:f.write('\n'+stamp()+' '+t+'\n')
def save(j):OUT.write_text(json.dumps(j,indent=2,allow_nan=False)+'\n')
def mac(a,b):return float(abs(np.vdot(a,b))**2/(np.vdot(a,a).real*np.vdot(b,b).real))
def cn(z):return dict(real=float(z.real),imag=float(z.imag))
def norm(m):return float(sparse.linalg.norm(m))

def matrix(path):
    b=path.read_bytes();cls,m,n,nz=struct.unpack_from('>4i',b,0)
    assert cls==1211216 and m==n==2000 and nz>=0,(cls,m,n,nz)
    assert len(b)==16+4*m+4*nz+16*nz,(len(b),m,n,nz)
    lens=np.frombuffer(b,dtype='>i4',count=m,offset=16).astype(np.int64)
    assert int(lens.sum())==nz
    col=np.frombuffer(b,dtype='>i4',count=nz,offset=16+4*m).astype(np.int32)
    val=np.frombuffer(b,dtype='>c16',count=nz,offset=16+4*m+4*nz).astype(np.complex128)
    assert np.all((col>=0)&(col<n)) and np.all(np.isfinite(val))
    ptr=np.r_[0,np.cumsum(lens)]
    a=sparse.csr_matrix((val,col,ptr),shape=(m,n));a.sum_duplicates();a.sort_indices()
    return a

def export(label,omega,source):
    d=HERE/'runs'/label;d.mkdir(parents=True,exist_ok=True)
    resultfile=d/'export.json'
    if resultfile.exists():
        r=json.loads(resultfile.read_text());assert r['omega']==omega
        return r
    assert not (d/'solver.log').exists(),'preserve partial export, use new label'
    cfg=configparser.ConfigParser();cfg.read(source['config_path'])
    cfg['lst']['omega']='%.17g, 0'%omega
    cfg['neutral']['outer_min']=str(omega);cfg['neutral']['outer_max']=str(omega)
    with (d/'config.ini').open('w') as f:cfg.write(f)
    cmd=[str(BIN),'-f','config.ini','-pep_tol','1e-13','-pep_max_it','1000','-pep_nev','1','-pep_ncv','8',
        '-st_pc_factor_mat_solver_type','superlu','-pep_view_vectors','binary:out/all_vectors.petsc',
        '-pep_view_values','ascii:out/all_values.txt','-options_left']
    for k in range(3):cmd+=['-pep_view_mat%d'%k,'binary:out/P%d.petsc'%k]
    t=time.time()
    with (d/'solver.log').open('w') as f:
        f.write('COMMAND '+json.dumps(cmd)+'\n');f.flush()
        try:p=subprocess.run(cmd,cwd=d,env=driver.ENV,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,timeout=45);code=p.returncode
        except subprocess.TimeoutExpired:code=124
    text=(d/'solver.log').read_text(errors='replace')
    pairs=[];blocks=[]
    for m in re.finditer(r'^\s*(\d+)\s+\(\s*('+vf.FLOAT+r')\s*,\s*('+vf.FLOAT+r')\s*\)\s+The error is\s+('+vf.FLOAT+r')',text,re.M):
        if int(m[1])==0 and pairs:blocks.append(pairs);pairs=[]
        pairs.append(dict(index=int(m[1]),alpha_r=float(m[2].replace('D','E')),alpha_i=float(m[3].replace('D','E')),residual=float(m[4].replace('D','E'))))
    if pairs:blocks.append(pairs)
    files=[d/('out/P%d.petsc'%k) for k in range(3)]
    r=dict(label=label,omega=omega,beta=source['fixed_parameters']['beta_r'],started=stamp(),seconds=time.time()-t,
        directory=str(d),config=str(d/'config.ini'),command=cmd,exit=code,petsc_errors=text.count('PETSC ERROR'),
        unused_options=text.count('Option left:'),final_pairs=blocks[-1] if blocks else [],
        coefficient_files=[dict(path=str(p),sha256=driver.sha(p),bytes=p.stat().st_size) for p in files if p.exists()],
        native_vectors=str(d/'out/all_vectors.petsc'))
    r['export_passed']=len(r['coefficient_files'])==3 and r['unused_options']==0 and r['petsc_errors']==0
    resultfile.write_text(json.dumps(r,indent=2)+'\n')
    print('MATRIX EXPORT '+json.dumps({k:r[k] for k in ['label','seconds','exit','unused_options','export_passed']}),flush=True)
    return r

def main():
    assert not OUT.exists(),'this bounded probe is already recorded; inspect before retry'
    LOG.write_text('# PEP矩阵重建temporal小试\n\n目标：不修改core或原H5，从同一wall_local spatial PEP直接导出P0/P1/P2，核验P0对omega线性、P1/P2不依赖omega和原空间向量回代，再在M1c/i11已验中性点用稀疏shift-invert重建有限temporal频率根。完成标准：矩阵/原始向量/残差绑定、频率与MAC定量对照及失败界限；本轮不作全谱、中性优化或新鼻部物理结论。采用单线程，最多3原生导出+2个shift试算。原H5/bin/core保持。\n\n资料：实际查阅[SLEPc PEPSolve source](https://slepc.upv.es/release/src/pep/interface/pepsolve.c.html)确认-pep_view_matk；[SciPy1.10.1 eigs](https://docs.scipy.org/doc/scipy-1.10.1/reference/generated/scipy.sparse.linalg.eigs.html)确认LinearOperator接口与M厄米/正定限制。本试验将一般复数G留在自建算子中，不向eigs传M。installed库能力以实际导出结果为准，不把最新文档版本当本机版本。日期Asia/Shanghai。\n')
    driver.snapshot();assert driver.sha(BIN)==SHA
    source_path=RUN/'outputs/nose_neutral_native_vectors_20261023.json'
    source=next(v for v in json.loads(source_path.read_text())['cases'] if v['family']=='M1c')
    native=source['actual_native'];point=source['actual_point']
    assert native['passed'] and driver.sha(Path(native['native_vector_file']))==native['vector_sha256']
    report=dict(created=stamp(),status='running',solver=str(BIN),solver_sha256=SHA,source_path=str(source_path),
        source_sha256=driver.sha(source_path),source=source,versions=dict(numpy=np.__version__,scipy=scipy.__version__),
        exports=[],formulation='P(alpha,omega)=A+omega*G+alpha*P1+alpha**2*P2; OP=(K+sigma*G)^-1*(-G); omega=sigma+1/mu',
        references=['https://slepc.upv.es/release/src/pep/interface/pepsolve.c.html','https://docs.scipy.org/doc/scipy-1.10.1/reference/generated/scipy.sparse.linalg.eigs.html'])
    save(report)
    try:
        for label,omega in [('native_omega',point['omega_r']),('omega0',0.0),('omega1',1.0)]:
            r=export(label,omega,native);report['exports'].append(r);save(report)
            assert r['export_passed'],'matrix view or solver failed; preserve actual results'
        matrices=[[matrix(Path(f['path'])) for f in r['coefficient_files']] for r in report['exports']]
        ps,p0,p1=matrices
        A=p0[0];G=p1[0]-p0[0];omega=point['omega_r'];ar=point['alpha_r'];ac=complex(ar,point['alpha_i'])
        linear=ps[0]-(A+omega*G)
        checks=dict(shape=list(A.shape),nnz=[m.nnz for m in ps],P0_linearity_absolute_frobenius=norm(linear),
            P0_linearity_relative_frobenius=norm(linear)/max(norm(ps[0]),1),
            P1_changes=[norm(ms[1]-ps[1]) for ms in [p0,p1]],P2_changes=[norm(ms[2]-ps[2]) for ms in [p0,p1]],
            G_norm=norm(G),G_nnz=G.nnz)
        q0=vf.vectors(Path(native['native_vector_file']),400)[native['raw_index']].ravel();q0/=np.linalg.norm(q0)
        pairs=report['exports'][0]['final_pairs'];iq=min(range(len(pairs)),key=lambda k:abs(complex(pairs[k]['alpha_r'],pairs[k]['alpha_i'])-ac))
        aq=complex(pairs[iq]['alpha_r'],pairs[iq]['alpha_i']);qs=vf.vectors(Path(report['exports'][0]['native_vectors']),400)[iq].ravel();qs/=np.linalg.norm(qs)
        original_polynomial=ps[0]+aq*ps[1]+aq*aq*ps[2]
        checks.update(exported_source_alpha=cn(aq),exported_source_alpha_difference=abs(aq-ac),
            source_raw_MAC=mac(q0,qs),exported_source_reported_residual=pairs[iq]['residual'],
            exported_source_matrix_residual=float(np.linalg.norm(original_polynomial@qs)),
            prior_source_q_same_matrix_residual=float(np.linalg.norm((ps[0]+ac*ps[1]+ac*ac*ps[2])@q0)))
        checks['passed']=bool(checks['P0_linearity_relative_frobenius']<1e-12 and max(checks['P1_changes']+checks['P2_changes'])<1e-10 and checks['exported_source_matrix_residual']<1e-7 and checks['source_raw_MAC']>1-1e-10)
        report['matrix_checks']=checks;save(report);log('实际矩阵核验 '+json.dumps(checks))
        assert checks['passed'],'coefficient/source-vector parity failed; do not solve replacement'
        K=(A+ar*ps[1]+ar*ar*ps[2]).tocsc();report['temporal_probes']=[]
        for ds in [1e-3,1e-4]:
            sigma=omega+ds;t=time.time();lu=splu((K+sigma*G).tocsc())
            op=LinearOperator(K.shape,matvec=lambda v:lu.solve(-(G@v)),dtype=np.complex128)
            mu,qq=eigs(op,k=2,which='LM',v0=qs,ncv=12,tol=1e-12,maxiter=1000)
            values=sigma+1/mu;rows=[]
            for j,w in enumerate(values):
                q=qq[:,j];q/=np.linalg.norm(q)
                residual=float(np.linalg.norm((K+w*G)@q))
                direct_residual=float(np.linalg.norm((ps[0]+ar*ps[1]+ar*ar*ps[2]+(w-omega)*G)@q))
                np.save(HERE/('temporal_sigma%.0e_pair%d.npy'%(ds,j)),q)
                rows.append(dict(omega=cn(w),frequency_difference_to_spatial_source=abs(w-omega),
                    absolute_temporal_pencil_residual=residual,absolute_original_PEP_residual_at_real_alpha_complex_omega=direct_residual,
                    raw_MAC_to_source=mac(qs,q),vector_file=str(HERE/('temporal_sigma%.0e_pair%d.npy'%(ds,j)))))
            chosen=max(range(len(rows)),key=lambda k:rows[k]['raw_MAC_to_source'])
            r=dict(sigma=sigma,seconds=time.time()-t,results=rows,associated_index=chosen)
            r['passed']=rows[chosen]['frequency_difference_to_spatial_source']<1e-8 and abs(rows[chosen]['omega']['imag'])<1e-8 and rows[chosen]['raw_MAC_to_source']>1-1e-10 and rows[chosen]['absolute_original_PEP_residual_at_real_alpha_complex_omega']<1e-7
            report['temporal_probes'].append(r);save(report)
            log('实际shift-invert '+json.dumps(r));print('TEMPORAL PROBE '+json.dumps(r),flush=True)
        assert all(r['passed'] for r in report['temporal_probes']),'temporal association failed; no optimization authorized from this probe'
        report['status']='exact_PEP_temporal_reconstruction_local_probe_passed'
    except Exception:
        report['status']='bounded_probe_failed';report['exception']=traceback.format_exc();log(report['exception']);raise
    finally:
        driver.snapshot();assert driver.sha(BIN)==SHA and driver.sha(source_path)==report['source_sha256']
        report.update(completed=stamp(),source_and_original_files_unchanged=True,total_native_exports=len(report['exports']),
            native_export_seconds=sum(r['seconds'] for r in report['exports']),qualification='Finite local matrix reconstruction at one already verified native point only; no temporal core support claim, all-spectrum proof, optimized-beta neutral continuation or physical absolute-instability inference')
        save(report);log('结束审查：'+report['status']+'；原始文件/bin/source结果hash保持，必要配置/三系数矩阵/原生向量/脚本均为复现证据，保留；没有已确认不再需要的文件。')

if __name__=='__main__':main()
