"""Six fresh original-v10 spatial PEP calls at complex frequencies near saddle."""
import configparser,json,re,subprocess,time,traceback
from pathlib import Path
import numpy as np
import saddle_probe as sa
tp=sa.tp;RUN=sa.RUN;HERE=sa.HERE
OUT=RUN/'outputs/complex_alpha_saddle_native10_neighbors_20261024.json'

def native_pair(directory,source,alpha_target,omega,beta,qref,model):
    directory.mkdir(parents=True,exist_ok=True);saved=directory/'result.json'
    assert not saved.exists() and not (directory/'solver.log').exists()
    cfg=configparser.ConfigParser();cfg.read(source['config_path']);cfg['neutral']['enabled']='F'
    cfg['lst']['alpha']='%.17g, %.17g'%(alpha_target.real+1e-6,alpha_target.imag)
    cfg['lst']['omega']='%.17g, %.17g'%(omega.real,omega.imag);cfg['lst']['beta']='%.17g, 0'%beta
    with (directory/'config.ini').open('w') as f:cfg.write(f)
    cmd=[str(tp.BIN),'-f','config.ini','-pep_tol','1e-13','-pep_max_it','1000','-pep_nev','1','-pep_ncv','8',
        '-st_pc_factor_mat_solver_type','superlu','-options_left','-pep_view_vectors','binary:out/eigenvectors.petsc:default:append']
    stdin='-1 0 0\n';(directory/'stdin.txt').write_text(stdin);t=time.time()
    with (directory/'solver.log').open('w') as f:
        try:p=subprocess.run(cmd,cwd=directory,env=tp.driver.ENV,input=stdin,text=True,stdout=f,stderr=subprocess.STDOUT,timeout=30);code=p.returncode
        except subprocess.TimeoutExpired:code=124
    text=(directory/'solver.log').read_text(errors='replace');pairs=[]
    for m in re.finditer(r'^\s*(\d+)\s+\(\s*('+tp.vf.FLOAT+r')\s*,\s*('+tp.vf.FLOAT+r')\s*\)\s+The error is\s+('+tp.vf.FLOAT+r')',text,re.M):
        pairs.append(dict(index=int(m[1]),alpha_r=float(m[2].replace('D','E')),alpha_i=float(m[3].replace('D','E')),residual=float(m[4].replace('D','E'))))
    group_count=text.count('ncov=');vv=tp.vf.vectors(directory/'out/eigenvectors.petsc',400)
    assert group_count==1 and len(vv)==len(pairs),(group_count,len(vv),len(pairs))
    usable=[k for k,v in enumerate(pairs) if v['residual']<1e-7]
    assert usable,'no genuinely residual-qualified spatial root; preserve failed Ritz values'
    k=min(usable,key=lambda k:abs(complex(pairs[k]['alpha_r'],pairs[k]['alpha_i'])-alpha_target))
    alpha=complex(pairs[k]['alpha_r'],pairs[k]['alpha_i']);q=vv[k].ravel();q/=np.linalg.norm(q)
    coeff=model.at_beta(beta);poly=coeff[0]+omega*model.G+alpha*coeff[1]+alpha*alpha*coeff[2]
    rawmac=tp.mac(qref,q);matrixres=float(np.linalg.norm(poly@q))
    result=dict(created=tp.stamp(),directory=str(directory),command=cmd,stdin=stdin,seconds=time.time()-t,exit=code,
        omega=tp.cn(omega),beta=beta,alpha_target=tp.cn(alpha_target),actual_alpha=tp.cn(alpha),pairs=pairs,
        raw_vector_file=str(directory/'out/eigenvectors.petsc'),raw_vector_sha256=tp.driver.sha(directory/'out/eigenvectors.petsc'),
        raw_index=k,source_MAC=rawmac,reported_polynomial_residual=pairs[k]['residual'],exported_same_PEP_matrix_residual=matrixres,
        actual_pep_solves=group_count,petsc_errors=text.count('PETSC ERROR'),unused_options=text.count('Option left:'),
        passed=code==0 and group_count==1 and rawmac>=.95 and matrixres<1e-7 and pairs[k]['residual']<1e-7 and text.count('PETSC ERROR')==0 and text.count('Option left:')==0)
    saved.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');assert result['passed'],result
    return result,q

def main():
    assert not OUT.exists()
    path=RUN/'outputs/complex_alpha_saddle_native10_20261024.json';src=json.loads(path.read_text());assert src['status']=='complex_alpha_saddle_candidate' and src['local_second_derivative_verification_passed']
    old=json.loads(Path(src['matrix_provenance_path']).read_text());p=old['source_temporal_refinement']['point']
    model=tp.Pencil(RUN/'work/pep_matrix_temporal_probe_20261023/M1c_native/i10_matrices',old['source']['actual_native'],10,p['beta'],p['omega']['real']);model.build()
    point=src['result']['point'];a0=sa.z(point,'alpha');w0=sa.z(point,'omega');beta=src['fixed_beta'];q0=np.load(src['saddle_vector_file'])
    waa=sa.z(src['saddle_h_h2_verification'][1],'omega_alphalpha_real_direction')
    report=dict(created=tp.stamp(),source_path=str(path),source_sha256=tp.driver.sha(path),saddle=point,omega_alphalpha=tp.cn(waa),records=[],
        complex_frequency_interface_source=str(RUN/'outputs/complex_omega_interface_v10_20261024/complex_omega_probe.json'),
        pinch_origin_verified=False,physical_absolute_instability_claim=False,
        qualification='Two local spatial branches near a finite saddle, not their high-Imomega half-plane origin, not a proven pinch')
    sa.log('恢复原PEP邻支验证：root反馈WSL可能重启；此前鞍点已完成保存，无在途求解可声称持续。核对来源与bin后从文件断点开始6个原生interactive短调用。')
    tp.driver.snapshot();assert tp.driver.sha(tp.BIN)==tp.SHA
    OUT.write_text(json.dumps(report,indent=2)+'\n')
    for eps in [1e-4,1e-5,1e-6]:
        omega=w0+1j*eps;da=np.sqrt(2j*eps/waa);rr=[];qq=[]
        for sign in [-1,1]:
            r,q=native_pair(HERE/'runs'/('neighbor_eps%.0e_s%+d'%(eps,sign)),src['source_actual'],a0+sign*da,omega,beta,q0,model);rr.append(r);qq.append(q)
        roots=[sa.z(v,'actual_alpha') for v in rr];assert abs(roots[0]-roots[1])>1e-5,'two shifts returned the same root'
        record=dict(epsilon=eps,omega=tp.cn(omega),pairs=rr,root_separation=abs(roots[1]-roots[0]),
            separation_over_sqrt_epsilon=abs(roots[1]-roots[0])/np.sqrt(eps),midpoint=tp.cn((roots[0]+roots[1])/2),
            midpoint_difference_to_saddle=abs((roots[0]+roots[1])/2-a0),mutual_raw_MAC=tp.mac(qq[0],qq[1]),
            individual_square_law_relative_errors=[abs(.5*waa*(a-a0)**2-1j*eps)/eps for a in roots])
        report['records'].append(record);OUT.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
        print('NATIVE NEIGHBORS '+json.dumps({k:v for k,v in record.items() if k!='pairs'}),flush=True)
    tp.driver.snapshot();assert tp.driver.sha(tp.BIN)==tp.SHA and tp.driver.sha(path)==report['source_sha256']
    report.update(status='six_fresh_complex_frequency_PEP_roots_verified',completed=tp.stamp(),original_files_unchanged=True,
        total_native_calls=6,total_native_seconds=sum(v['seconds'] for r in report['records'] for v in r['pairs']))
    OUT.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');sa.log('原PEP6根通过，残差/MAC/矩阵回代/平方根收拢实测已保存；仍无两半平面起源，pinch未证。')

if __name__=='__main__':main()
