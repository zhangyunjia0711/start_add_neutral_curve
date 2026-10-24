"""Independent native spatial beta-curvature check of temporal stationary root."""
import configparser,json,subprocess,time
from pathlib import Path
import numpy as np
import temporal_pencil as tp
RUN=tp.RUN;HERE=Path(__file__).resolve().parent/'M1c_native'
OUT=RUN/'outputs/matrix_temporal_M1c_native10_curvature_20261023.json'
LOG=RUN/'outputs/工作日志_2026-10-23_矩阵频率换元M1c原生鼻部延拓.md'
def save(r):OUT.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
def main():
    assert not OUT.exists()
    sourcepath=RUN/'outputs/matrix_temporal_M1c_native10_secant_20261023.json';s=json.loads(sourcepath.read_text());r=s['result']['point'];native=s['native_fresh']
    first=json.loads((RUN/'outputs/matrix_temporal_M1c_native_20261023.json').read_text());b0=first['source_temporal_refinement']['point']['beta'];w0=first['source_temporal_refinement']['point']['omega']['real']
    model=tp.Pencil(HERE/'i10_matrices',first['source']['actual_native'],10,b0,w0);model.build()
    q=np.load(HERE/'i10_secant_temporal_candidate_vector.npy');wa=complex(r['omega_alpha']['real'],r['omega_alpha']['imag']);ab=r['implicit_alpha_beta']['real']
    report=dict(created=tp.stamp(),station=10,source_path=str(sourcepath),source_sha256=tp.driver.sha(sourcepath),point=r,
        chain_formula='g_betabeta=partial_beta[Im(omega_beta/omega_alpha)] + Re(alpha_beta)*partial_alpha[Im(omega_beta/omega_alpha)] at neutral simple spatial beta stationary point',
        qualification='Our implicit-function chain derivation, checked independently by actual native spatial beta samples',chain_checks=[],native_spatial_checks=[])
    save(report)
    for h in [1e-3,5e-4,1e-4,5e-5]:
        vals=[]
        for da,db in [(-h,0),(h,0),(0,-h),(0,h)]:
            v,_,_=model.evaluate(r['alpha']+da,r['beta']+db,q,r['omega']['real']+1e-3,left=True);vals.append(v)
        f=[v['cross_divided_abs_omega_alpha_squared'] for v in vals]
        Ja=(f[1]-f[0])/(2*h);Jb=(f[3]-f[2])/(2*h);curv=Jb+ab*Ja
        report['chain_checks'].append(dict(h=h,partial_alpha=Ja,partial_beta=Jb,g_betabeta=curv,samples=vals));save(report)
    for h in [1e-6,5e-7,2.5e-7]:
        label='M1c_i10_beta_curvature_h%.1e'%h;d=HERE/'runs'/label;d.mkdir(parents=True)
        cfg=configparser.ConfigParser();cfg.read(native['config_path'])
        cfg['neutral']['inner_min']=str(r['beta']-h);cfg['neutral']['inner_max']=str(r['beta']+h)
        cfg['neutral']['scan_points']='3'
        with (d/'config.ini').open('w') as f:cfg.write(f)
        cmd=list(native['command']);t=time.time()
        with (d/'solver.log').open('w') as f:
            try:p=subprocess.run(cmd,cwd=d,env=tp.driver.ENV,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,timeout=45);code=p.returncode
            except subprocess.TimeoutExpired:code=124
        field=tp.driver.table(d/'out/neutral_field.plt');text=(d/'solver.log').read_text(errors='replace')
        g=[v['growth'] for v in field];runlog=(d/'out/neutral_run.log').read_text();overlaps=[float(x.replace('D','E')) for x in __import__('re').findall(r'vector_overlap=\s*('+tp.vf.FLOAT+r')',runlog)]
        check=dict(h=h,seconds=time.time()-t,command=cmd,directory=str(d),exit=code,petsc_errors=text.count('PETSC ERROR'),
            field=field,minimum_MAC=min(overlaps) if overlaps else None)
        check['passed']=len(field)==3 and all(v['status']==0 for v in field) and code in (0,1) and check['petsc_errors']==0 and check['minimum_MAC']>=.95
        if check['passed']:check.update(g_beta=(g[2]-g[0])/(2*h),g_betabeta=(g[2]-2*g[1]+g[0])/h**2)
        report['native_spatial_checks'].append(check);save(report)
    native_pass=all(v['passed'] for v in report['native_spatial_checks'])
    negative=native_pass and all(v['g_betabeta']<0 for v in report['native_spatial_checks']) and all(v['g_betabeta']<0 for v in report['chain_checks'])
    positive=native_pass and all(v['g_betabeta']>0 for v in report['native_spatial_checks']) and all(v['g_betabeta']>0 for v in report['chain_checks'])
    report.update(status='native_beta_local_maximum_verified' if negative else ('native_beta_local_minimum_verified_not_envelope_peak' if positive else 'native_curvature_unresolved'),
        physical_beta_maximum_candidate=negative,completed=tp.stamp())
    tp.driver.snapshot();assert tp.driver.sha(tp.BIN)==tp.SHA and tp.driver.sha(sourcepath)==report['source_sha256'];report['original_files_unchanged']=True;save(report)
    with LOG.open('a') as f:f.write('\n'+tp.stamp()+' native10曲率审查 '+report['status']+'；chain='+str([(x['h'],x['g_betabeta']) for x in report['chain_checks']])+'；原生='+str([(x['h'],x.get('g_betabeta')) for x in report['native_spatial_checks']])+'。实际驻点若为beta最小值不可并入beta最大增长包络。\n')
    print('CURVATURE '+json.dumps(dict(status=report['status'],chain=[(x['h'],x['g_betabeta']) for x in report['chain_checks']],native=[(x['h'],x.get('g_betabeta')) for x in report['native_spatial_checks']])) ,flush=True)

if __name__=='__main__':main()
