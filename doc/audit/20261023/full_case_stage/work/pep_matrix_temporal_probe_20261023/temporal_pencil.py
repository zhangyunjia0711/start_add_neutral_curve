"""Private no-core spatial-PEP polynomial reconstruction and temporal derivatives.

All beta reconstructions require actual independent-matrix checks before use.
Artificial frequency eigenvalues are always validated against the original pencil.
"""
import configparser,datetime,json,re,subprocess,sys,time
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import LinearOperator,eigs,splu,norm
import probe

RUN=probe.RUN;driver=probe.driver;vf=probe.vf;BIN=probe.BIN;SHA=probe.SHA
stamp=probe.stamp;mac=probe.mac;cn=probe.cn;matrix=probe.matrix

def export(directory,source,station,beta,omega):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True);saved=directory/'export.json'
    params=dict(station=station,beta=beta,omega=omega)
    if saved.exists():
        r=json.loads(saved.read_text());assert r['parameters']==params and r['solver_sha256']==SHA
        assert all(driver.sha(Path(x['path']))==x['sha256'] for x in r['coefficient_files'])
        return r
    assert not (directory/'solver.log').exists(),'preserve partial export; choose new owned directory'
    assert driver.sha(BIN)==SHA
    cfg=configparser.ConfigParser();cfg.read(source['config_path'])
    cfg['domain']['iloc']=str(station)
    cfg['lst']['omega']='%.17g, 0'%omega;cfg['lst']['beta']='%.17g, 0'%beta
    cfg['neutral']['outer_parameter']='omega';cfg['neutral']['inner_parameter']='beta'
    for k in ['outer_min','outer_max']:cfg['neutral'][k]=str(omega)
    cfg['neutral']['inner_min']=str(beta-1e-8);cfg['neutral']['inner_max']=str(beta+1e-8)
    with (directory/'config.ini').open('w') as f:cfg.write(f)
    cmd=[str(BIN),'-f','config.ini','-pep_tol','1e-13','-pep_max_it','1000','-pep_nev','1','-pep_ncv','8',
        '-st_pc_factor_mat_solver_type','superlu','-pep_view_vectors','binary:out/all_vectors.petsc',
        '-pep_view_values','ascii:out/all_values.txt','-options_left']
    for k in range(3):cmd+=['-pep_view_mat%d'%k,'binary:out/P%d.petsc'%k]
    t=time.time()
    with (directory/'solver.log').open('w') as f:
        f.write('COMMAND '+json.dumps(cmd)+'\n');f.flush()
        try:p=subprocess.run(cmd,cwd=directory,env=driver.ENV,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,timeout=45);code=p.returncode
        except subprocess.TimeoutExpired:code=124
    text=(directory/'solver.log').read_text(errors='replace');blocks=[];pairs=[]
    for m in re.finditer(r'^\s*(\d+)\s+\(\s*('+vf.FLOAT+r')\s*,\s*('+vf.FLOAT+r')\s*\)\s+The error is\s+('+vf.FLOAT+r')',text,re.M):
        if int(m[1])==0 and pairs:blocks.append(pairs);pairs=[]
        pairs.append(dict(index=int(m[1]),alpha_r=float(m[2].replace('D','E')),alpha_i=float(m[3].replace('D','E')),residual=float(m[4].replace('D','E'))))
    if pairs:blocks.append(pairs)
    paths=[directory/('out/P%d.petsc'%k) for k in range(3)]
    r=dict(parameters=params,solver=str(BIN),solver_sha256=SHA,started=stamp(),seconds=time.time()-t,
        directory=str(directory),config=str(directory/'config.ini'),command=cmd,exit=code,
        petsc_errors=text.count('PETSC ERROR'),unused_options=text.count('Option left:'),
        final_pairs=blocks[-1] if blocks else [],native_vectors=str(directory/'out/all_vectors.petsc'),
        coefficient_files=[dict(path=str(p),sha256=driver.sha(p),bytes=p.stat().st_size) for p in paths if p.exists()])
    r['export_passed']=len(r['coefficient_files'])==3 and r['petsc_errors']==0 and r['unused_options']==0
    saved.write_text(json.dumps(r,indent=2)+'\n')
    print('TEMPORAL MATRIX '+json.dumps(dict(parameters=params,seconds=r['seconds'],passed=r['export_passed'])),flush=True)
    assert r['export_passed'],'matrix export failure retained'
    return r

class Pencil:
    def __init__(self,directory,source,station,beta0,known_omega):
        self.directory=Path(directory);self.source=source;self.station=station;self.beta0=beta0
        self.known_omega=known_omega;self.exports=[];self.calls=0
    def fetch(self,label,beta,omega):
        r=export(self.directory/label,self.source,self.station,float(beta),float(omega));self.exports.append(r)
        return [matrix(Path(v['path'])) for v in r['coefficient_files']]
    def build(self):
        c=self.fetch('b0_w0',self.beta0,0);lo=self.fetch('bm1_w0',self.beta0-1,0);hi=self.fetch('bp1_w0',self.beta0+1,0)
        self.coef=[(c[k],(hi[k]-lo[k])*.5,(hi[k]+lo[k]-2*c[k])*.5) for k in range(3)]
        w1=self.fetch('b0_w1',self.beta0,1);self.G=w1[0]-c[0];self.M=-self.G
        test=self.fetch('b037_w0',self.beta0+.37,0);test1=self.fetch('b037_w1',self.beta0+.37,1)
        actualG=test1[0]-test[0];Gdiff=actualG-self.G
        rebuilt=self.at_beta(self.beta0+.37)
        checks=dict(independent_beta=self.beta0+.37,
            independent_matrix_absolute_errors=[float(norm(rebuilt[k]-test[k])) for k in range(3)],
            independent_matrix_relative_errors=[float(norm(rebuilt[k]-test[k]))/max(float(norm(test[k])),1) for k in range(3)],
            G_beta_independence_absolute_error=float(norm(Gdiff)),
            G_beta_independence_relative_error=float(norm(Gdiff))/max(float(norm(self.G)),1),
            omega1_P1_P2_errors=[float(norm(test1[k]-test[k])) for k in [1,2]])
        checks['passed']=max(checks['independent_matrix_relative_errors'])<1e-12 and checks['G_beta_independence_absolute_error']<1e-10 and max(checks['omega1_P1_P2_errors'])<1e-10
        self.checks=checks
        (self.directory/'polynomial_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
        assert checks['passed'],'actual independent beta matrix parity failed'
        return checks
    def at_beta(self,beta,derivative=0):
        d=beta-self.beta0
        if derivative==0:return [a+d*b+d*d*c for a,b,c in self.coef]
        if derivative==1:return [b+2*d*c for a,b,c in self.coef]
        if derivative==2:return [2*c for a,b,c in self.coef]
        raise ValueError(derivative)
    def matrices(self,alpha,beta):
        a=self.at_beta(beta);b=self.at_beta(beta,1)
        K=(a[0]+alpha*a[1]+alpha*alpha*a[2]).tocsc()
        Ka=a[1]+2*alpha*a[2];Kb=b[0]+alpha*b[1]+alpha*alpha*b[2]
        return K,Ka,Kb
    def evaluate(self,alpha,beta,reference_q,omega_shift=None,left=True):
        self.calls+=1;t=time.time();K,Ka,Kb=self.matrices(alpha,beta)
        shift=complex(self.known_omega+1e-3 if omega_shift is None else omega_shift)
        S=(K-shift*self.M).tocsc();lu=splu(S)
        OP=LinearOperator(K.shape,matvec=lambda v:lu.solve(self.M@v),dtype=np.complex128)
        mu,qv=eigs(OP,k=2,which='LM',v0=reference_q,ncv=12,tol=1e-12,maxiter=1000)
        omegas=shift+1/mu;idx=max(range(len(mu)),key=lambda k:mac(reference_q,qv[:,k]));q=qv[:,idx];q/=np.linalg.norm(q)
        omega=omegas[idx];res=float(np.linalg.norm((K-omega*self.M)@q));overlap=mac(reference_q,q)
        r=dict(alpha=float(alpha),beta=float(beta),omega=cn(omega),raw_MAC=overlap,right_residual=res,shift=cn(shift),seconds=0)
        assert res<1e-7 and overlap>=.95,('temporal right eigenpair/branch failure',r)
        if left:
            HOP=LinearOperator(K.shape,matvec=lambda v:self.M.conj().T@lu.solve(v,trans='H'),dtype=np.complex128)
            mv,wv=eigs(HOP,k=2,which='LM',ncv=12,tol=1e-12,maxiter=1000)
            j=int(np.argmin(abs(mv-np.conj(mu[idx]))));w=wv[:,j]
            p=lu.solve(w,trans='H');p/=np.linalg.norm(p)
            leftres=float(np.linalg.norm((K-omega*self.M).conj().T@p))
            den=np.vdot(p,self.M@q);wa=np.vdot(p,Ka@q)/den;wb=np.vdot(p,Kb@q)/den
            assert leftres<1e-7 and abs(den)>1e-14,('temporal left eigenpair failure',leftres,abs(den))
            assert abs(wa)>1e-9,'omega_alpha too small to resolve implicit branch'
            cross=float(np.imag(wb*np.conj(wa)))
            r.update(left_residual=leftres,left_transformed_mu_error=float(abs(mv[j]-np.conj(mu[idx]))),
                pH_M_q=cn(den),omega_alpha=cn(wa),omega_beta=cn(wb),cross_numerator=cross,
                cross_divided_abs_omega_alpha_squared=cross/abs(wa)**2,
                implicit_alpha_beta=cn(-wb/wa),eigenvalue_condition_estimate=float(np.linalg.norm(p)*np.linalg.norm(q)/abs(den)))
        r['seconds']=time.time()-t
        return r,q,(p if left else None)

def objective(r):
    # Equivalent simple-branch stationarity, scaled so a small omega_alpha
    # cannot make the raw cross numerator look converged.
    wa=complex(r['omega_alpha']['real'],r['omega_alpha']['imag'])
    return np.array([r['omega']['imag']/abs(wa),r['cross_divided_abs_omega_alpha_squared']],float)

def correct(pencil,alpha,beta,q,omega_shift,maxiter=20,backtracks=4):
    trace=[];currentq=q;wshift=omega_shift;last=None
    for it in range(maxiter):
        r,v,p=pencil.evaluate(alpha,beta,currentq,wshift,left=True);F=objective(r);record=dict(iteration=it,center=r,F=F.tolist(),trials=[]);trace.append(record)
        if max(abs(F))<1e-7:
            return dict(status='temporal_candidate',point=r,trace=trace,iterations=it+1),v
        # Smooth temporal derivatives permit substantially wider steps than
        # the 1e-7 spatial FD windows; h/h2 diagnostics are retained separately.
        hs=[max(1e-4,abs(alpha)*1e-5),max(1e-5,abs(beta)*1e-5)]
        J=np.empty((2,2));stencil=[]
        for k,h in enumerate(hs):
            pair=[]
            for s in [-1,1]:
                args=[alpha,beta];args[k]+=s*h
                rr,_,_=pencil.evaluate(*args,v,complex(r['omega']['real'],r['omega']['imag'])+1e-3,left=True);pair.append(rr)
            J[:,k]=(objective(pair[1])-objective(pair[0]))/(2*h);stencil.append(pair)
        record.update(J=J.tolist(),J_condition=float(np.linalg.cond(J)),h=hs,stencil=stencil)
        delta=np.linalg.solve(J,-F);scale=min(1,.5/max(abs(delta[0]),1e-300),.2/max(abs(delta[1]),1e-300));delta*=scale
        record['delta']=delta.tolist();accepted=False
        for bt in range(backtracks):
            factor=.5**bt;a=alpha+factor*delta[0];b=beta+factor*delta[1]
            try:
                rr,vv,_=pencil.evaluate(a,b,v,complex(r['omega']['real'],r['omega']['imag'])+1e-3,left=True);ff=objective(rr)
                trial=dict(factor=factor,point=rr,F=ff.tolist());record['trials'].append(trial)
                if np.linalg.norm(ff)<np.linalg.norm(F):
                    alpha,beta,currentq,wshift=a,b,vv,complex(rr['omega']['real'],rr['omega']['imag'])+1e-3;accepted=True;last=rr;break
            except Exception as exc:record['trials'].append(dict(factor=factor,error=str(exc)))
        if not accepted:return dict(status='unresolved_no_descent',trace=trace,iterations=it+1,last_actual=last),currentq
    return dict(status='unresolved_max_iterations',trace=trace,iterations=maxiter,last_actual=last),currentq
