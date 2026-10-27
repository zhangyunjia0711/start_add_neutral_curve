"""Frozen START polynomial at real omega=0; all beta coefficients exported.

Bordered Newton corrects the original 2000 x 2000 polynomial and its physical
left null vector. An explicit companion shift inverse is the fallback; neither
method inverts the singular leading matrix. Between-station operator matrices
are guides only; exported native matrices are the numerical evidence.
"""
import configparser, hashlib, json, re, struct, subprocess, time
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix, save_npz, load_npz, bmat, csc_matrix
from scipy.sparse.linalg import LinearOperator, eigs, splu, norm, ArpackNoConvergence
from common import ROOT, CASES, SOLVER, SOLVER_HASH, driver, sha, run, table

def mac(a,b):
    return float(abs(np.vdot(a,b))**2/(np.vdot(a,a).real*np.vdot(b,b).real))

def matrix(path):
    data=path.read_bytes();cls,m,n,nz=struct.unpack_from('>iiii',data);assert cls==1211216
    offset=16;counts=np.frombuffer(data,dtype='>i4',count=m,offset=offset).astype(np.int32);offset+=4*m
    cols=np.frombuffer(data,dtype='>i4',count=nz,offset=offset).astype(np.int32);offset+=4*nz
    values=np.frombuffer(data,dtype='>c16',count=nz,offset=offset).astype(complex);offset+=16*nz
    assert offset==len(data) and counts.sum()==nz
    return csr_matrix((values,cols,np.r_[0,np.cumsum(counts)]),shape=(m,n)).tocsc()

def vectors(path):
    data=path.read_bytes();offset=0;out=[]
    while offset<len(data):
        cls,n=struct.unpack_from('>ii',data,offset);offset+=8;assert cls==1211214 and n==2000
        out.append(np.frombuffer(data,dtype='>c16',count=n,offset=offset).astype(complex));offset+=16*n
    assert offset==len(data);return out

def export(case,i,beta,label=None):
    desc=dict(case=str(case),station=int(i),beta=float(beta),omega=0.,solver_sha256=SOLVER_HASH)
    token=hashlib.sha256(json.dumps(desc,sort_keys=True).encode()).hexdigest()[:12]
    directory=ROOT/'matrix_exports'/((label or str(case).split('/')[-1])+f'_i{i:03d}_'+token)
    report=directory/'export.json'
    if report.exists():
        result=json.loads(report.read_text());assert result['descriptor']==desc
        return [matrix(Path(p['path'])) for p in result['files']],result
    directory.mkdir(parents=True,exist_ok=True);assert not (directory/'solver.log').exists()
    template=ROOT/'cases'/str(case);cfg=configparser.ConfigParser();cfg.read(template/'config.ini')
    cfg['filenames']['dir']=str(template/'in')+'/'
    cfg['domain']['iloc']=str(i);cfg['lst']['profile_coordinates']='wall_local'
    cfg['lst']['omega']='0, 0';cfg['lst']['beta']='%.17g, 0'%beta
    cfg['lst']['alpha']='35.93418168, -4.02320790';cfg['neutral']={'enabled':'F'}
    with (directory/'config.ini').open('w') as f:cfg.write(f)
    (directory/'out').mkdir(exist_ok=True)
    command=[str(SOLVER),'-f','config.ini','-pep_tol','1e-12','-pep_max_it','1','-pep_nev','1','-pep_ncv','8','-st_pc_factor_mat_solver_type','superlu']
    for k in range(3):command+=['-pep_view_mat'+str(k),'binary:out/P'+str(k)+'.petsc']
    with (directory/'solver.log').open('w') as f:
        r=subprocess.run(command,cwd=directory,env=driver.ENV,input='-1\n',text=True,stdout=f,stderr=subprocess.STDOUT,timeout=40)
    paths=[directory/'out'/f'P{k}.petsc' for k in range(3)]
    assert all(p.exists() for p in paths) and 'PETSC ERROR' not in (directory/'solver.log').read_text()
    arrays=[matrix(p) for p in paths];assert all(np.all(np.isfinite(a.data)) for a in arrays)
    result=dict(descriptor=desc,command=command,matrix_only_export=True,exit=r.returncode,files=[dict(path=str(p),sha256=sha(p)) for p in paths])
    report.write_text(json.dumps(result,indent=2)+'\n');return arrays,result

def original(case,i,state,beta,label):
    h=1e-6
    spec=dict(station=i,mode=str(case),seed={k:float(state[k]) for k in ['omega_r','beta_r','alpha_r','alpha_i']},
        seed_shift=1e-5,solver=str(SOLVER),solver_sha256=SOLVER_HASH,coordinates='wall_local',task=33,
        outer_parameter='iloc',outer_min=i,outer_max=i,outer_points=1,
        inner_parameter='beta',inner_min=beta-h,inner_max=beta+h,inner_points=3,
        sample_only=True,anchor_row_seed=True,step=.05,jump=4,adaptive=True,
        vector_tracking=True,minimum_overlap=.95,minimum_step=1e-7,retry_eigenpairs=8,
        options=['-pep_tol','1e-13','-pep_max_it','1000','-pep_nev','1','-pep_ncv','8',
        '-st_pc_factor_mat_solver_type','superlu','-pep_view_vectors','binary:out/all_vectors.petsc:default:append'],timeout=180)
    token=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()[:12];spec['name']=label+'_'+token
    s=run(spec);directory=ROOT/'runs'/spec['name']
    assert s['exit']!=124 and not any(s[k] for k in ['petsc_errors','unused_options','failed_field_records','coordinate_metadata_missing'])
    values=sorted(table(directory/'out/neutral_field.plt'),key=lambda q:q['beta_r']);assert len(values)==3
    center=values[1];target=complex(center['alpha_r'],center['alpha_i'])
    floatpat=r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?'
    pairs=[(complex(float(m[1].replace('D','E')),float(m[2].replace('D','E'))),float(m[3].replace('D','E'))) for m in re.finditer(r'^\s*\d+\s+\(\s*('+floatpat+r')\s*,\s*('+floatpat+r')\s*\)\s+The error is\s+('+floatpat+r')',(directory/'solver.log').read_text(),re.M)]
    qs=vectors(directory/'out/all_vectors.petsc');assert len(qs)==len(pairs)
    k=min(range(len(pairs)),key=lambda k:abs(pairs[k][0]-target));assert abs(pairs[k][0]-target)<1e-8 and pairs[k][1]<1e-7
    return center,qs[k]/np.linalg.norm(qs[k]),dict(run=spec['name'],original_residual=pairs[k][1])

def qep_arnoldi(K,target,qref,left=False):
    start=time.monotonic();n=K[0].shape[0];sigma=complex(target)+.003+.001j
    S=(K[0]+sigma*K[1]+sigma*sigma*K[2]).tocsc();lu=splu(S);B=K[1]+sigma*K[2]
    def action(y):
        u,v=y[:n],y[n:];a=-lu.solve(B@u+K[2]@v);return np.r_[a,u+sigma*a]
    op=LinearOperator((2*n,2*n),matvec=action,dtype=complex)
    initial=np.r_[qref,target*qref];initial/=np.linalg.norm(initial)
    chosen=None;overlap=0.
    for k,ncv in [(1,12),(4,20),(8,28)]:
        try:mu,z=eigs(op,k=k,which='LM',v0=initial,ncv=ncv,tol=1e-12,maxiter=500)
        except ArpackNoConvergence as e:mu,z=e.eigenvalues,e.eigenvectors
        if z is None or len(mu)==0:continue
        j=max(range(len(mu)),key=lambda j:mac(qref,z[:n,j]));overlap=mac(qref,z[:n,j])
        if overlap>=.95:chosen=j;break
    assert chosen is not None,('Unresolved eigenvector continuity',overlap,target)
    lam=sigma+1/mu[j];q=z[:n,j];q/=np.linalg.norm(q)
    Q=K[0]+lam*K[1]+lam*lam*K[2];res=float(np.linalg.norm(Q@q));assert res<1e-7,('QEP residual',res,lam)
    r=dict(value_r=float(lam.real),value_i=float(lam.imag),raw_MAC=overlap,right_residual=res,nev=k,seconds=time.monotonic()-start)
    p=None
    if left:
        # Recenter after mode selection: the tracked mode need not have been
        # among the nearest few adjoint Ritz values at the old predictor shift.
        sigma=lam+.003+.001j;lu=splu((K[0]+sigma*K[1]+sigma*sigma*K[2]).tocsc());B=K[1]+sigma*K[2]
        reference_mu=1/(lam-sigma)
        def c_inverse_h(y):
            p2=-lu.solve(y[:n]+np.conj(sigma)*y[n:],trans='H');p1=y[n:]+B.conj().T@p2
            return p1,p2
        def adj_action(y):
            p1,p2=c_inverse_h(y);return np.r_[p1,K[2].conj().T@p2]
        hop=LinearOperator((2*n,2*n),matvec=adj_action,dtype=complex);ok=False
        for kl,nl in [(1,12),(4,20),(8,28)]:
            try:mm,w=eigs(hop,k=kl,which='LM',ncv=nl,tol=1e-12,maxiter=500)
            except ArpackNoConvergence as e:mm,w=e.eigenvalues,e.eigenvectors
            if w is None or len(mm)==0:continue
            jl=int(np.argmin(abs(mm-np.conj(reference_mu))));_,p=c_inverse_h(w[:,jl]);p/=np.linalg.norm(p)
            lr=float(np.linalg.norm(Q.conj().T@p));den=np.vdot(p,(K[1]+2*lam*K[2])@q)
            if lr<1e-7 and abs(den)>1e-14:ok=True;break
        assert ok,('Adjoint unresolved',lr,abs(den))
        r.update(left_residual=lr,adjoint_nev=kl)
    return lam,q,r,p

def qep(K,target,qref,left=False):
    """Bordered Newton correction on the original polynomial; Arnoldi fallback.

    A bordered physical null-vector solve avoids crowded constraint eigenvalues
    in the companion adjoint spectrum; close-shift inverse iteration is a backup.
    """
    start=time.monotonic();lam=complex(target);q=qref/np.linalg.norm(qref);ref=q.copy();ok=False
    for it in range(18):
        Q=(K[0]+lam*K[1]+lam*lam*K[2]).tocsc();res=float(np.linalg.norm(Q@q))
        if res<5e-10 and it>0:ok=True;break
        derivative=(K[1]+2*lam*K[2])@q
        bordered=bmat([[Q,csc_matrix(derivative[:,None])],[csc_matrix(ref.conj()[None,:]),csc_matrix((1,1),dtype=complex)]],format='csc')
        update=splu(bordered).solve(np.r_[-Q@q,0j]);dlam=update[-1]
        if abs(dlam)>max(5.,abs(lam)) or not np.isfinite(dlam):break
        q+=update[:-1];lam+=dlam;q/=np.linalg.norm(q)
    overlap=mac(ref,q)
    if not ok or overlap<.95:return qep_arnoldi(K,target,qref,left)
    Q=(K[0]+lam*K[1]+lam*lam*K[2]).tocsc();res=float(np.linalg.norm(Q@q));assert res<1e-7
    r=dict(value_r=float(lam.real),value_i=float(lam.imag),raw_MAC=overlap,right_residual=res,nev=0,newton_iterations=it,seconds=time.monotonic()-start,method='original_polynomial_bordered_Newton')
    p=None
    if left:
        derivative=(K[1]+2*lam*K[2])@q
        bordered_h=bmat([[Q.conj().T,csc_matrix(ref[:,None])],[csc_matrix(derivative.conj()[None,:]),csc_matrix((1,1),dtype=complex)]],format='csc')
        p=splu(bordered_h).solve(np.r_[np.zeros(len(q),dtype=complex),1j*0+1])[:-1];p/=np.linalg.norm(p)
        lr=float(np.linalg.norm(Q.conj().T@p));den=np.vdot(p,derivative)
        found=lr<1e-9 and abs(den)>1e-14;k=-1
        for shift in ([] if found else [1e-10+3e-11j,1e-12+3e-13j]):
            sigma=lam+shift;lu=splu((K[0]+sigma*K[1]+sigma*sigma*K[2]).tocsc());p=q.copy()
            for k in range(6):
                p=lu.solve(p,trans='H');p/=np.linalg.norm(p)
                lr=float(np.linalg.norm(Q.conj().T@p));den=np.vdot(p,(K[1]+2*lam*K[2])@q)
                if lr<1e-9 and abs(den)>1e-14:found=True;break
            if found:break
        if not found:return qep_arnoldi(K,lam,q,left)
        r.update(left_residual=lr,adjoint_nev=0,adjoint_inverse_iterations=k+1,adjoint_method='bordered_physical_null_vector')
    return lam,q,r,p

class Pencil:
    def __init__(self,case,i,validate=False):
        self.case,self.i=case,i;directory=ROOT/'native_pencils'/str(case)/f'i{i:03d}'
        directory.mkdir(parents=True,exist_ok=True);report=directory/'coefficients.json'
        if report.exists():
            self.C=[[load_npz(directory/f'C{k}{j}.npz') for j in range(3)] for k in range(3)]
            self.report=json.loads(report.read_text())
        else:
            H=16.;minus,rm=export(case,i,-H);zero,r0=export(case,i,0.);plus,rp=export(case,i,H)
            self.C=[[zero[k],(plus[k]-minus[k])/(2*H),(plus[k]+minus[k]-2*zero[k])/(2*H*H)] for k in range(3)]
            for k in range(3):
                for j in range(3):save_npz(directory/f'C{k}{j}.npz',self.C[k][j])
            self.report=dict(case=str(case),station=i,beta_nodes=[-H,0,H],exports=[rm,r0,rp],validations=[])
            report.write_text(json.dumps(self.report,indent=2)+'\n')
        if validate and not self.report['validations']:
            for beta in [-7.,.37,-35.]:
                actual,r=export(case,i,beta);pred=self.spatial(beta)
                errors=[float(norm(a-b))/max(1.,float(norm(a))) for a,b in zip(actual,pred)]
                assert max(errors)<2e-12,errors
                self.report['validations'].append(dict(beta=beta,relative_matrix_errors=errors,export=r))
            report.write_text(json.dumps(self.report,indent=2)+'\n')
    def spatial(self,beta):return [c[0]+beta*c[1]+beta*beta*c[2] for c in self.C]
    def transverse(self,alpha):return [self.C[0][j]+alpha*self.C[1][j]+alpha*alpha*self.C[2][j] for j in range(3)]
    def solve_alpha(self,beta,target,qref,left=False):
        lam,q,r,p=qep(self.spatial(beta),target,qref,left)
        r.update(alpha_r=float(lam.real),alpha_i=float(lam.imag),beta_r=float(beta),omega_r=0.,growth=float(-lam.imag),iloc=self.i)
        if left:
            Pb=sum(lam**k*(self.C[k][1]+2*beta*self.C[k][2]) for k in range(3));Pa=self.spatial(beta)[1]+2*lam*self.spatial(beta)[2]
            ab=-np.vdot(p,Pb@q)/np.vdot(p,Pa@q);r.update(alpha_beta_r=float(ab.real),alpha_beta_i=float(ab.imag))
        return lam,q,r

if __name__=='__main__':
    from common import snapshot
    snapshot();results=[]
    for case,c in CASES.items():
        i=c['seed_station'];state=c['seed'];p=Pencil(case,i,True)
        point,q,original_receipt=original(case,i,state,-7.,case+'_seed')
        a,q,r=p.solve_alpha(-7.,complex(point['alpha_r'],point['alpha_i']),q,True)
        np.savez(ROOT/'outputs'/f'{case}_seed.npz',q=q)
        record=dict(case=case,original=point,receipt=original_receipt,spatial=r,scan=[])
        for direction in [1,-1]:
            b=-7.;aa=a;qq=q;prev=None;record['scan'].append(dict(r,direction=direction))
            while (-45<b<-.049):
                step=min(.2,max(.02,abs(b)*.04));nb=b+direction*step
                if direction==1:nb=min(nb,-.05)
                while True:
                    target=aa if prev is None else aa+(nb-b)*(aa-prev[1])/(b-prev[0])
                    try:na,nq,rr=p.solve_alpha(nb,target,qq,False);break
                    except AssertionError as e:
                        nb=(nb+b)/2
                        if abs(nb-b)<1e-5:
                            record['scan'].append(dict(direction=direction,stopped_at_beta=b,reason=str(e)));break
                if abs(nb-b)<1e-5:break
                prev=(b,aa);aa,qq=na,nq;record['scan'].append(dict(rr,direction=direction));b=nb
                if b==-.05:break
        # Independently solve the positive-beta conjugate through the original PEP.
        positive=dict(omega_r=0.,beta_r=7.,alpha_r=-a.real,alpha_i=a.imag)
        pp,qp,pr=original(case,i,positive,7.,case+'_positive_conjugate')
        record['symmetry']=dict(alpha_difference=float(abs(complex(pp['alpha_r'],pp['alpha_i'])+np.conj(a))),raw_MAC=mac(q.conj(),qp),receipt=pr)
        record['matrix_symmetry']=[[float(norm(((-1)**(k+j))*p.C[k][j]-p.C[k][j].conj()))/max(1.,float(norm(p.C[k][j]))) for j in range(3)] for k in range(3)]
        results.append(record);(ROOT/'outputs/pilot.json').write_text(json.dumps(results,indent=2)+'\n')
        print(json.dumps(dict(case=case,seed=r,symmetry=record['symmetry'],scan=[{k:v[k] for k in ['beta_r','alpha_r','growth']} for v in record['scan'][::10]])),flush=True)
