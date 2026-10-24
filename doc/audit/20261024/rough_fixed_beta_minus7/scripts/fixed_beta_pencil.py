"""Exact fixed-beta frequency pencil reconstructed from frozen spatial PEP matrices.

P(alpha,omega)=A0+alpha*A1+alpha^2*A2+omega*G. With real alpha,
K q = omega M q, M=-G; a zero Im(omega) is a neutral solution of
the same original spatial polynomial, without a Gaster approximation.
"""
import configparser,hashlib,json,re,struct,time,subprocess
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import LinearOperator,eigs,splu,norm,ArpackNoConvergence
from common import ROOT,SOLVER,SOLVER_HASH,OPTIONS,run,sha,driver
from profile_interpolation import mac,Pair

def matrix(path):
    data=path.read_bytes();cls,m,n,nz=struct.unpack_from('>iiii',data);assert cls==1211216
    offset=16;counts=np.frombuffer(data,dtype='>i4',count=m,offset=offset).astype(np.int32);offset+=4*m
    columns=np.frombuffer(data,dtype='>i4',count=nz,offset=offset).astype(np.int32);offset+=4*nz
    values=np.frombuffer(data,dtype='>c16',count=nz,offset=offset).astype(complex);offset+=16*nz
    assert offset==len(data) and counts.sum()==nz
    return csr_matrix((values,columns,np.r_[0,np.cumsum(counts)]),shape=(m,n)).tocsc()

def export(mode,station,state,omega,label):
    # Assemble through one ordinary spatial PEP call, then leave its interactive
    # loop. Ritz convergence is irrelevant to coefficient export; never reuse a
    # failed Ritz value as the next field-survey target.
    descriptor=dict(mode=str(mode),station=station,omega=float(omega),beta=-7,
        alpha_target=[state['alpha_r']+1e-4,state['alpha_i']],solver_sha256=SOLVER_HASH)
    token=hashlib.sha256(json.dumps(descriptor,sort_keys=True).encode()).hexdigest()[:12]
    directory=ROOT/'matrix_exports'/(label+'_'+token);directory.mkdir(parents=True,exist_ok=True)
    report=directory/'export.json'
    if report.exists():
        result=json.loads(report.read_text());assert result['descriptor']==descriptor
        assert all(sha(Path(p['path']))==p['sha256'] for p in result['files'])
        return [matrix(Path(p['path'])) for p in result['files']],result
    assert not (directory/'solver.log').exists(),'Partial matrix export retained'
    template=ROOT/str(mode);cfg=configparser.ConfigParser();cfg.read(template/'config.ini')
    cfg['filenames']['dir']=str(template/'in')+'/'
    cfg['domain']['iloc']=str(station);cfg['lst']['profile_coordinates']='wall_local'
    cfg['lst']['omega']='%.17g, 0'%omega;cfg['lst']['beta']='-7, 0'
    cfg['lst']['alpha']='%.17g, %.17g'%(state['alpha_r']+1e-4,state['alpha_i'])
    cfg['neutral']={'enabled':'F'}
    with (directory/'config.ini').open('w') as f:cfg.write(f)
    (directory/'out').mkdir(exist_ok=True)
    command=[str(SOLVER),'-f','config.ini','-pep_tol','1e-12','-pep_max_it','1','-pep_nev','1','-pep_ncv','8',
        '-st_pc_factor_mat_solver_type','superlu']
    for k in range(3):command+=['-pep_view_mat'+str(k),'binary:out/P'+str(k)+'.petsc']
    with (directory/'solver.log').open('w') as f:
        result=subprocess.run(command,cwd=directory,env=driver.ENV,input='-1\n',universal_newlines=True,stdout=f,stderr=subprocess.STDOUT,timeout=30)
    text=(directory/'solver.log').read_text();paths=[directory/'out'/('P'+str(k)+'.petsc') for k in range(3)]
    assert all(p.exists() for p in paths) and 'PETSC ERROR' not in text
    arrays=[matrix(p) for p in paths];assert all(np.all(np.isfinite(a.data)) for a in arrays)
    result=dict(descriptor=descriptor,command=command,matrix_only_export=True,
        files=[dict(path=str(p),sha256=sha(p)) for p in paths])
    report.write_text(json.dumps(result,indent=2)+'\n');return arrays,result

class Pencil:
    def __init__(self,mode,station,state,label,validate=False):
        self.mode,self.station,self.label=mode,station,label
        base=state['omega_r'];a,r0=export(mode,station,state,base,label+'_wa');b,r1=export(mode,station,state,base+1,label+'_wb')
        self.G=b[0]-a[0];self.M=-self.G;self.A=[a[0]-base*self.G,a[1],a[2]];self.exports=[r0,r1]
        assert norm(b[1]-a[1])<1e-10 and norm(b[2]-a[2])<1e-10
        self.validation=None
        if validate:
            c,rc=export(mode,station,state,base+.37,label+'_independent');self.exports.append(rc)
            errors=[float(norm(c[k]-(self.A[k]+(base+.37)*self.G if k==0 else self.A[k])))/max(1,float(norm(c[k]))) for k in range(3)]
            assert max(errors)<1e-12;self.validation=errors

    def evaluate(self,alpha,reference_q,omega_shift,left=True):
        start=time.monotonic();K=self.A[0]+alpha*self.A[1]+alpha*alpha*self.A[2]
        Ka=self.A[1]+2*alpha*self.A[2];shift=complex(omega_shift)+.001+.0001j
        lu=splu((K-shift*self.M).tocsc())
        op=LinearOperator(K.shape,matvec=lambda v:lu.solve(self.M@v),dtype=complex)
        chosen=None
        for k,ncv in [(1,12),(4,20)]:
            try:mu,q=eigs(op,k=k,which='LM',v0=reference_q,ncv=ncv,tol=1e-12,maxiter=300)
            except ArpackNoConvergence as error:mu,q=error.eigenvalues,error.eigenvectors
            if q is None or len(mu)==0:continue
            j=max(range(len(mu)),key=lambda j:mac(reference_q,q[:,j]));overlap=mac(reference_q,q[:,j])
            if overlap>=.95:chosen=j;break
        assert chosen is not None,('Temporal eigenvector continuity unresolved',overlap)
        v=q[:,j];v/=np.linalg.norm(v);omega=shift+1/mu[j]
        residual=float(np.linalg.norm((K-omega*self.M)@v));assert residual<1e-7
        result=dict(alpha_r=float(np.real(alpha)),alpha_i=float(np.imag(alpha)),omega_r=float(omega.real),omega_i=float(omega.imag),beta_r=-7.,
            raw_MAC=overlap,right_residual=residual,nev=k,seconds=time.monotonic()-start)
        if left:
            hop=LinearOperator(K.shape,matvec=lambda v:self.M.conj().T@lu.solve(v,trans='H'),dtype=complex)
            adjoint_ok=False
            for kl,nl in [(1,12),(4,20),(8,28)]:
                try:mm,w=eigs(hop,k=kl,which='LM',ncv=nl,tol=1e-12,maxiter=400)
                except ArpackNoConvergence as error:mm,w=error.eigenvalues,error.eigenvectors
                if w is None or len(mm)==0:continue
                jleft=int(np.argmin(abs(mm-np.conj(mu[j]))));p=lu.solve(w[:,jleft],trans='H');p/=np.linalg.norm(p)
                lr=float(np.linalg.norm((K-omega*self.M).conj().T@p));den=np.vdot(p,self.M@v)
                if lr<1e-7 and abs(den)>1e-14:adjoint_ok=True;break
            assert adjoint_ok,('No matching validated adjoint eigenpair',lr,abs(den))
            wa=np.vdot(p,Ka@v)/den
            result.update(left_residual=lr,adjoint_nev=kl,omega_alpha_r=float(wa.real),omega_alpha_i=float(wa.imag),
                original_real_frequency_polynomial_residual=float(np.linalg.norm((K-omega.real*self.M)@v)))
            if hasattr(self,'dA'):
                Kx=self.dA[0]+alpha*self.dA[1]+alpha*alpha*self.dA[2]
                wx=np.vdot(p,(Kx-omega*self.dM)@v)/den
                result.update(omega_xi_r=float(wx.real),omega_xi_i=float(wx.imag))
        return result,v

class OperatorLinearFamily:
    """C0 interpolation of exact native discrete pencils, with analytic gradients.

    Native points are unchanged; between-node points explicitly use the
    declared discrete-operator interpolation and are not extra CFD stations.
    """
    def __init__(self):
        from profile_interpolation import AXIS
        from functools import lru_cache
        self.axis=AXIS
        self.native=lru_cache(maxsize=16)(self._native)
    def _native(self,i):
        return Pencil('M3',i,self.axis[i],f'M3b7_native_matrix_i{i:03d}')
    def pencil(self,xi):
        assert 1<=xi<=361
        i=min(int(np.floor(xi)),360);theta=xi-i
        a,b=self.native(i),self.native(i+1)
        p=object.__new__(Pencil)
        p.A=[(1-theta)*aa+theta*bb for aa,bb in zip(a.A,b.A)]
        p.M=(1-theta)*a.M+theta*b.M;p.G=-p.M
        p.dA=[bb-aa for aa,bb in zip(a.A,b.A)];p.dM=b.M-a.M
        p.xi=xi;p.native_bracket=[i,i+1];p.theta=theta
        return p

class ProfileFamily:
    def __init__(self):self.pairs={};self.cache={}
    def pencil(self,xi,state,validate=False):
        assert 1<=xi<=361
        rounded=round(xi)
        if abs(xi-rounded)<1e-10:
            key=('native',rounded);mode='M3';station=rounded;x=None
        else:
            i=int(np.floor(xi));theta=xi-i
            if i not in self.pairs:self.pairs[i]=Pair(i,i+1)
            directory=self.pairs[i].build(theta);key=('interpolated',i,float(theta).hex());mode=directory;station=5
        if key not in self.cache:
            label='M3b7_temporal_'+hashlib.sha256(str(key).encode()).hexdigest()[:12]
            self.cache[key]=Pencil(mode,station,state,label,validate)
        return self.cache[key]
