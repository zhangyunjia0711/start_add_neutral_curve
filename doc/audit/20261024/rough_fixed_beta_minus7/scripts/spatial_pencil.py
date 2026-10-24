"""Exact shift-invert companion of the original quadratic spatial polynomial."""
import time
import numpy as np
from scipy.sparse.linalg import LinearOperator,eigs,splu,ArpackNoConvergence
from profile_interpolation import mac

def solve(pencil,omega,target,qref,minimum_mac=.95):
    start=time.monotonic();C0=pencil.A[0]+omega*pencil.G;C1,C2=pencil.A[1:]
    n=C0.shape[0];sigma=complex(target)+.003+.001j
    lu=splu((C0+sigma*C1+sigma*sigma*C2).tocsc())
    def action(y):
        u,v=y[:n],y[n:]
        z=-lu.solve((C1+sigma*C2)@u+C2@v)
        return np.r_[z,u+sigma*z]
    op=LinearOperator((2*n,2*n),matvec=action,dtype=complex)
    v0=np.r_[qref,complex(target)*qref];v0/=np.linalg.norm(v0)
    for k,ncv in [(1,12),(4,20),(8,28)]:
        try:mu,v=eigs(op,k=k,which='LM',v0=v0,ncv=ncv,tol=1e-12,maxiter=500)
        except ArpackNoConvergence as e:mu,v=e.eigenvalues,e.eigenvectors
        if v is None or len(mu)==0:continue
        j=max(range(len(mu)),key=lambda j:mac(qref,v[:n,j]));overlap=mac(qref,v[:n,j])
        if overlap<minimum_mac:continue
        alpha=sigma+1/mu[j];q=v[:n,j];q=q/np.linalg.norm(q)
        residual=float(np.linalg.norm((C0+alpha*C1+alpha*alpha*C2)@q))
        assert residual<1e-7,('Original spatial polynomial residual',residual)
        return dict(alpha_r=float(alpha.real),alpha_i=float(alpha.imag),omega_r=float(omega),omega_i=0.,beta_r=-7.,
            growth=-float(alpha.imag),raw_MAC=overlap,original_spatial_polynomial_residual=residual,nev=k,seconds=time.monotonic()-start),q
    raise RuntimeError('No converged spatial eigenvector passed MAC gate')
