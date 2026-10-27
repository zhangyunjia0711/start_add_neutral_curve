"""Original bordered equations and exact Hessian for frequency/beta folds.

At fixed x a simple spatial eigenmode is neutral when the two real bordered
equations vanish. A merger of beta neutral roots satisfies det(J_alpha,beta)=0.
This avoids dividing by the eigenvalue derivative. Physical growth-sign tests
remain mandatory, including the condition that alpha is a simple eigenvalue.
"""
import numpy as np
from scipy.sparse import bmat,csc_matrix
from scipy.sparse.linalg import splu,norm
from solver import mac
W=6.0727867021443425/50
S=np.array([20.,10.,100.])
class Equations:
    def __init__(self,p,qref,pref):self.p,self.qref,self.pref=p,qref,pref
    def evaluate(self,z,hessian=True):
        a,b,f=z*S;p=self.p;C=p.baseC;Q=sum(a**k*b**j*C[k][j] for k in range(3) for j in range(3))+W*f*p.G
        A=sum(b**j*(C[1][j]+2*a*C[2][j]) for j in range(3))*S[0]
        B=sum(a**k*(C[k][1]+2*b*C[k][2]) for k in range(3))*S[1]
        F=W*p.G*S[2];D=[A,B,F];n=len(self.qref)
        M=bmat([[Q,csc_matrix(self.pref[:,None])],[csc_matrix(self.qref.conj()[None,:]),csc_matrix((1,1),dtype=complex)]],format='csc');lu=splu(M);rhs=np.r_[np.zeros(n,dtype=complex),1.+0j];v=lu.solve(rhs);w=lu.solve(rhs,trans='H')
        g=np.array([-np.vdot(w[:-1],d@v[:-1]) for d in D]);q=v[:-1]/np.linalg.norm(v[:-1]);pvec=w[:-1]/np.linalg.norm(w[:-1]);e=v[-1];H=np.zeros((3,3),dtype=complex)
        if hessian:
            AA=sum(2*b**j*C[2][j] for j in range(3))*S[0]**2
            AB=sum(2*b*(C[1][2]+2*a*C[2][2])+(C[1][1]+2*a*C[2][1]) for unused in range(1))*S[0]*S[1]
            BB=sum(2*a**k*C[k][2] for k in range(3))*S[1]**2
            second={(0,0):AA,(0,1):AB,(1,0):AB,(1,1):BB}
            for j in range(3):
                vt=lu.solve(np.r_[-D[j]@v[:-1],0j]);wt=lu.solve(np.r_[-D[j].conj().T@w[:-1],0j],trans='H')
                for i in range(3):
                    H[i,j]=-np.vdot(wt[:-1],D[i]@v[:-1])-np.vdot(w[:-1],D[i]@vt[:-1])
                    if (i,j) in second:H[i,j]-=np.vdot(w[:-1],second[i,j]@v[:-1])
        determinant=float(np.imag(np.conj(g[0])*g[1]));third=np.imag(np.conj(H[0])*g[1]+np.conj(g[0])*H[1])
        value=np.array([e.real,e.imag,determinant]);J=np.vstack([g.real,g.imag,third]);res=float(np.linalg.norm(Q@q));coeff=res/max(1.,sum(abs(a)**k*norm(K) for k,K in enumerate([row[0]+b*row[1]+b*b*row[2] for row in C])))
        return dict(F=value,J=J,q=q,p=pvec,right_residual=res,scaled_residual=coeff,gradient=g,hessian=H,MAC=mac(self.qref,q))
    def scale(self,z):
        r=self.evaluate(z);g=r['gradient'];return np.r_[np.maximum(np.linalg.norm(np.array([g.real,g.imag]),axis=1),1e-8),max(abs(g[0])*abs(g[1]),1e-8)]
    def derivative_check(self,z):
        r=self.evaluate(z);errors=[]
        for j in range(3):
            h=2e-6;delta=np.eye(3)[j]*h;plus=self.evaluate(z+delta);minus=self.evaluate(z-delta);fd=(plus['F']-minus['F'])/(2*h);error=np.linalg.norm(fd-r['J'][:,j])/max(1.,np.linalg.norm(fd));errors.append(float(error))
        assert max(errors)<2e-4,errors
        return errors
