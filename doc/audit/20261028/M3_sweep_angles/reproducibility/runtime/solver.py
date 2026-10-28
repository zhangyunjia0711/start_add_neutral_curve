"""Joint beta-quadratic, omega-linear polynomial exported from frozen START."""
import configparser,hashlib,json,subprocess,time
from pathlib import Path
from functools import lru_cache
import numpy as np
from scipy.sparse import save_npz,load_npz
from scipy.sparse.linalg import norm
from common import ROOT,CASES,SOLVER,SOLVER_HASH,driver,sha,snapshot
from core import matrix,qep,mac,original
PREVIOUS=ROOT

def export(case,i,beta,omega,label=None):
    assert str(case)!='s75', 'User requested read-only reuse of existing 75-degree results'
    desc=dict(case=str(case),station=int(i),beta=float(beta),omega=float(omega),solver_sha256=SOLVER_HASH)
    token=hashlib.sha256(json.dumps(desc,sort_keys=True).encode()).hexdigest()[:12]
    directory=ROOT/'matrix_exports'/((label or str(case).split('/')[-1])+f'_i{i:03d}_'+token);report=directory/'export.json'
    if report.exists():
        result=json.loads(report.read_text());assert result['descriptor']==desc
        return [matrix(Path(p['path'])) for p in result['files']],result
    directory.mkdir(parents=True,exist_ok=True);assert not (directory/'solver.log').exists()
    from slice_input import template as get_template
    template,local_station,slice_receipt=get_template(case,i)
    cfg=configparser.ConfigParser();cfg.read(template/'config.ini')
    cfg['filenames']['dir']=str(template/'in')+'/'
    cfg['domain']['iloc']=str(local_station);cfg['lst']['profile_coordinates']='wall_local_consistent'
    cfg['lst']['omega']='%.17g, 0'%omega;cfg['lst']['beta']='%.17g, 0'%beta
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
    result=dict(descriptor=desc,command=command,matrix_only_export=True,exit=r.returncode,local_station=local_station,native_slice=slice_receipt,files=[dict(path=str(p),sha256=sha(p)) for p in paths])
    report.write_text(json.dumps(result,indent=2)+'\n');return arrays,result

@lru_cache(maxsize=12)
def _coefficients_unlocked(case,i):
    directory=PREVIOUS/'native_pencils'/case/f'i{i:03d}'
    from core import Pencil as BetaPencil
    BetaPencil(case,i)
    prior=json.loads((directory/'coefficients.json').read_text())
    assert prior['case']==case and prior['station']==i
    C=[[load_npz(directory/f'C{k}{j}.npz') for j in range(3)] for k in range(3)]
    new=ROOT/'native_pencils'/case/f'i{i:03d}';new.mkdir(parents=True,exist_ok=True);report=new/'frequency_coefficient.json'
    if report.exists():
        receipt=json.loads(report.read_text());G=load_npz(new/'G.npz')
    else:
        a,source=export(case,i,0.,1.)
        differences=[a[k]-C[k][0] for k in range(3)]
        assert norm(differences[1])<1e-10 and norm(differences[2])<1e-10
        G=differences[0];save_npz(new/'G.npz',G)
        receipt=dict(case=case,station=i,prior_beta_coefficients=str(directory/'coefficients.json'),prior_report_sha256=sha(directory/'coefficients.json'),
            prior_coefficient_files=[dict(path=str(directory/f'C{k}{j}.npz'),sha256=sha(directory/f'C{k}{j}.npz')) for k in range(3) for j in range(3)],
            export=source,frequency_independent_higher_coefficients=[float(norm(differences[1])),float(norm(differences[2]))],validations=[])
        report.write_text(json.dumps(receipt,indent=2)+'\n')
    return C,G,receipt

def coefficients(case,i):
    import fcntl
    lockdir=ROOT/'work/cache_locks';lockdir.mkdir(exist_ok=True)
    token=hashlib.sha256((str(case)+':'+str(i)).encode()).hexdigest()
    with (lockdir/(token+'.lock')).open('a') as f:
        fcntl.flock(f,fcntl.LOCK_EX)
        return _coefficients_unlocked(case,i)


class Pencil:
    def __init__(self,case,i,omega=0.):
        assert str(case)!='s75', 'User requested read-only reuse of existing 75-degree results'
        self.case,self.i,self.omega=case,i,float(omega);self.baseC,self.G,self.receipt=coefficients(case,i)
        self.C=[list(row) for row in self.baseC];self.C[0][0]=self.C[0][0]+self.omega*self.G
    def spatial(self,beta):return [row[0]+beta*row[1]+beta*beta*row[2] for row in self.C]
    def transverse(self,alpha):return [sum(alpha**k*self.C[k][j] for k in range(3)) for j in range(3)]
    def solve_alpha(self,beta,target,qref,left=False):
        K=self.spatial(beta);a,q,r,p=qep(K,target,qref,left)
        r.update(alpha_r=float(a.real),alpha_i=float(a.imag),beta_r=float(beta),omega_r=self.omega,growth=float(-a.imag),iloc=self.i)
        if left:
            Pb=sum(a**k*(self.C[k][1]+2*beta*self.C[k][2]) for k in range(3));Pa=K[1]+2*a*K[2];den=np.vdot(p,Pa@q)
            ab=-np.vdot(p,Pb@q)/den;aw=-np.vdot(p,self.G@q)/den
            r.update(alpha_beta_r=float(ab.real),alpha_beta_i=float(ab.imag),alpha_omega_r=float(aw.real),alpha_omega_i=float(aw.imag))
        return a,q,r
    def validate(self):
        checks=[]
        for beta,w in [(-7.,50.),(16.,100.),(-35.,-100.),(.37,-50.)]:
            actual,receipt=export(self.case,self.i,beta,w);p=Pencil(self.case,self.i,w)
            errors=[float(norm(a-b))/max(1.,float(norm(a))) for a,b in zip(actual,p.spatial(beta))]
            assert max(errors)<2e-12,errors;checks.append(dict(beta=beta,omega=w,relative_matrix_errors=errors,export=receipt))
        return checks

def frequency_track(case,i,beta,a,q,target,step=.5):
    w=0.;records=[];last=None;minimum_step=1e-6
    while abs(w-target)>1e-12:
        dw=np.sign(target-w)*min(step,abs(target-w));p=Pencil(case,i,w)
        a,q,r=p.solve_alpha(beta,a,q,True)
        derivative=complex(r['alpha_omega_r'],r['alpha_omega_i'])
        try:aa,qq,rr=Pencil(case,i,w+dw).solve_alpha(beta,a+dw*derivative,q,True)
        except Exception as e:
            step/=2
            if step<minimum_step:return dict(finished=False,stopped_at=w,reason=repr(e),records=records),a,q
            continue
        w+=dw;a,q=aa,qq;records.append(rr)
        if rr['raw_MAC']>.999:step=min(2.,step*1.15)
    return dict(finished=True,omega=target,records=records),a,q

if __name__=='__main__':
    snapshot();results=[]
    for case,c in CASES.items():
        i=c['seed_station'];p=Pencil(case,i);validation=p.validate();pilot=next(r for r in json.loads((PREVIOUS/'outputs/pilot.json').read_text()) if r['case']==case)
        seed=pilot['spatial'];a0=complex(seed['alpha_r'],seed['alpha_i']);q0=np.load(PREVIOUS/'outputs'/f'{case}_seed.npz')['q']
        result=dict(case=case,station=i,validations=validation,tracks=[])
        for sign in [1,-1]:
            for omega in [50.,100.]:
                beta=-7.*sign;a=a0 if sign==1 else -np.conj(a0);q=q0 if sign==1 else q0.conj()
                track,a,q=frequency_track(case,i,beta,a,q,omega)
                path=ROOT/'outputs'/f'{case}_w{int(omega)}_seed_b{int(beta)}.npz';np.savez(path,q=q)
                track.update(beta=beta,vector_file=str(path));result['tracks'].append(track)
                print('PILOT',case,omega,beta,track['finished'],a,flush=True)
        results.append(result);(ROOT/'outputs/pilot.json').write_text(json.dumps(results,indent=2)+'\n')
