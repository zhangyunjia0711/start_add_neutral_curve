"""Track real neutral alpha/beta in signed physical frequency on native profiles."""
import concurrent.futures,json,sys
from pathlib import Path
import numpy as np
from scipy.sparse import bmat,csc_matrix
from scipy.sparse.linalg import splu
ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/neutral_xbeta_extended_frequencies_20261027')
sys.path.insert(0,str(ROOT/'work'))
from common import CASES,snapshot,sha
from solver import Pencil,qep,mac
SOURCE=ROOT.parent/'neutral_xbeta_four_omega_20261026/outputs'
S=np.array([20.,10.]);W=6.0727867021443425/50

def start(case,sign):
    d=json.loads((SOURCE/f'{case}_f100kHz_{"A" if sign>0 else "B"}_final.json').read_text())
    pts=[p for p in d['points'] if p['kind']==0]
    stations=sorted(set(int(p['xi']) for p in pts))
    if sign>0:i=min(stations,key=lambda i:abs(CASES[case]['x_mm'][i-1]))
    else:
        i=max(stations,key=lambda i:np.ptp([p['beta_r'] for p in pts if p['xi']==i]))
    rows=sorted([p for p in pts if p['xi']==i],key=lambda p:sign*p['beta_r'])
    row=rows[0];q=np.load(row['vector_file']);a=row['alpha_r'];b=row['beta_r']
    if sign<0:a=-a;b=-b;q=q.conj()
    return i,a,b,q,str(SOURCE/f'{case}_f100kHz_{"A" if sign>0 else "B"}_final.json')

def solve(p,z,q0,p0):
    a,b=z*S;C=p.C;Q=sum(a**k*b**j*C[k][j] for k in range(3) for j in range(3))
    n=len(q0);M=bmat([[Q,csc_matrix(p0[:,None])],[csc_matrix(q0.conj()[None,:]),csc_matrix((1,1),dtype=complex)]],format='csc')
    lu=splu(M);rhs=np.r_[np.zeros(n,dtype=complex),1.+0j];v=lu.solve(rhs);w=lu.solve(rhs,trans='H')
    A=sum(b**j*(C[1][j]+2*a*C[2][j]) for j in range(3))
    B=sum(a**k*(C[k][1]+2*b*C[k][2]) for k in range(3))
    g=np.array([-np.vdot(w[:-1],D@v[:-1])*scale for D,scale in zip([A,B],S)])
    gw=-np.vdot(w[:-1],p.G@v[:-1])*W
    q=v[:-1]/np.linalg.norm(v[:-1]);left=w[:-1]/np.linalg.norm(w[:-1])
    return v[-1],np.array([g.real,g.imag]),np.array([gw.real,gw.imag]),q,left,float(np.linalg.norm(Q@q))

def one(job):
    case,sign=job;i,a,b,q,source=start(case,sign);f=100.*sign
    p=Pencil(case,i,W*f);_,q,_,left=qep(p.spatial(b),a,q,True);z=np.array([a,b])/S
    result=dict(case=case,sign=sign,station=i,x_mm=CASES[case]['x_mm'][i-1],source=source,source_sha256=sha(source),records=[],failures=[],status='running')
    path=ROOT/'outputs'/f'{case}_{"positive" if sign>0 else "negative"}_frequency_extension.json'
    outdir=ROOT/'work'/f'{case}_{sign}_frequency_vectors';outdir.mkdir(exist_ok=True)
    step=1.;target=200.*sign
    while abs(target-f)>1e-9:
        df=sign*min(step,abs(target-f));p=Pencil(case,i,W*f)
        e,J,gw,_,_,_=solve(p,z,q,left)
        try:
            pred=z+np.linalg.solve(J,-gw*df);trial=pred.copy();new=Pencil(case,i,W*(f+df));tr=[]
            for iteration in range(25):
                err,jac,_,v,w,res=solve(new,trial,q,left);tr.append(dict(iteration=iteration,border_abs=float(abs(err)),residual=res,alpha_r=float(trial[0]*20),beta_r=float(trial[1]*10)))
                if abs(err)<1e-12 and res<1e-8:break
                weights=1/np.maximum(np.linalg.norm(jac,axis=1),1e-30);delta=np.linalg.solve(jac*weights[:,None],-np.array([err.real,err.imag])*weights)
                if np.linalg.norm(delta)>.15:delta*=.15/np.linalg.norm(delta)
                trial+=delta
            assert abs(err)<1e-12 and res<1e-8 and mac(q,v)>.97,(f+df,tr[-1],mac(q,v))
            assert np.linalg.norm(trial-z)<.5
        except Exception as error:
            result['failures'].append(dict(f_kHz=f,attempt_df=df,reason=repr(error)));step/=2
            if step<1e-3:result['status']='fixed_station_frequency_limit_or_unresolved';break
            continue
        overlap=mac(q,v);f+=df;z=trial;q,left=v,w;file=outdir/f'q{len(result["records"]):04d}.npy';np.save(file,q)
        result['records'].append(dict(case=case,f_kHz=f,station=i,omega_r=W*f,alpha_r=float(z[0]*20),alpha_i=0.,beta_r=float(z[1]*10),growth=0.,right_residual=res,vector_file=str(file),vector_sha256=sha(file),raw_MAC=overlap,trace=tr))
        path.write_text(json.dumps(result,indent=2)+'\n')
        if abs(f-round(f/25)*25)<1e-8:print('FREQUENCY',case,f,i,z*S,flush=True)
        if len(tr)<=4:step=min(1.,step*1.25)
    else:result['status']='completed_to_200kHz'
    result['stopped_frequency_kHz']=f;path.write_text(json.dumps(result,indent=2)+'\n')
    print('SEED_EXTENSION',case,sign,result['status'],f,len(result['records']),flush=True);return result

def main():
    snapshot()
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:results=list(pool.map(one,[(case,sign) for case in CASES for sign in [-1,1]]))
    (ROOT/'outputs/seed_extension_summary.json').write_text(json.dumps(dict(paths=results),indent=2)+'\n')
if __name__=='__main__':main()
