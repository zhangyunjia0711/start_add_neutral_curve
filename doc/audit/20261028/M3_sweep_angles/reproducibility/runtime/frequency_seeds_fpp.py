"""Continue native neutral roots in physical frequency at several x anchors.

An unsuccessful fixed-x continuation is recorded as unresolved, never stable.
The bordered two-real-equation correction avoids division by alpha projection.
"""
import argparse,json,sys
from pathlib import Path
import numpy as np
ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/neutral_m3_sweep_consistent_20261028')
sys.path.insert(0,str(ROOT/'work'))
from common import CASES,snapshot
from solver import Pencil
from core import qep,mac
from contour3d import bordered,SCALE
FACTOR=2*np.pi*1000*.035/1810.6307411853716
def track(case,rec,target):
    i=int(round(rec['point']['xi']));alpha=rec['point']['alpha_r'];beta=rec['point']['beta_r'];q=np.load(rec['vector_file']);pencil=Pencil(case,i,0.)
    aa,q,r,p=qep(pencil.spatial(beta),complex(alpha),q,True);w=0.;step=.35;history=[];failure=None
    C0=pencil.baseC;G=pencil.G;zero=[[C*0 for C in row] for row in C0]
    def solve(omega,a,b,q,p):
        z=np.array([i,a,b])/SCALE;q0=q.copy();iters=0
        C=[list(row) for row in C0];C[0][0]=C[0][0]+omega*G
        for it in range(20):
            e,J,t,v,left,res,lr=bordered(C,zero,z,q,p)
            if abs(e)<1e-12 and res<1e-8:
                assert mac(q0,v)>.95;return float(z[1]*SCALE[1]),float(z[2]*SCALE[2]),v,left,res,it
            A=J[:,1:];weights=1/np.maximum(np.linalg.norm(A,axis=1),1e-30)
            dz=np.linalg.solve(A*weights[:,None],-np.array([e.real,e.imag])*weights)
            if np.linalg.norm(dz)>.1:dz*=.1/np.linalg.norm(dz)
            z[1:]+=dz
        raise RuntimeError('Neutral frequency corrector did not converge')
    while abs(w-target)>1e-12:
        new=w+np.sign(target-w)*min(step,abs(target-w))
        try:a,b,v,l,res,it=solve(new,alpha,beta,q,p)
        except Exception as e:
            step/=2
            if step<1e-5:failure=repr(e);break
            continue
        history.append(dict(omega=new,alpha_r=a,beta_r=b,right_residual=res,raw_MAC=mac(q,v)));w=new;alpha,beta,q,p=a,b,v,l
        if it<=4:step=min(.75,step*1.2)
    return dict(finished=bool(abs(w-target)<1e-12),station=i,omega_reached=float(w),history=history,failure=failure),alpha,beta,q
def main(case):
    guide=json.loads((ROOT/'outputs'/f'{case}_f0.0kHz_S_v3_guide.json').read_text());assert guide['status'] in ['closed_verified','open_at_domain_boundaries']
    records=[r for r in guide['records'] if abs(r['point']['xi']-round(r['point']['xi']))<1e-9]
    path=ROOT/'outputs'/f'{case}_neutral_seed_scan.json';data=json.loads(path.read_text());frequencies=sorted([r['f_kHz'] for r in json.loads((ROOT/'outputs/frequency_conversion.json').read_text())['frequencies']],key=lambda f:(abs(f),f));attempts=[]
    for f in frequencies:
        if any(r['f_kHz']==f and r['side']=='lower' for r in data['records']):continue
        xs=[-7.,-5.,-3.,-15.,-25.,0.,25.,100.,-30.,-33.,-34.5] if f>0 else [-7.,-3.,0.,5.,25.,-15.,-25.,50.,-30.,-33.,-34.5]
        done=False;used=[]
        for x in xs:
            distance=min(abs(r['point']['x_mm']-x) for r in records)
            candidates=[r for r in records if abs(abs(r['point']['x_mm']-x)-distance)<1e-8]
            candidates=sorted(candidates,key=lambda r:r['point']['beta_r'])
            for rec in candidates[:1]:
                if rec['index'] in used:continue
                used.append(rec['index']);trial,a,b,q=track(case,rec,f*FACTOR);trial.update(case=case,f_kHz=f,anchor_index=rec['index']);attempts.append(trial)
                if trial['finished']:
                    vf=ROOT/'outputs'/f'{case}_f{float(f)}_frequency_seed.npy';np.save(vf,q)
                    r=dict(case=case,f_kHz=float(f),station=trial['station'],alpha_r=a,alpha_i=0.,beta_r=b,beta_i=0.,omega_r=f*FACTOR,omega_i=0.,growth=0.,side='lower',kind='native_seed',vector_file=str(vf),frequency_seed_audit_index=len(attempts)-1)
                    data['records'].append(r);tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(data,indent=2)+'\n');tmp.replace(path);done=True;break
            if done:break
        audit=ROOT/'outputs'/f'{case}_frequency_seed_search.json';audit.write_text(json.dumps(dict(case=case,attempts=attempts,interpretation='Failure to continue a fixed-x seed is unresolved; no global stability claim.'),indent=2)+'\n')
        print('FREQUENCY_SEED',case,f,done,len(used),flush=True)
    snapshot()
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case');a=p.parse_args();main(a.case)
