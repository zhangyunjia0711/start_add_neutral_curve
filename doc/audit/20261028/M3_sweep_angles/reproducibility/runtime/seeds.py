"""Native angle seeds at five common physical frequencies."""
import argparse,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'work'))
from common import CASES,snapshot
from solver import Pencil,frequency_track
from core import mac
FREQUENCIES=[0.,-50.,50.,-100.,100.]
FACTOR=2*np.pi*1000*.035/1810.6307411853716
def one(case):
    pilot=next(r for r in json.loads((ROOT/'outputs/angle_mode_pilot.json').read_text()) if r['case']==case)
    a0=complex(pilot['alpha_r'],pilot['alpha_i']);q0=np.load(pilot['vector_file']);i=pilot['station'];records=[];failures=[];scans=[]
    for f in FREQUENCIES:
        omega=f*FACTOR;state,a,q=frequency_track(case,i,-7.,a0,q0,omega,step=.25)
        if not state['finished']:failures.append(dict(f_kHz=f,phase='frequency',state=state));continue
        pencil=Pencil(case,i,omega);roots=[]
        for direction,end in [(-1,-64.),(1,-.01 if f==0 else 40.)]:
            b=-7.;aa,qq,rr=pencil.solve_alpha(b,a,q,True);step=.5;series=[rr]
            while direction*(end-b)>1e-10:
                h=direction*min(step,abs(end-b));bb=b+h;der=complex(rr['alpha_beta_r'],rr['alpha_beta_i'])
                try:ab,qb,rb=pencil.solve_alpha(bb,aa+h*der,qq,True)
                except Exception as error:
                    step/=2
                    if step<1e-5:failures.append(dict(case=case,f_kHz=f,phase='beta_scan',beta=b,direction=direction,reason=repr(error)));break
                    continue
                if aa.imag*ab.imag<0:
                    lo=(b,aa,qq,rr);hi=(bb,ab,qb,rb)
                    for it in range(55):
                        mid=(lo[0]+hi[0])/2;ref=lo if abs(mid-lo[0])<abs(mid-hi[0])+1e-15 else hi
                        ar,qr,r=pencil.solve_alpha(mid,ref[1]+(mid-ref[0])*complex(ref[3]['alpha_beta_r'],ref[3]['alpha_beta_i']),ref[2],True)
                        if abs(ar.imag)<1e-10:break
                        if ar.imag*lo[1].imag>0:lo=(mid,ar,qr,r)
                        else:hi=(mid,ar,qr,r)
                    assert abs(ar.imag)<1e-8
                    if not any(abs(mid-v['beta_r'])<1e-5 for v in roots):
                        vfile=ROOT/'outputs'/f'{case}_f{f}_seed{len(roots)}.npy';np.save(vfile,qr)
                        r.update(case=case,f_kHz=f,station=i,vector_file=str(vfile),kind='native_seed',side='candidate');roots.append(r)
                b,aa,qq,rr=bb,ab,qb,rb;series.append(rr)
                if rr['raw_MAC']>.997:step=min(1.,step*1.2)
            scans.append(dict(case=case,f_kHz=f,direction=direction,series=series))
        if roots:
            lower=min(roots,key=lambda r:r['beta_r']);lower['side']='lower';records+=roots
        print('SEEDS',case,f,[(r['beta_r'],r['growth']) for r in roots],flush=True)
        path=ROOT/'outputs'/f'{case}_neutral_seed_scan.json';tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(dict(case=case,records=records,failures=failures,scans=scans),indent=2)+'\n');tmp.replace(path)
    snapshot()
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case');a=p.parse_args();one(a.case)
