"""Identify new-angle M3 seeds with physical normal-grid and velocity-basis transfer."""
import csv,json,sys
from pathlib import Path
import h5py,numpy as np
ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/neutral_m3_sweep_consistent_20261028')
sys.path.insert(0,str(ROOT/'work'))
from common import CASES,snapshot
from solver import Pencil
from core import mac,qep

def profile(case,station):
    Pencil(case,station,0.)
    report=json.loads((ROOT/'native_pencils'/case/f'i{station:03d}/coefficients.json').read_text());export=report['exports'][0];directory=Path(export['files'][0]['path']).parent.parent
    with (directory/'wall_local_profile_metadata.csv').open() as f:meta=next(r for r in csv.DictReader(line for line in f if not line.startswith('#')) if int(r['iloc'])==export.get('local_station',station))
    tangent=np.array([float(meta['t_x']),float(meta['t_y'])]);normal=np.array([float(meta['n_x']),float(meta['n_y'])]);basis=np.row_stack([tangent,normal])
    with h5py.File(ROOT/'cases'/case/'in/grid_vertical2d_N.h5') as h:g=h['grid'][:,station-1,:,0]
    n=(g-g[0])@normal;n[0]=0.;assert np.all(np.diff(n)>0)
    return n,basis

def main():
    path=ROOT/'outputs/angle_mode_pilot.json';rows=json.loads(path.read_text());snapshot()
    for case in ['s50','s55','s65','s70','s73']:
        if any(r['case']==case for r in rows):continue
        angle=CASES[case]['angle'];source=min([r for r in rows if r['case']!='s75'],key=lambda r:(abs(CASES[r['case']]['angle']-angle),CASES[r['case']]['angle']>angle));i=CASES[case]['seed_station'];na,ba=profile(source['case'],source['station']);nb,bb=profile(case,i)
        q=np.load(source['vector_file']).reshape((-1,5));reference=np.column_stack([np.interp(nb,na,q[:,k].real,right=0)+1j*np.interp(nb,na,q[:,k].imag,right=0) for k in range(5)]);reference[:,1:3]=reference[:,1:3]@ba@bb.T;reference=reference.ravel();a=complex(source['alpha_r'],source['alpha_i'])
        trace=[];direct_failure=None
        try:
            aa,qq,r=Pencil(case,i,0.).solve_alpha(-7.,a,reference,True);assert mac(reference,qq)>.95
        except Exception as e:
            direct_failure=repr(e);Kold=Pencil(source['case'],source['station'],0.).spatial(-7.);Knew=Pencil(case,i,0.).spatial(-7.);dK=[y-x for x,y in zip(Kold,Knew)];lam=0.;step=.025;q=np.load(source['vector_file'])
            while lam<1-1e-14:
                target=min(1.,lam+step);K=[x+lam*y for x,y in zip(Kold,dK)];a,q,rr,left=qep(K,a,q,True);den=np.vdot(left,(K[1]+2*a*K[2])@q);der=-np.vdot(left,sum(a**k*dK[k] for k in range(3))@q)/den
                try:
                    aa,qq,r,_=qep([x+target*y for x,y in zip(Kold,dK)],a+(target-lam)*der,q,True);assert mac(q,qq)>.97
                except Exception as err:
                    step/=2
                    if step<1e-6:raise RuntimeError('Mode homotopy unresolved: '+repr(err))
                    continue
                trace.append(dict(lambda_value=target,alpha_r=aa.real,alpha_i=aa.imag,adjacent_raw_MAC=mac(q,qq),right_residual=r['right_residual']));lam=target;a,q=aa,qq;step=min(.08,step*1.2)
            aa,qq,r=Pencil(case,i,0.).solve_alpha(-7.,a,q,True)
        overlap=mac(reference,qq);assert r['right_residual']<1e-8
        vector=ROOT/'outputs'/f'{case}_mode_seed.npy';np.save(vector,qq)
        rows.append(dict(case=case,station=i,alpha_r=aa.real,alpha_i=aa.imag,beta_r=-7.,omega_r=0.,growth=-aa.imag,vector_file=str(vector),finished=True,receipt=r,identification='Nearest available non-75-degree actual M3 eigenfunction transferred in physical normal coordinates; local u/v basis rotated, then original-polynomial correction',source_case=source['case'],source_station=source['station'],transferred_vector_MAC=overlap,direct_transfer_failure=direct_failure,operator_homotopy_trace=trace,operator_homotopy_is_identification_only=True,source_basis=ba.tolist(),target_basis=bb.tolist()))
        tmp=path.with_suffix('.pilot.tmp');tmp.write_text(json.dumps(rows,indent=2)+'\n');tmp.replace(path);print('ADDITIONAL_PILOT',case,'from',source['case'],'station',i,'alpha',aa,'MAC',overlap,'residual',r['right_residual'],flush=True)
    snapshot()
if __name__=='__main__':main()
