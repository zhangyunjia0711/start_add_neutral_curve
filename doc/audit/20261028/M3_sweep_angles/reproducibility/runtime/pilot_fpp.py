"""Identify M3 on current FPP fields using physical eigenfunction transfer."""
import json,sys
from pathlib import Path
import h5py,numpy as np
from scipy.sparse import load_npz
ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/neutral_m3_sweep_consistent_20261028')
sys.path.insert(0,str(ROOT/'work'))
from common import CASES,snapshot
from solver import Pencil
from core import qep,mac
OLD=ROOT.parent/'75°/neutral_xbeta_omega0_20261026'
REF=ROOT.parent/'75°/neutral_xbeta_four_omega_20261026/cases/bluntplate'
def normal_grid(folder):
    with h5py.File(folder/'in/grid_vertical2d_N.h5') as h:g=h['grid'][:,118,:,0]
    normal=np.array([g[0,0],g[0,1]]);normal/=np.linalg.norm(normal)
    n=(g-g[0])@normal;assert np.all(np.diff(n)>0);return n
def transfer(q,na,nb):
    q=q.reshape((-1,5));return np.column_stack([np.interp(nb,na,q[:,k].real,right=0)+1j*np.interp(nb,na,q[:,k].imag,right=0) for k in range(5)]).ravel()
def main():
    snapshot();pilot=next(r for r in json.loads((OLD/'outputs/pilot.json').read_text()) if r['case']=='bluntplate');r=pilot['spatial'];a=complex(r['alpha_r'],r['alpha_i']);q=np.load(OLD/'outputs/bluntplate_seed.npz')['q'];na=normal_grid(REF);rows=[]
    for case in ['s75','s60','s45','s30']:
        if case!='s75':
            trial=ROOT.parent/'neutral_m3_sweep_angles_20261027'
            source=next(v for v in json.loads((trial/'outputs/angle_mode_pilot.json').read_text()) if v['case']==case)
            a=complex(source['alpha_r'],source['alpha_i']);q=np.load(source['vector_file']);na=normal_grid(trial/'cases'/case)
        p=Pencil(case,119,0.);nb=normal_grid(ROOT/'cases'/case);reference=transfer(q,na,nb)
        try:
            aa,qq,r,_=qep(p.spatial(-7.),a,reference,True)
            method='Physical normal-coordinate eigenfunction transfer and original-polynomial correction'
        except Exception as e:
            raise RuntimeError('FPP physical mode transfer needs explicit repair at '+case+': '+repr(e))
        path=ROOT/'outputs'/f'{case}_mode_seed.npy';np.save(path,qq)
        rows.append(dict(case=case,station=119,alpha_r=aa.real,alpha_i=aa.imag,beta_r=-7.,omega_r=0.,growth=-aa.imag,vector_file=str(path),finished=True,receipt=r,identification=method,transferred_vector_MAC=mac(reference,qq)))
        (ROOT/'outputs/angle_mode_pilot.json').write_text(json.dumps(rows,indent=2)+'\n');print('FPP_PILOT',case,aa,r['right_residual'],mac(reference,qq),flush=True)
        a,q,na=aa,qq,nb
    snapshot()
if __name__=='__main__':main()
