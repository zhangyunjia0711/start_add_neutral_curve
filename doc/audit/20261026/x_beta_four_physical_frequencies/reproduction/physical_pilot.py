"""User-confirmed cyclic frequencies in kHz, converted to solver angular omega."""
import json,numpy as np
from common import ROOT,CASES,snapshot
from solver import Pencil,PREVIOUS,frequency_track
if __name__=='__main__':
    snapshot();conversion=json.loads((ROOT/'outputs/frequency_conversion.json').read_text());results=[]
    for case,c in CASES.items():
        prior=next(p for p in json.loads((PREVIOUS/'outputs/pilot.json').read_text()) if p['case']==case);a0=complex(prior['spatial']['alpha_r'],prior['spatial']['alpha_i']);q0=np.load(PREVIOUS/'outputs'/f'{case}_seed.npz')['q']
        result=dict(case=case,station=c['seed_station'],tracks=[])
        for frequency in conversion['frequencies']:
            if frequency['f_kHz']<0:continue
            for sign in [1,-1]:
                beta=-7.*sign;a=a0 if sign==1 else -np.conj(a0);q=q0 if sign==1 else q0.conj();w=frequency['omega_nd']
                track,a,q=frequency_track(case,c['seed_station'],beta,a,q,w,step=.2);assert track['finished'],track
                path=ROOT/'outputs'/f"{case}_f{frequency['f_kHz']}kHz_b{int(beta)}_seed.npz";np.savez(path,q=q)
                track.update(f_kHz=frequency['f_kHz'],beta=beta,vector_file=str(path));result['tracks'].append(track)
                print('PHYSICAL_PILOT',case,frequency['f_kHz'],beta,w,a,flush=True)
        results.append(result);(ROOT/'outputs/physical_pilot.json').write_text(json.dumps(results,indent=2)+'\n')
