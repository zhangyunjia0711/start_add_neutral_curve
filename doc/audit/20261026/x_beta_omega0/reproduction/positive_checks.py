import json,numpy as np
from pathlib import Path
from common import ROOT,CASES,snapshot
from solver import original,mac

results=[]
for case,c in CASES.items():
    data=json.loads((ROOT/'outputs'/f'{case}_final.json').read_text());points=data['points']
    for label,row in [('nose_tip',data['physical_tips'][0]),('tail_high_abs_beta',points[0]),('tail_low_abs_beta',points[-1])]:
        if row['kind']==1:
            mode=str(Path(row['generation']).parent);station=5
        else:mode=case;station=int(row['xi'])
        state=dict(omega_r=0.,beta_r=-row['beta_r'],alpha_r=-row['alpha_r'],alpha_i=row['alpha_i'])
        fresh,q,receipt=original(mode,station,state,state['beta_r'],f'{case}_positive_{label}')
        expected=complex(-row['alpha_r'],row['alpha_i']);actual=complex(fresh['alpha_r'],fresh['alpha_i']);overlap=mac(np.load(row['vector_file']).conj(),q)
        assert abs(actual-expected)<1e-7 and overlap>1-1e-7 and abs(fresh['growth'])<1e-8
        results.append(dict(case=case,point=label,alpha_difference=abs(actual-expected),raw_MAC=overlap,growth=fresh['growth'],receipt=receipt))
(ROOT/'outputs/positive_beta_checks.json').write_text(json.dumps(dict(checks=results,inputs=snapshot()),indent=2)+'\n')
print('POSITIVE_CHECK_PASS',len(results),flush=True)
