import concurrent.futures,json
import numpy as np
from common import ROOT,CASES,snapshot
from survey import track_beta

def one(item):
    case,f,b,end=item;pilot=next(p for p in json.loads((ROOT/'outputs/physical_pilot.json').read_text()) if p['case']==case)
    track=next(t for t in pilot['tracks'] if t['f_kHz']==f and t['beta']==b);r=track['records'][-1];q=np.load(track['vector_file'])['q']
    result=track_beta(case,track['omega'],b,complex(r['alpha_r'],r['alpha_i']),q,end,step=.2);result['f_kHz']=f
    path=ROOT/'outputs'/f'{case}_f{f}kHz_b{int(b)}_to{int(end)}_survey.json';path.write_text(json.dumps(result,indent=2)+'\n')
    print('PHYSICAL_SURVEY',case,f,b,end,result['reached_end'],'MAXG',max(r['growth'] for r in result['records']),'ROOTS',[(q['beta_r'],q['alpha_r']) for q in result['roots']],flush=True)
    return dict(path=str(path),case=case,f_kHz=f,omega=track['omega'],seed_beta=b,end=end,reached_end=result['reached_end'],maximum_growth=max(r['growth'] for r in result['records']),root_count=len(result['roots']))

if __name__=='__main__':
    snapshot();jobs=[(case,f,b,end) for case in CASES for f in [50,100] for b in [-7.,7.] for end in [-60.,60.]]
    with concurrent.futures.ProcessPoolExecutor(max_workers=8) as pool:results=list(pool.map(one,jobs))
    (ROOT/'outputs/physical_survey_summary.json').write_text(json.dumps(results,indent=2)+'\n')
