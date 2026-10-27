"""Continue verified new-frequency seeds into closed neutral contours."""
import concurrent.futures,datetime,json,subprocess,sys
from pathlib import Path
import numpy as np
ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/neutral_xbeta_extended_frequencies_20261027')
sys.path.insert(0,str(ROOT/'work'))
from common import CASES,driver,snapshot,sha
from solver import Pencil,mac
FS=[-112,-110,125,150,175,200]

def prepare():
    path=ROOT/'outputs/frequency_conversion.json';conversion=json.loads(path.read_text())
    conversion['initial_requested_probe_frequencies_kHz']=[-200,-150,-125,125,150,200]
    conversion['frequencies']=[dict(f_kHz=f,omega_nd=6.0727867021443425*f/50) for f in FS]
    path.write_text(json.dumps(conversion,indent=2)+'\n')
    for name in ['contour3d.py','finalize.py']:
        path=ROOT/'work'/name;s=path.read_text()
        s=s.replace('choices=[-200, -150, -125, 125, 150, 200]','choices='+repr(FS))
        s=s.replace('choices=[-110, 125, 150, 175, 200]','choices='+repr(FS));path.write_text(s)
    sources=json.loads((ROOT/'outputs/seed_extension_summary.json').read_text())['paths'];records=[]
    for case in CASES:
        for f in FS:
            source=next(p for p in sources if p['case']==case and p['sign']==int(np.sign(f)))
            row=min(source['records'],key=lambda r:abs(r['f_kHz']-f));assert abs(row['f_kHz']-f)<1e-10
            q=np.load(row['vector_file']);assert sha(row['vector_file'])==row['vector_sha256']
            p=Pencil(case,row['station'],row['omega_r']);a,q,r=p.solve_alpha(row['beta_r'],row['alpha_r'],q,True)
            assert abs(a.imag)<1e-8 and r['raw_MAC']>1-1e-7
            out=ROOT/'outputs'/f'{case}_f{f:+d}_seed.npy';np.save(out,q)
            records.append(dict(row,side='lower',alpha_r=a.real,alpha_i=a.imag,right_residual=r['right_residual'],vector_file=str(out),vector_sha256=sha(out),frequency_path_source=str(ROOT/'outputs'/f'{case}_{"positive" if f>0 else "negative"}_frequency_extension.json')))
    (ROOT/'outputs/neutral_seeds.json').write_text(json.dumps(dict(records=records,selection='Direct signed-frequency continuation of the existing representative; negative trial at a fixed station is not a whole-domain no-root proof'),indent=2)+'\n')

def one(job):
    case,f=job;tag=f'{case}_f{f}kHz_S';result=dict(case=case,f_kHz=f,stages=[])
    for script,suffix in [('contour3d.py','v3_guide'),('finalize.py','final')]:
        path=ROOT/'outputs'/f'{tag}_{suffix}.json';log=ROOT/'outputs'/f'{tag}_{suffix}.log'
        if not path.exists():
            with log.open('w') as stream:p=subprocess.run([sys.executable,str(ROOT/'work'/script),case,str(f),'S'],env=driver.ENV,stdout=stream,stderr=subprocess.STDOUT)
            code=p.returncode
        else:code=0
        d=json.loads(path.read_text()) if path.exists() else {}
        result['stages'].append(dict(script=script,exit=code,status=d.get('status','NO_RESULT'),records=len(d.get('records',d.get('points',[]))),log=str(log)))
        print('STAGE',case,f,result['stages'][-1],flush=True)
        if code or d.get('status') not in ['closed_verified','open_at_domain_boundaries']:break
    return result

def main():
    snapshot()
    if not (ROOT/'outputs/neutral_seeds.json').exists():prepare()
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:results=list(pool.map(one,[(case,f) for case in CASES for f in FS]))
    (ROOT/'outputs/run_summary.json').write_text(json.dumps(dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),results=results,protected_inputs=snapshot()),indent=2)+'\n')
    print('COMPLETE',[(r['case'],r['f_kHz'],r['stages'][-1]['status']) for r in results],flush=True)
if __name__=='__main__':main()
