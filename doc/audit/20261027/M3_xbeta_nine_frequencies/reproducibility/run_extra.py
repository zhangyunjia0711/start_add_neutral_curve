"""Compute added-frequency neutral seeds, guides and corrected delivered contours."""
import concurrent.futures,datetime,json,os,shutil,subprocess,sys
from pathlib import Path
ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/neutral_xbeta_extra_frequencies_20261027')
SOURCE=ROOT.parent/'neutral_frequency_connection_20261026_f100/outputs/connection_report.json'
sys.path.insert(0,str(ROOT/'work'))
from common import CASES,driver,snapshot,sha
from solver import Pencil,mac
import numpy as np
FS=[-75,-25,25,75]

def seeds():
    assert not (ROOT/'outputs/neutral_seeds.json').exists()
    data=json.loads(SOURCE.read_text());out=[]
    for path in data['paths']:
        case,side,i=path['case'],path['side'],path['station']
        for f in FS:
            r=min(path['records'],key=lambda r:abs(r['f_kHz']-f));q=np.load(r['vector_file']);assert sha(r['vector_file'])==r['vector_sha256']
            beta=r['beta'];a=complex(r['alpha_r'],r['alpha_i']);p=Pencil(case,i,6.0727867021443425*f/50)
            for it in range(12):
                aa,v,s=p.solve_alpha(beta,a,q,True)
                assert mac(q,v)>.95
                a,q=aa,v
                if abs(a.imag)<1e-11:break
                beta-=a.imag/s['alpha_beta_i']
            assert abs(a.imag)<1e-8
            filename=ROOT/'outputs'/f'{case}_{side}_f{f:+d}_seed.npy';np.save(filename,q)
            out.append(dict(case=case,side=side,station=i,f_kHz=f,omega_r=p.omega,beta_r=beta,alpha_r=a.real,alpha_i=a.imag,
                vector_file=str(filename),vector_sha256=sha(filename),right_residual=s['right_residual'],source_frequency=r['f_kHz'],source_vector=r['vector_file'],source_sha256=r['vector_sha256']))
            print('SEED',case,side,f,beta,a,flush=True)
    (ROOT/'outputs/neutral_seeds.json').write_text(json.dumps(dict(source=str(SOURCE),source_sha256=sha(SOURCE),records=out),indent=2)+'\n')

def one(job):
    case,f=job;tag=f'{case}_f{f}kHz_S';result=dict(case=case,f_kHz=f,stages=[])
    for script,suffix in [('contour3d.py','v3_guide'),('finalize.py','final')]:
        path=ROOT/'outputs'/f'{tag}_{suffix}.json';log=ROOT/'outputs'/f'{tag}_{suffix}.log'
        if not path.exists():
            with log.open('w') as stream:
                p=subprocess.run([sys.executable,str(ROOT/'work'/script),case,str(f),'S'],env=driver.ENV,stdout=stream,stderr=subprocess.STDOUT)
            code=p.returncode
        else:code=0
        d=json.loads(path.read_text()) if path.exists() else {}
        stage=dict(script=script,exit=code,status=d.get('status','NO_RESULT'),records=len(d.get('records',d.get('points',[]))),log=str(log))
        result['stages'].append(stage);print('STAGE',json.dumps(dict(case=case,f_kHz=f,**stage)),flush=True)
        if code or d.get('status') not in ['closed_verified','open_at_domain_boundaries']:break
    return result

def main():
    snapshot()
    shutil.copy2(Path(__file__),ROOT/'work/run_extra.py')
    if not (ROOT/'outputs/neutral_seeds.json').exists():seeds()
    jobs=[(case,f) for case in CASES for f in FS]
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:
        results=list(pool.map(one,jobs))
    (ROOT/'outputs/run_summary.json').write_text(json.dumps(dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),results=results,protected_inputs=snapshot()),indent=2)+'\n')
    print('RUN_COMPLETE',[(r['case'],r['f_kHz'],r['stages'][-1]['status']) for r in results],flush=True)
if __name__=='__main__':main()
