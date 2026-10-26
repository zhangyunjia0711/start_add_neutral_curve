import concurrent.futures,json,subprocess,sys
from common import ROOT,CASES,driver
def one(job):
    case,f,c=job;name=f'{case}_f{f}kHz_{c}';guide=json.loads((ROOT/'outputs'/(name+'_v3_guide.json')).read_text())
    if guide['status'] not in ['closed_verified','open_at_domain_boundaries']:return dict(case=case,f_kHz=f,component=c,status=guide['status'],skipped=True)
    log=ROOT/'outputs'/(name+'_finalize.log')
    with log.open('w') as stream:p=subprocess.run([sys.executable,str(ROOT/'work/finalize.py'),case,str(f),c],env=driver.ENV,stdout=stream,stderr=subprocess.STDOUT)
    path=ROOT/'outputs'/(name+'_final.json');data=json.loads(path.read_text()) if path.exists() else None
    result=dict(case=case,f_kHz=f,component=c,exit=p.returncode,points=len(data['points']) if data else 0,log=str(log));print('FINALIZE',json.dumps(result),flush=True);return result
if __name__=='__main__':
    jobs=[(case,f,c) for case in CASES for f in [50,100] for c in ['A','B']]
    with concurrent.futures.ProcessPoolExecutor(max_workers=6) as pool:results=list(pool.map(one,jobs))
    (ROOT/'outputs/finalize_summary.json').write_text(json.dumps(results,indent=2)+'\n')
