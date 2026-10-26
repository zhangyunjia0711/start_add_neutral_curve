import concurrent.futures,json,subprocess,sys
from common import ROOT,CASES,driver
def one(job):
    case,f,component=job;name=f'{case}_f{f}kHz_{component}_v3';log=ROOT/'outputs'/(name+'_contour.log')
    existing=ROOT/'outputs'/(name+'_guide.json')
    if existing.exists():
        data=json.loads(existing.read_text());assert data['status']=='closed_verified',data['status']
        return dict(case=case,f_kHz=f,component=component,exit=0,status=data['status'],points=len(data['records']),reused=True)
    command=[sys.executable,str(ROOT/'work/contour3d.py'),case,str(f),component]
    with log.open('w') as stream:p=subprocess.run(command,env=driver.ENV,stdout=stream,stderr=subprocess.STDOUT)
    output=ROOT/'outputs'/(name+'_guide.json');data=json.loads(output.read_text()) if output.exists() else None
    result=dict(case=case,f_kHz=f,component=component,exit=p.returncode,status=data['status'] if data else 'NO_RESULT',points=len(data['records']) if data else 0,log=str(log))
    print('CONTOUR',json.dumps(result),flush=True);return result
if __name__=='__main__':
    jobs=[(case,f,c) for case in CASES for f in [50,100] for c in ['A','B']]
    with concurrent.futures.ProcessPoolExecutor(max_workers=8) as pool:results=list(pool.map(one,jobs))
    (ROOT/'outputs/contour_summary.json').write_text(json.dumps(results,indent=2)+'\n')
