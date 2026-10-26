import concurrent.futures,json,time
from solver import coefficients
from common import ROOT,CASES,snapshot
def one(job):coefficients(*job);return job
if __name__=='__main__':
    snapshot();jobs=[(case,i) for case,c in CASES.items() for i in range(1,c['stations']+1) if i!=c['seed_station']];start=time.monotonic()
    with concurrent.futures.ProcessPoolExecutor(max_workers=6) as pool:
        for n,item in enumerate(pool.map(one,jobs),1):
            if n%100==0:print('PREP',n,len(jobs),item,round(time.monotonic()-start,1),flush=True)
    (ROOT/'outputs/preparation.json').write_text(json.dumps(dict(stations=len(jobs),seed_stations_prepared_in_pilot=True,seconds=time.monotonic()-start,inputs=snapshot()),indent=2)+'\n')
