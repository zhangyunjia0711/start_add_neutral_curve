import concurrent.futures, datetime, json, time
from common import ROOT,CASES,snapshot
from solver import Pencil

def one(item):
    case,i=item;p=Pencil(case,i)
    return case,i
if __name__=='__main__':
    snapshot();start=time.monotonic();jobs=[(case,i) for case,c in CASES.items() for i in range(1,c['stations']+1)]
    with concurrent.futures.ProcessPoolExecutor(max_workers=6) as pool:
        for n,item in enumerate(pool.map(one,jobs,chunksize=1),1):
            if n%50==0:print('PREP',n,len(jobs),item,round(time.monotonic()-start,1),flush=True)
    (ROOT/'outputs/preparation.json').write_text(json.dumps(dict(stations=len(jobs),seconds=time.monotonic()-start,completed=datetime.datetime.now().isoformat(),original_inputs=snapshot()),indent=2)+'\n')
