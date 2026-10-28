"""Bounded process parallelism for independent native matrix exports."""
import argparse,concurrent.futures,json,os,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'work'))
from common import CASES,snapshot
from solver import Pencil
def one(job):
    case,i=job;p=Pencil(case,i,0.);return i
def main(case):
    assert case!='s75', '75 degrees must only reuse existing data'
    assert (ROOT/'outputs/stencil_validation.json').exists()
    count=CASES[case]['stations'];order=list(range(400,count+1))+list(range(1,min(400,count+1)))
    jobs=[(case,i) for i in order if not (ROOT/'native_pencils'/case/f'i{i:03d}/frequency_coefficient.json').exists()]
    start=time.monotonic()
    with concurrent.futures.ProcessPoolExecutor(max_workers=int(os.environ.get("M3_PRECOMPUTE_WORKERS","2"))) as pool:
        for n,i in enumerate(pool.map(one,jobs),1):
            if n%100==0:print('PREP',case,n,len(jobs),i,round(time.monotonic()-start,1),flush=True)
    snapshot();(ROOT/'outputs'/f'{case}_matrix_preparation.json').write_text(json.dumps(dict(case=case,stations=count,new_exports=len(jobs),seconds=time.monotonic()-start),indent=2)+'\n');print('PREP_DONE',case,count,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case');a=p.parse_args();main(a.case)
