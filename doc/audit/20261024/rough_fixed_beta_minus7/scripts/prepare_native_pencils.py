"""Export native fixed-beta pencils in parallel; native physics is unmodified."""
import argparse,concurrent.futures,json,datetime
from common import ROOT,snapshot
from fixed_beta_pencil import Pencil
from profile_interpolation import AXIS

def main():
    p=argparse.ArgumentParser();p.add_argument('--first',type=int,default=1);p.add_argument('--last',type=int,default=361)
    p.add_argument('--workers',type=int,default=6);a=p.parse_args();snapshot()
    out=ROOT/'outputs/native_pencil_preparation.json'
    result=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),
        first=a.first,last=a.last,results=[],finished=False)
    def one(i):
        try:
            q=Pencil('M3',i,AXIS[i],f'M3b7_native_matrix_i{i:03d}')
            return dict(station=i,passed=True,exports=q.exports)
        except Exception as error:return dict(station=i,passed=False,error=repr(error))
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
        for r in pool.map(one,range(a.first,a.last+1)):
            result['results'].append(r);temp=out.with_suffix('.tmp');temp.write_text(json.dumps(result,indent=2)+'\n');temp.replace(out)
            if r['station']%25==0 or not r['passed']:print('MATRIX_NATIVE',r['station'],r['passed'],flush=True)
    result['finished']=True;result['original_inputs']=snapshot();out.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
