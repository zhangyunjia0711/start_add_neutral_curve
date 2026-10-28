"""70-degree guide from actual local-profile matrices; native/fold evidence stays exact.

Adaptive cubic matrix approximations are guides only. Each segment is checked
at an independent physical-profile midpoint; final points use native START
matrices or the actual local-profile interpolant, never this approximation.
"""
import argparse,concurrent.futures,fcntl,hashlib,json,pickle,time
from functools import lru_cache
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import norm
from common import ROOT,CASES,sha
from solver import Pencil
from physical import Pair
TOL=5e-7
NODES=np.array([0.,.25,.75,1.]);INVERSE=np.linalg.inv(np.vander(NODES,4,increasing=True))

def build_cell(case,i):
    assert case=='s70';dest=ROOT/'physical_guides'/case/f'i{i:03d}';dest.mkdir(parents=True,exist_ok=True);report=dest/'receipt.json';cache=dest/'cell.pkl'
    with (dest/'build.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if report.exists():return dest
        pair=Pair(case,i,0.,f'physical_guide_{case}_i{i:03d}');samples={};checks=[]
        def sample(theta):
            key=float(theta)
            if key in samples:return samples[key]
            if key in [0.,1.]:p=Pencil(case,i+int(key),0.)
            else:p=pair.pencil(key)
            samples[key]=[p.baseC[k][j] for k in range(3) for j in range(3)]+[p.G]
            return samples[key]
        def segment(lo,hi,depth):
            values=[sample(lo+(hi-lo)*t) for t in NODES];power=[[sum(INVERSE[p,n]*values[n][j] for n in range(4)) for j in range(10)] for p in range(4)]
            mid=sample((lo+hi)/2);fit=[sum(.5**p*power[p][j] for p in range(4)) for j in range(10)];errors=[float(norm(a-b))/max(1.,float(norm(a))) for a,b in zip(mid,fit)];checks.append(dict(lo=lo,hi=hi,depth=depth,midpoint_coefficient_errors=errors,accepted=max(errors)<=TOL))
            if max(errors)>TOL:
                assert depth<7,('Physical guide needs further refinement',case,i,lo,hi,errors)
                return segment(lo,(lo+hi)/2,depth+1)+segment((lo+hi)/2,hi,depth+1)
            return [dict(lo=lo,hi=hi,power=power,maximum_midpoint_error=max(errors))]
        segments=segment(0.,1.,0);tmp=cache.with_suffix('.tmp')
        with tmp.open('wb') as f:pickle.dump(segments,f,protocol=4)
        tmp.replace(cache);data=dict(case=case,left_station=i,right_station=i+1,model='Adaptive cubic guide of actual matched-j local primitive-profile interpolants; native nodes exact; not CFD refinement',independent_midpoint_tolerance=TOL,segments=len(segments),physical_samples=sorted(samples),checks=checks,cache_sha256=sha(cache));tmp=report.with_suffix('.tmp');tmp.write_text(json.dumps(data,indent=2)+'\n');tmp.replace(report)
    return dest

@lru_cache(maxsize=8)
def cell(case,i):
    dest=build_cell(case,i)
    with (dest/'cell.pkl').open('rb') as f:return pickle.load(f)

def evaluate(case,i,theta,omega):
    segments=cell(case,i);s=next(q for q in segments if q['lo']-1e-12<=theta<=q['hi']+1e-12);width=s['hi']-s['lo'];u=(theta-s['lo'])/width;power=s['power'];v=[sum(u**p*power[p][j] for p in range(4)) for j in range(10)];dv=[sum(p*u**(p-1)*power[p][j] for p in range(1,4))/width for j in range(10)]
    if abs(theta)<1e-12 or abs(theta-1)<1e-12:
        native=Pencil(case,i+int(theta>.5),0.);v=[native.baseC[k][j] for k in range(3) for j in range(3)]+[native.G]
    C=[[v[3*k+j] for j in range(3)] for k in range(3)];dC=[[dv[3*k+j] for j in range(3)] for k in range(3)];C[0][0]=C[0][0]+omega*v[9];dC[0][0]=dC[0][0]+omega*dv[9];return C,dC

def endpoint_check():
    out=[]
    for i in [1,25,50,90,150]:
        pair=Pair('s70',i,0.,f'physical_guide_s70_i{i:03d}')
        for theta in [0.,1.]:
            a=pair.pencil(theta);b=Pencil('s70',i+int(theta),0.);errors=[float(norm(a.baseC[k][j]-b.baseC[k][j]))/max(1.,float(norm(b.baseC[k][j]))) for k in range(3) for j in range(3)]+[float(norm(a.G-b.G))/max(1.,float(norm(b.G)))];assert max(errors)<2e-12;out.append(dict(left_station=i,theta=theta,coefficient_errors=errors))
    (ROOT/'outputs/s70_physical_guide_endpoint_validation.json').write_text(json.dumps(dict(passed=True,checks=out),indent=2)+'\n');print('PHYSICAL_GUIDE_ENDPOINTS_PASSED',flush=True)

def main():
    endpoint_check();start=time.monotonic();count=CASES['s70']['stations']-1
    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as pool:
        jobs=[pool.submit(build_cell,'s70',i) for i in range(1,count+1)]
        for n,future in enumerate(concurrent.futures.as_completed(jobs),1):
            future.result()
            if n%20==0:print('PHYSICAL_GUIDE_PREP',n,count,round(time.monotonic()-start,1),flush=True)
    print('PHYSICAL_GUIDE_COMPLETE',count,round(time.monotonic()-start,1),flush=True)
if __name__=='__main__':main()
