"""Re-solve every original station crossing of the verified M3 closed contour."""
import argparse,concurrent.futures,datetime,json,math
from pathlib import Path
import numpy as np
from common import ROOT,sha
from fixed_beta_pencil import OperatorLinearFamily
FAMILY=None

def initialize():
    global FAMILY
    FAMILY=OperatorLinearFamily()

def solve(job):
    path=ROOT/'outputs/native_points'/('i%03d_p%d.json'%(job['station'],job['ordinal']))
    if path.exists():return json.loads(path.read_text())
    pencil=FAMILY.pencil(job['station']);q=np.load(job['vector_file']);alpha=job['alpha_guess'];omega=job['omega_guess'];history=[]
    result=dict(job=job,status='unknown',trace=history)
    try:
        for k in range(18):
            r,v=pencil.evaluate(alpha,q,omega);history.append(r);q=v;omega=r['omega_r']
            if abs(r['omega_i'])<1e-11 and r['original_real_frequency_polynomial_residual']<1e-8:
                assert abs(alpha-job['alpha_guess'])<max(.3,2*job['segment_alpha_span'])
                assert r['omega_alpha_i']*job['derivative_sign']>0
                vector=ROOT/'work/native_point_vectors'/path.with_suffix('.npy').name;np.save(vector,q)
                result.update(status='converged',point=dict(r,station=job['station'],x_mm=job['x_mm'],alpha_i=0.,
                    original_native_operator=True,beta_r=-7.),vector_file=str(vector),vector_sha256=sha(vector));break
            if abs(r['omega_alpha_i'])<1e-12:raise RuntimeError('Near native turning point; not a scalar alpha root')
            step=-r['omega_i']/r['omega_alpha_i'];alpha+=max(-1.,min(1.,step))
        else:raise RuntimeError('Native temporal-neutral scalar iteration limit')
    except Exception as error:result['reason']=repr(error)
    path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');return result

def main():
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=4);a=p.parse_args()
    arc=ROOT/'outputs/M3b7_closed_arc_v2.json';validation=ROOT/'outputs/closed_contour_validation.json'
    d=json.loads(arc.read_text());v=json.loads(validation.read_text());assert v['closure']['passed']
    folds={x['arc_segment'][0]:x for x in v['folds']};path=[]
    for item in d['records']:
        path.append(item)
        if item['index'] in folds:
            f=folds[item['index']];path.append(dict(index='fold',point=f['point'],vector_file=f['vector_file']))
    grid=np.array(json.loads((ROOT.parent/'neutral_beta_full_20261022/outputs/seed_catalog.json').read_text())['grids']['main_841']['x_nondim'])*35
    jobs=[];known={}
    for segment,(x,y) in enumerate(zip(path,path[1:])):
        ax,bx=x['point']['xi'],y['point']['xi']
        if abs(ax-bx)<1e-10:continue
        for station in range(math.ceil(min(ax,bx)-1e-10),math.floor(max(ax,bx)+1e-10)+1):
            theta=(station-ax)/(bx-ax);alpha=(1-theta)*x['point']['alpha_r']+theta*y['point']['alpha_r']
            if any(abs(alpha-old)<1e-5 for old in known.get(station,[])):continue
            ordinal=len(known.setdefault(station,[]));known[station].append(alpha)
            q=x if abs(station-ax)<abs(station-bx) else y
            jobs.append(dict(station=station,ordinal=ordinal,alpha_guess=alpha,
                omega_guess=(1-theta)*x['point']['omega_r']+theta*y['point']['omega_r'],
                vector_file=q['vector_file'],arc_segment=segment,arc_endpoints=[x['index'],y['index']],
                derivative_sign=1 if q['point']['omega_alpha_i']>0 else -1,
                segment_alpha_span=abs(y['point']['alpha_r']-x['point']['alpha_r']),x_mm=float(grid[station-1])))
    for f in [ROOT/'outputs/native_points',ROOT/'work/native_point_vectors']:f.mkdir(exist_ok=True)
    output=ROOT/'outputs/native_neutral_points.json';report=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),
        arc_sha256=sha(arc),validation_sha256=sha(validation),jobs=jobs,results=[],finished=False)
    with concurrent.futures.ProcessPoolExecutor(max_workers=a.workers,initializer=initialize) as pool:
        for index,r in enumerate(pool.map(solve,jobs)):
            report['results'].append(dict(station=r['job']['station'],ordinal=r['job']['ordinal'],status=r['status'],
                file=str(ROOT/'outputs/native_points'/('i%03d_p%d.json'%(r['job']['station'],r['job']['ordinal'])))))
            if (index+1)%50==0 or r['status']!='converged':
                temp=output.with_suffix('.tmp');temp.write_text(json.dumps(report,indent=2)+'\n');temp.replace(output)
                print('NATIVE',index+1,'/',len(jobs),'last',r['job']['station'],r['status'],flush=True)
    report['finished']=True;output.write_text(json.dumps(report,indent=2)+'\n')
    print('NATIVE_FINAL',len(jobs),'unknown',sum(r['status']!='converged' for r in report['results']),flush=True)

if __name__=='__main__':main()
