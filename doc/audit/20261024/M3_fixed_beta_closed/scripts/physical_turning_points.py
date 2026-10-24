"""Re-converge six turning points on actual interpolated local CFD profiles.

Operator-linear contour points are predictors only. Final turning points use
the frozen original solver's matrices assembled from interpolated n/primitive
profiles, with endpoint coefficient parity checked against native H5.
"""
import concurrent.futures,datetime,json
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import norm
from common import ROOT,sha
from fixed_beta_pencil import Pencil
from profile_interpolation import Pair,mac

def solve(item):
    index,source=item;out=ROOT/'outputs'/('physical_fold_%d.json'%index)
    if out.exists():return json.loads(out.read_text())
    pair_indices=source['point']['native_bracket'];pair=Pair(*pair_indices)
    theta=source['point']['theta'];alpha=source['point']['alpha_r'];omega=source['point']['omega_r'];q=np.load(source['vector_file'])
    cache={};history=[]
    def ev(t,a,ref,w):
        assert 0<=t<=1
        key=float(t).hex()
        state=dict(omega_r=w,beta_r=-7.,alpha_r=a,alpha_i=0.)
        if key not in cache:
            directory=pair.build(t);cache[key]=Pencil(directory,5,state,'M3b7_physical_fold%d_%s'%(index,key.replace('.','').replace('+','p').replace('-','m')))
        p=cache[key];r,v=p.evaluate(a,ref,w)
        r.update(theta=t,xi=pair.left+t,x_mm=json.loads((pair.build(t)/'generation.json').read_text())['x_mm'],alpha_i=0.)
        return r,v
    for k in range(18):
        r,q=ev(theta,alpha,q,omega);omega=r['omega_r'];F=np.array([r['omega_i'],r['omega_alpha_i']]);history.append(dict(point=r,F=F.tolist()))
        if max(abs(F))<1e-10:break
        ht=min(.002,theta/3,(1-theta)/3);ha=1e-4
        rp,_=ev(theta+ht,alpha,q,omega);rm,_=ev(theta-ht,alpha,q,omega)
        ap,_=ev(theta,alpha+ha,q,omega);am,_=ev(theta,alpha-ha,q,omega)
        J=np.array([[(rp['omega_i']-rm['omega_i'])/(2*ht),r['omega_alpha_i']],
                    [(rp['omega_alpha_i']-rm['omega_alpha_i'])/(2*ht),(ap['omega_alpha_i']-am['omega_alpha_i'])/(2*ha)]])
        delta=np.linalg.solve(J,-F);scale=min(1,.2/max(abs(delta[0]),1e-100),.5/max(abs(delta[1]),1e-100))
        theta+=scale*delta[0];alpha+=scale*delta[1]
    assert max(abs(F))<1e-8 and r['original_real_frequency_polynomial_residual']<1e-7
    checks=[]
    for t,station in [(0.,pair.left),(1.,pair.right)]:
        state=dict(omega_r=omega,beta_r=-7.,alpha_r=alpha,alpha_i=0.)
        native=Pencil('M3',station,state,f'M3b7_fold{index}_endpoint_native{station}')
        derived=Pencil(pair.build(t),5,state,f'M3b7_fold{index}_endpoint_derived{station}')
        relative=[float(norm(x-y))/max(1,float(norm(x))) for x,y in zip(native.A+[native.M],derived.A+[derived.M])]
        assert max(relative)<1e-11,('Endpoint coefficient mismatch',station,relative)
        checks.append(dict(station=station,relative_coefficient_errors=relative,passed=True))
    path=ROOT/'work'/('M3b7_physical_fold_q%d.npy'%index);np.save(path,q)
    result=dict(index=index,point=r,vector_file=str(path),vector_sha256=sha(path),trace=history,
        operator_predictor=source['point'],delta_x_mm=r['x_mm']-source['point']['x_mm'],
        delta_omega=r['omega_r']-source['point']['omega_r'],MAC_to_predictor=mac(q,np.load(source['vector_file'])),
        endpoint_coefficient_checks=checks,model='Linear interpolation of adjacent native frozen-local n and primitive profiles; actual new LST matrices',passed=True)
    out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');return result

def main():
    source=ROOT/'outputs/closed_contour_validation.json';d=json.loads(source.read_text())
    with concurrent.futures.ProcessPoolExecutor(max_workers=3) as pool:
        results=list(pool.map(solve,enumerate(d['folds'])))
    report=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),
        source_sha256=sha(source),results=results,finished=True,all_passed=all(r['passed'] for r in results))
    (ROOT/'outputs/physical_turning_points.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print('PHYSICAL_FOLDS',[(r['point']['xi'],r['point']['x_mm'],r['point']['omega_r'],r['delta_x_mm']) for r in results])

if __name__=='__main__':main()
