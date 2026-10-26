"""Independent signed-frequency checks and local neutral-boundary checks."""
import concurrent.futures,json
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import norm
from common import ROOT,snapshot,sha
from solver import Pencil,original,mac
from physical import Pair

def one(job):
    case,f,c=job;tag=f'{case}_f{f}kHz_{c}'
    path=ROOT/'outputs'/f'{tag}_final.json';data=json.loads(path.read_text())
    guide=ROOT/'outputs'/f'{tag}_v3_guide.json'
    assert data['status']=='closed_verified'
    pts=data['points'];assert pts[0]['x_mm']==pts[-1]['x_mm'] and pts[0]['beta_r']==pts[-1]['beta_r']
    native=[p for p in pts if p['kind']==0];assert len(native)==data['native_point_count']
    counts={i:sum(p['xi']==i for p in native) for i in {p['xi'] for p in native}}
    assert all(n>=2 and n%2==0 for n in counts.values())
    # The same native physical root must not be emitted twice in a closed loop.
    for i in counts:
        rows=[p for p in native if p['xi']==i]
        for j,a in enumerate(rows):
            for b in rows[j+1:]:
                assert abs(a['alpha_r']-b['alpha_r'])>1e-6 or abs(a['beta_r']-b['beta_r'])>1e-6
    checks=[];sides=[]
    picks=[native[len(native)//4],native[3*len(native)//4]]+data['physical_tips']
    for k,p in enumerate(picks):
        q=np.load(p['vector_file']);state=dict(p);state.update(omega_r=-p['omega_r'],beta_r=-p['beta_r'],alpha_r=-p['alpha_r'])
        source=case if p['kind']==0 else Path(p['generation']).parent
        station=int(p['xi']) if p['kind']==0 else 5
        actual,v,receipt=original(source,station,state,state['beta_r'],f'{tag}_negative_check{k}')
        difference=abs(complex(actual['alpha_r'],actual['alpha_i'])-complex(state['alpha_r'],state['alpha_i']))
        overlap=mac(v,q.conj());assert difference<1e-7 and overlap>1-1e-7 and abs(actual['growth'])<1e-8
        checks.append(dict(kind=p['kind'],xi=p['xi'],x_mm=p['x_mm'],beta_r=state['beta_r'],omega_r=state['omega_r'],alpha_difference=difference,raw_MAC=overlap,growth=actual['growth'],receipt=receipt))
    # A regular native point on each side of the loop: central finite difference
    # checks the adjoint derivative and opposite signs across the neutral root.
    for p in picks[:2]:
        pencil=Pencil(case,int(p['xi']),p['omega_r']);a=complex(p['alpha_r'],p['alpha_i']);q=np.load(p['vector_file'])
        aa,q,r=pencil.solve_alpha(p['beta_r'],a,q,True);h=1e-3
        am,_,rm=pencil.solve_alpha(p['beta_r']-h,aa,q)
        ap,_,rp=pencil.solve_alpha(p['beta_r']+h,aa,q)
        fd=(ap-am)/(2*h);adj=complex(r['alpha_beta_r'],r['alpha_beta_i'])
        error=abs(fd-adj)/max(1.,abs(adj));assert error<1e-4 and am.imag*ap.imag<0
        sides.append(dict(xi=p['xi'],beta_r=p['beta_r'],delta_beta=h,growth_minus=-am.imag,growth_plus=-ap.imag,derivative_relative_error=error,raw_MAC_min=min(rm['raw_MAC'],rp['raw_MAC'])))
    data.update(final_guide=dict(path=str(guide),sha256=sha(guide)),negative_frequency_checks=checks,local_boundary_checks=sides,native_counts={str(k):v for k,v in sorted(counts.items())})
    path.write_text(json.dumps(data,indent=2)+'\n')
    out=dict(case=case,f_kHz=f,component=c,points=len(pts),x_range_mm=data['x_range_mm'],beta_range=data['negative_beta_range'],native_counts=sorted(set(counts.values())),negative_checks=len(checks),maximum_growth_error=data['maximum_growth_error'],maximum_right_residual=data['maximum_right_residual'],minimum_negative_MAC=min(p['raw_MAC'] for p in checks),maximum_negative_alpha_difference=max(p['alpha_difference'] for p in checks))
    print('AUDIT',json.dumps(out),flush=True);return out

if __name__=='__main__':
    snapshot();jobs=[(case,f,c) for case in ['bluntplate','rough'] for f in [50,100] for c in ['A','B']]
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:results=list(pool.map(one,jobs))
    (ROOT/'outputs/physical_audit.json').write_text(json.dumps(dict(results=results,protected_files=snapshot()),indent=2)+'\n')
