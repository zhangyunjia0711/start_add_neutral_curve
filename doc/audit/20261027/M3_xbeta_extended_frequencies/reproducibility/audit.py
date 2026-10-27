"""Audit new actual contours, original equations, coverage and mode continuity."""
import concurrent.futures,json,sys
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import norm
ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/neutral_xbeta_extended_frequencies_20261027')
sys.path.insert(0,str(ROOT/'work'))
from common import CASES,snapshot,sha
from solver import Pencil,mac
FS=[-112,-110,125,150,175,200]
def one(job):
    case,f=job;tag=f'{case}_f{f}kHz_S';file=ROOT/'outputs'/f'{tag}_final.json';d=json.loads(file.read_text());guidefile=ROOT/'outputs'/f'{tag}_v3_guide.json';guide=json.loads(guidefile.read_text())
    assert d['status']==guide['status']=='closed_verified' and guide['closure']['passed']
    points=d['points'];assert points[-1]['kind']==2 and all(points[0][k]==points[-1][k] for k in ['xi','x_mm','alpha_r','beta_r'])
    native=[p for p in points if p['kind']==0];stations=set(int(p['xi']) for p in native)
    tips=d['physical_tips'];expected=set(range(int(np.ceil(min(t['xi'] for t in tips))),int(np.floor(max(t['xi'] for t in tips)))+1))
    assert stations==expected,('COVERAGE',case,f,sorted(expected-stations),sorted(stations-expected))
    counts={i:sum(p['xi']==i for p in native) for i in stations};assert all(n>=2 and n%2==0 for n in counts.values())
    for i in stations:
        rows=[p for p in native if p['xi']==i]
        for j,p in enumerate(rows):
            for q in rows[j+1:]:assert abs(p['beta_r']-q['beta_r'])>1e-6 or abs(p['alpha_r']-q['alpha_r'])>1e-6
    checks=[]
    for index in sorted(set(np.linspace(0,len(native)-1,16,dtype=int).tolist())):
        p=native[index];q=np.load(p['vector_file']);K=Pencil(case,int(p['xi']),p['omega_r']).spatial(p['beta_r']);a=complex(p['alpha_r'],p['alpha_i']);res=float(np.linalg.norm(sum(a**k*K[k] for k in range(3))@q));assert res<1e-7
        checks.append(dict(index=index,xi=p['xi'],beta=p['beta_r'],original_residual=res,scaled_residual=res/max(1.,sum(abs(a)**k*norm(K[k]) for k in range(3)))))
    seed=next(s for s in json.loads((ROOT/'outputs/neutral_seeds.json').read_text())['records'] if s['case']==case and s['f_kHz']==f);i=seed['station'];rows=sorted([p for p in native if p['xi']==i],key=lambda p:p['beta_r']);assert len(rows)==2
    boundary=[]
    for side,p in zip(['lower','upper'],rows):
        pencil=Pencil(case,i,p['omega_r']);a=complex(p['alpha_r'],p['alpha_i']);q=np.load(p['vector_file']);aa,q,r=pencil.solve_alpha(p['beta_r'],a,q,True);h=1e-3
        am,qm,rm=pencil.solve_alpha(p['beta_r']-h,aa,q);ap,qp,rp=pencil.solve_alpha(p['beta_r']+h,aa,q);fd=(ap-am)/(2*h);der=complex(r['alpha_beta_r'],r['alpha_beta_i']);error=abs(fd-der)/max(1.,abs(der))
        assert error<1e-4 and am.imag*ap.imag<0
        assert (am.imag>0>ap.imag) if side=='lower' else (am.imag<0<ap.imag)
        boundary.append(dict(side=side,xi=i,beta=p['beta_r'],growth_minus=-am.imag,growth_plus=-ap.imag,derivative_error=error,minimum_MAC=min(rm['raw_MAC'],rp['raw_MAC'])))
    assert all(abs(p['growth'])<1e-8 and p['right_residual']<1e-7 for p in points)
    assert d['original_checks'] and all(c['alpha_difference']<1e-7 and c['raw_MAC']>1-1e-7 for c in d['original_checks'])
    assert guide['closure']['raw_MAC']>1-1e-7
    result=dict(case=case,f_kHz=f,status=d['status'],final_file=str(file),final_sha256=sha(file),guide_file=str(guidefile),guide_sha256=sha(guidefile),points=len(points),native_points=len(native),native_stations=len(stations),x_range_mm=d['x_range_mm'],beta_range=d['negative_beta_range'],maximum_growth_error=max(abs(p['growth']) for p in points),maximum_right_residual=max(p['right_residual'] for p in points),closure=guide['closure'],minimum_adjacent_guide_MAC=min(p['point']['raw_MAC'] for p in guide['records']),matrix_checks=checks,boundary_checks=boundary,physical_tips=[{k:p[k] for k in ['xi','x_mm','beta_r','alpha_r','original_receipt']} for p in tips],original_checks=d['original_checks'])
    print('AUDIT',case,f,len(points),d['x_range_mm'],flush=True);return result
def main():
    snapshot()
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:results=list(pool.map(one,[(case,f) for case in CASES for f in FS]))
    paths=[]
    for case in CASES:
        for sign in ['negative','positive']:
            file=ROOT/'outputs'/f'{case}_{sign}_frequency_extension.json';d=json.loads(file.read_text());records=d['records'];overlaps=[mac(np.load(a['vector_file']),np.load(b['vector_file'])) for a,b in zip(records[:-1],records[1:])];assert min(overlaps)>.97
            paths.append(dict(path=str(file),sha256=sha(file),status=d['status'],stopped_frequency_kHz=d['stopped_frequency_kHz'],minimum_adjacent_MAC=min(overlaps)))
    (ROOT/'outputs/new_frequency_audit.json').write_text(json.dumps(dict(results=results,frequency_paths=paths,protected_inputs=snapshot()),indent=2)+'\n')
    print('AUDIT_COMPLETE',len(results),flush=True)
if __name__=='__main__':main()
