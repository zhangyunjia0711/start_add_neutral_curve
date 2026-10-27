"""Numerical audit of the eight additional directly computed contours."""
import concurrent.futures,json,sys
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import norm
ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/neutral_xbeta_extra_frequencies_20261027')
sys.path.insert(0,str(ROOT/'work'))
from common import snapshot,sha,CASES
from solver import Pencil,mac

def one(job):
    case,f=job;tag=f'{case}_f{f}kHz_S';file=ROOT/'outputs'/f'{tag}_final.json';d=json.loads(file.read_text());assert d['status']=='closed_verified'
    guidefile=ROOT/'outputs'/f'{tag}_v3_guide.json';guide=json.loads(guidefile.read_text())
    assert guide['closure']['passed'] and guide['closure']['raw_MAC']>1-1e-7
    points=d['points'];assert points[-1]['kind']==2 and points[0]['x_mm']==points[-1]['x_mm'] and points[0]['beta_r']==points[-1]['beta_r']
    native=[p for p in points if p['kind']==0];stations={int(p['xi']) for p in native}
    assert stations==set(range(min(stations),max(stations)+1))
    counts={i:sum(p['xi']==i for p in native) for i in stations};assert all(n%2==0 and n>=2 for n in counts.values())
    # Physical folds can expose new native crossings; require exact physical coverage.
    guide_folds=[r for r in d['physical_tips']]
    expected=set(range(int(np.ceil(min(r['xi'] for r in guide_folds))),int(np.floor(max(r['xi'] for r in guide_folds)))+1))
    assert stations==expected,('Native coverage differs from physical folds',tag,sorted(expected-stations),sorted(stations-expected))
    for r in guide_folds:
        if int(np.floor(r['xi']))!=int(np.floor(r['guide_xi'])):
            assert d.get('native_crossing_repairs'),('Physical fold crosses an unrepaired native station',tag,r['xi'],r['guide_xi'])
            assert all(len(t['records'])>=20 for repair in d['native_crossing_repairs'] for t in repair['traces'])
    for i in stations:
        rows=[p for p in native if p['xi']==i]
        for j,a in enumerate(rows):
            for b in rows[j+1:]:assert abs(a['alpha_r']-b['alpha_r'])>1e-6 or abs(a['beta_r']-b['beta_r'])>1e-6
    sources=[]
    for p in native:
        path=Path(p['vector_file']);assert path.exists();sources.append(dict(point_index=p['guide_order'],path=str(path),sha256=sha(path)))
    selected=sorted(set(np.linspace(0,len(native)-1,16,dtype=int).tolist()))
    matrix_checks=[]
    for index in selected:
        p=native[index];q=np.load(p['vector_file']);pencil=Pencil(case,int(p['xi']),p['omega_r']);alpha=complex(p['alpha_r'],p['alpha_i']);K=pencil.spatial(p['beta_r'])
        residual=float(np.linalg.norm(sum(alpha**k*K[k] for k in range(3))@q));assert residual<1e-7
        relative=residual/(sum(abs(alpha)**k*norm(K[k]) for k in range(3))*np.linalg.norm(q))
        matrix_checks.append(dict(xi=p['xi'],beta=p['beta_r'],absolute_original_residual=residual,coefficient_Frobenius_scaled_residual=relative))
    seed=CASES[case]['seed_station'];sides=sorted([p for p in native if p['xi']==seed],key=lambda p:p['beta_r']);assert len(sides)==2
    side_checks=[]
    for side,p in zip(['lower','upper'],sides):
        pencil=Pencil(case,seed,p['omega_r']);q=np.load(p['vector_file']);a,q,r=pencil.solve_alpha(p['beta_r'],complex(p['alpha_r'],p['alpha_i']),q,True);h=1e-3
        am,_,rm=pencil.solve_alpha(p['beta_r']-h,a,q);ap,_,rp=pencil.solve_alpha(p['beta_r']+h,a,q)
        derivative=complex(r['alpha_beta_r'],r['alpha_beta_i']);error=abs((ap-am)/(2*h)-derivative)/max(1,abs(derivative))
        assert error<1e-4 and am.imag*ap.imag<0
        if side=='lower':assert am.imag>0>ap.imag
        else:assert am.imag<0<ap.imag
        side_checks.append(dict(side=side,x_mm=p['x_mm'],beta=p['beta_r'],growth_minus=-am.imag,growth_plus=-ap.imag,derivative_relative_error=error,minimum_MAC=min(rm['raw_MAC'],rp['raw_MAC'])))
    checks=d['original_checks'];assert len(checks)>=8
    assert all(r['alpha_difference']<1e-7 and r['raw_MAC']>1-1e-7 for r in checks)
    result=dict(case=case,f_kHz=f,status=d['status'],final_file=str(file),final_sha256=sha(file),guide_file=str(guidefile),guide_sha256=sha(guidefile),
        closure=guide['closure'],points=len(points),native_points=len(native),native_stations=len(stations),native_counts=sorted(set(counts.values())),
        x_range_mm=d['x_range_mm'],beta_range=d['negative_beta_range'],maximum_growth_error=d['maximum_growth_error'],maximum_recorded_right_residual=d['maximum_right_residual'],
        minimum_adjacent_guide_MAC=min(r['point']['raw_MAC'] for r in guide['records']),matrix_checks=matrix_checks,local_boundary_checks=side_checks,
        original_independent_checks=checks,physical_tips=[{k:r[k] for k in ['xi','x_mm','beta_r','alpha_r','alpha_i','turning_condition','original_receipt']} for r in guide_folds],vector_files=sources)
    print('AUDITED',case,f,len(points),result['x_range_mm'],flush=True);return result

def main():
    snapshot()
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:results=list(pool.map(one,[(case,f) for case in CASES for f in [-75,-25,25,75]]))
    out=dict(results=results,protected_inputs=snapshot(),scope='Eight new signed-frequency curves; 16 original-matrix spot checks per curve plus local boundary sign tests and prior independent frozen START checks')
    (ROOT/'outputs/new_frequency_audit.json').write_text(json.dumps(out,indent=2)+'\n');print('AUDIT_COMPLETE',flush=True)
if __name__=='__main__':main()
