"""Final numerical audit; partial searches remain explicitly partial."""
import concurrent.futures,json,sys
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import norm
ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/neutral_m3_frequency_limits_20261027')
WIN=Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate/75°')
sys.path.insert(0,str(ROOT/'work'));sys.path.insert(0,str(WIN/'work/M3_frequency_limits_20261027'))
from common import snapshot,sha,CASES
from solver import Pencil,mac
from physical import Pair
from publish_limits import load,W
def one(source):
    file=Path(source['path']);d=json.loads(file.read_text());points=d['points'];assert d['status']=='closed_verified';assert points[-1]['kind']==2 and all(points[0][k]==points[-1][k] for k in ['xi','x_mm','beta_r','alpha_r'])
    native=[p for p in points if p['kind']==0];stations={int(p['xi']) for p in native};tips=[p for p in points if p['kind']==1];assert len(tips)==2
    expected=set(range(int(np.ceil(min(p['xi'] for p in tips))),int(np.floor(max(p['xi'] for p in tips)))+1));assert stations==expected,('COVERAGE',file.name,sorted(expected-stations),sorted(stations-expected))
    for i in stations:
        rows=[p for p in native if int(p['xi'])==i];assert len(rows)==2 and abs(rows[0]['beta_r']-rows[1]['beta_r'])>1e-7
    for p in points:assert abs(p['growth'])<1e-8 and p['right_residual']<1e-7
    overlaps=[mac(np.load(a['vector_file']),np.load(b['vector_file'])) for a,b in zip(points[:-1],points[1:])];assert min(overlaps)>.9,('Adjacent vector discontinuity',file.name,min(overlaps))
    original=d['original_checks'];assert original
    for r in original:assert r['alpha_difference']<1e-7 and r.get('MAC',r.get('raw_MAC',0))>1-1e-7
    checks=[]
    for j in sorted(set(np.linspace(0,len(native)-1,min(16,len(native)),dtype=int).tolist())):
        p=native[j];q=np.load(p['vector_file']);a=complex(p['alpha_r'],p['alpha_i']);K=Pencil(source['case'],int(p['xi']),p['omega_r']).spatial(p['beta_r']);res=float(np.linalg.norm(sum(a**k*K[k] for k in range(3))@q));assert res<1e-7;checks.append(dict(xi=p['xi'],beta=p['beta_r'],residual=res))
    return dict(**source,points=len(points),native_stations=len(stations),additional_profile_points=sum(p['kind']==3 for p in points),maximum_growth_error=max(abs(p['growth']) for p in points),maximum_right_residual=max(p['right_residual'] for p in points),minimum_adjacent_MAC=min(overlaps),matrix_checks=checks,original_checks=original)
def main():
    protected=snapshot();curves,limits,sources=load()
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:results=list(pool.map(one,[r for r in sources if not r['prior']]))
    critical=[]
    for r in limits:
        d=r['record'];assert d['status']=='actual_profile_extremum_verified'
        for t in d['frequency_side_checks']:assert (t['maximum_growth']>0) if t['side']=='inside' else (t['maximum_growth']<0)
        q=np.load(d['vector_file']);assert sha(d['vector_file'])==d['vector_sha256']
        critical.append(dict(case=r['case'],sign=r['sign'],f_kHz=r['f_kHz'],omega_nd=r['omega_nd'],x_mm=r['x_mm'],beta=r['beta'],right_residual=r['right_residual'],original_check=d['original_check'],frequency_side_checks=d['frequency_side_checks']))
    coverage=[]
    for case in CASES:
        negative=[]
        for kind in ['negative_native','negative_secondary_native','negative_tail']:
            d=json.loads((ROOT/'outputs'/f'{case}_{kind}_fold_curve.json').read_text());negative.extend(d['records'])
        stations={r['station'] for r in negative};extreme=min(negative,key=lambda r:r['f_kHz']);coverage.append(dict(case=case,negative_native_stations=len(stations),first_station=min(stations),last_station=max(stations),missing_stations=sorted(set(range(1,CASES[case]['stations']+1))-stations),negative_native_extreme=extreme['f_kHz']))
        for kind in ['secondary','tail']:
            d=json.loads((ROOT/'outputs'/f'{case}_negative_{kind}_arc.json').read_text());assert d['status'] in ['connected_to_known_branch','domain_boundary'];assert d['minimum_frequency']['point']['f_kHz']>=min(r['f_kHz'] for r in negative)-1e-4
    endpoints=[]
    for case,indices in [('rough',[41,42,43]),('bluntplate',[207,208,215,216])]:
        for i in indices:
            pair=Pair(case,i,0.,case+'_positive_frequency_limit')
            for t,j in [(0.,i),(1.,i+1)]:
                p=pair.pencil(t);native=Pencil(case,j,264.5*W);p.C=[list(c) for c in p.baseC];p.C[0][0]=p.C[0][0]+264.5*W*p.G;errors=[float(norm(a-b))/max(1.,float(norm(a))) for a,b in zip(native.spatial(18.),p.spatial(18.))];assert max(errors)<1e-9;endpoints.append(dict(case=case,pair=i,theta=t,station=j,relative_matrix_errors=errors))
    outside=json.loads((ROOT/'outputs/positive_outside_global_samples.json').read_text());outside_summary=[dict(case=r['case'],station=r['station'],x_mm=r['x_mm'],target_f_kHz=r['target_f_kHz'],status=r['status'],maximum_computed_growth=r.get('maximum_growth'),directions=r.get('directions')) for r in outside];assert all(r['maximum_computed_growth']['growth']<0 for r in outside_summary)
    out=dict(curves=results,limits=critical,negative_native_coverage=coverage,physical_endpoint_checks=endpoints,positive_outside_samples=outside_summary,outside_scope='Supplementary mode-continuation checks in beta=[-64,64], including explicit unresolved projection segments; not an exhaustive eigenmode/global spectral proof',scope='Tracked M3 complex representative and resolved connected pockets; no separately generated conjugate copy',protected_inputs=protected)
    (ROOT/'outputs/final_numerical_audit.json').write_text(json.dumps(out,indent=2)+'\n');print('AUDIT_COMPLETE',len(results),'new curves',sum(r['points'] for r in results),'new points',len(critical),'limits',flush=True)
if __name__=='__main__':main()
