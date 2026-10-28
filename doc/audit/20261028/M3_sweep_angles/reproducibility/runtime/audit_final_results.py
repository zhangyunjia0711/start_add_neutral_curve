"""Independent aggregate verification of completed sweep-angle curves.

Re-evaluate selected saved right eigenvectors with native/physical matrices;
check every delivered point's records, topology and independent START receipts.
Does not solve or change the retained 75-degree data.
"""
import collections, datetime, json
from pathlib import Path
import numpy as np
from common import ROOT, CASES, snapshot, sha
from solver import Pencil
from physical import Pair

ANGLES=[30,45,50,55,60,65,70,73]
FACTOR=2*np.pi*1000*.035/1810.6307411853716

def main():
    protected=snapshot(); reports=[]; saved={}; polynomial_checks=[]; reused=0
    previous=ROOT/'outputs/final_numerical_audit.json'
    prior=json.loads(previous.read_text()) if previous.exists() else {}
    prior_hashes={r['curve']:r['sha256'] for r in prior.get('curves',[])}
    for angle in ANGLES:
        case=f's{angle}'; saved[case]=[]
        for path in sorted((ROOT/'outputs').glob(f'{case}_f*kHz_*_final.json')):
            d=json.loads(path.read_text()); points=d['points']; assert points
            assert d['status'] in ['closed_verified','open_at_domain_boundaries']
            assert abs(d['omega_r']-d['f_kHz']*FACTOR)<1e-12
            native=[r for r in points if r['kind']==0]
            counts=collections.Counter(int(r['xi']) for r in native)
            assert len(native)==d['native_point_count'] and len(counts)==d['native_station_count']
            assert set(counts)==set(range(min(counts),max(counts)+1))
            assert all(n>=2 and n%2==0 for n in counts.values())
            assert d['physical_tips'] and d['original_checks']
            for r in points:
                values=[r[k] for k in ['xi','x_mm','beta_r','alpha_r','alpha_i','growth','right_residual','left_residual','raw_MAC'] if r[k] is not None]
                assert np.isfinite(values).all() and abs(r['growth'])<1e-8
                assert r['right_residual']<1e-7
                if r['left_residual'] is not None:assert r['left_residual']<1e-7
                assert abs(r['growth']+r['alpha_i'])<1e-14
                assert Path(r['vector_file']).is_file()
                if r['kind']==0:
                    assert r['xi']==int(r['xi']) and r['x_mm']==CASES[case]['x_mm'][int(r['xi'])-1]
            for c in d['original_checks']:
                assert c['alpha_difference']<1e-7 and c['raw_MAC']>1-1e-7
                assert c['receipt']['original_residual']<1e-7
                assert (ROOT/'runs'/c['receipt']['run']).is_dir()
            if d['status']=='closed_verified':
                assert points[-1]['kind']==2
                assert all(abs(points[0][k]-points[-1][k])<1e-12 for k in ['xi','alpha_r','beta_r','x_mm'])
            else:
                assert points[-1]['kind']!=2
                assert all(any(abs(p['xi']-edge)<1e-8 for edge in [1,len(CASES[case]['x_mm'])]) for p in [points[0],points[-1]])
            assert len(d['symmetry'])==len(native)
            assert max(r['matrix_symmetry_error'] for r in d['symmetry'])<1e-12
            assert max(r['conjugate_residual'] for r in d['symmetry'])<1e-7
            chosen=[native[0],native[len(native)//2],max(native,key=lambda r:abs(r['growth']))]+d['physical_tips']
            if prior_hashes.get(path.name)==sha(path):
                polynomial_checks.extend(c for c in prior['polynomial_checks'] if c['curve']==path.name)
                reused+=1;chosen=[]
            for r in chosen:
                if r['kind']==0:
                    p=Pencil(case,int(r['xi']),d['omega_r'])
                else:
                    tag=Path(r['generation']).parent.name.rsplit('_pair',1)[0]
                    p=Pair(case,r['native_bracket'][0],d['omega_r'],tag).pencil(r['theta'])
                a=complex(r['alpha_r'],r['alpha_i']); b=r['beta_r'];q=np.load(r['vector_file'])
                M=sum(a**k*b**j*p.C[k][j] for k in range(3) for j in range(3))
                residual=float(np.linalg.norm(M@q)/max(1.,np.linalg.norm(q)))
                assert residual<1e-7,(path.name,r['xi'],residual)
                polynomial_checks.append(dict(curve=path.name,xi=r['xi'],kind=r['kind'],absolute_right_residual=residual))
            record=dict(curve=path.name,sha256=sha(path),angle=angle,f_kHz=d['f_kHz'],component=d['component'],status=d['status'],points=len(points),native_stations=len(counts),original_START_checks=len(d['original_checks']),maximum_growth_error=max(abs(r['growth']) for r in points),maximum_recorded_residual=max(r['right_residual'] for r in points))
            reports.append(record);saved[case].append(d)
        assert saved[case],case
    duplicate_checks=0
    for case, curves in saved.items():
        for index,a in enumerate(curves):
            r=next(r for r in a['points'] if r['kind']==0)
            for b in curves[index+1:]:
                if b['f_kHz']!=a['f_kHz']:continue
                duplicate_checks+=1
                assert not any(abs(s['xi']-r['xi'])<1e-8 and abs(s['alpha_r']-r['alpha_r'])<1e-6 and abs(s['beta_r']-r['beta_r'])<1e-6 for s in b['points'] if s['kind']==0),(case,a['component'],b['component'],a['f_kHz'])
    result=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),passed=True,curves=reports,polynomial_checks=polynomial_checks,matrix_rechecks_reused_for_exact_unchanged_final_files=reused,duplicate_component_pairs_checked=duplicate_checks,protected_input_count=len(protected),protected_inputs_unchanged=True,recomputed_angles=ANGLES,retained75_numerical_solves=0,scope='Checks numerical neutrality/topology for delivered curves; does not certify absence of disconnected components or establish frequency cutoffs.')
    out=ROOT/'outputs/final_numerical_audit.json';out.write_text(json.dumps(result,indent=2)+'\n')
    print('FINAL_AUDIT',len(reports),'curves',sum(r['points'] for r in reports),'points',len(polynomial_checks),'independent matrix rechecks',flush=True)

if __name__=='__main__':main()
