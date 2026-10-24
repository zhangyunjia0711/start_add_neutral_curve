"""Independent final data/source audit, including saved spatial side samples."""
import ast,csv,datetime,json
from pathlib import Path
from common import ROOT,SOLVER,SOLVER_HASH,sha,snapshot,table

def main():
    out=ROOT/'outputs'
    summary=json.loads((out/'M3_closed_beta_minus7_summary.json').read_text())
    with (out/'M3_closed_beta_minus7_curve.csv').open() as f:rows=list(csv.DictReader(f))
    assert len(rows)==1571
    assert [int(r['order']) for r in rows]==list(range(1571))
    assert all(float(r['beta_r'])==-7 and float(r['alpha_i'])==0 for r in rows)
    assert all(rows[0][k]==rows[-1][k] for k in ['x_mm','omega_r','beta_r','alpha_r','alpha_i'])
    kinds={k:sum(int(r['point_kind'])==k for r in rows) for k in [0,1,2]}
    assert kinds=={0:1564,1:6,2:1}
    assert all(Path(r['source_file']).is_file() for r in rows)
    native={(int(r['native_station']),float(r['alpha_r']),float(r['omega_r'])) for r in rows if int(r['point_kind'])==0}
    assert len(native)==1564
    stations={r[0] for r in native};assert stations==set(range(10,712))
    spatial=json.loads((out/'native_spatial_verification.json').read_text())
    side_checks=[]
    for item in spatial['results']:
        path=ROOT/'runs'/item['spatial_run']/'out/neutral_field.plt'
        points=sorted(table(path),key=lambda p:p['omega_r'])
        assert len(points)==3 and all(p['status']==0 and p['beta_r']==-7 for p in points)
        crossing=points[0]['growth']*points[-1]['growth']<0
        side_checks.append(dict(station=item['station'],ordinal=item['ordinal'],file=str(path),sha256=sha(path),
            omega=[p['omega_r'] for p in points],growth=[p['growth'] for p in points],strict_sign_crossing=crossing))
    for p in (ROOT/'work').glob('*.py'):ast.parse(p.read_text(),filename=str(p))
    original=snapshot();assert original==summary['original_inputs']
    assert sha(SOLVER)==SOLVER_HASH
    report=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),
        rows=1571,point_kind_counts=kinds,unique_native_points=1564,original_station_count=702,
        first_last_physical_coordinates_exact=True,metadata_kind_intentionally_differs=True,
        all_beta_strictly_minus7=True,no_operator_guide_points_exported=True,
        curve_max_abs_temporal_imag=max(float(r['temporal_imag_residual']) for r in rows),
        curve_max_original_polynomial_residual=max(float(r['polynomial_residual']) for r in rows),
        spatial_side_checks=side_checks,strict_side_sign_crossings=sum(q['strict_sign_crossing'] for q in side_checks),
        spatial_side_check_count=len(side_checks),source_inputs_unchanged=True,original_inputs=original,
        solver_path=str(SOLVER),solver_sha256=sha(SOLVER),python_syntax_passed=True)
    (out/'M3_closed_beta_minus7_audit.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:report[k] for k in ['rows','point_kind_counts','strict_side_sign_crossings','spatial_side_check_count','source_inputs_unchanged']}))

if __name__=='__main__':main()
