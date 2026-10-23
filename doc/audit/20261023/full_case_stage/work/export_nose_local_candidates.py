#!/usr/bin/env python3
"""Collect independently accepted local-nose roots and explicit unknown intervals."""
import csv,json
from full_driver import ROOT,sha,snapshot
from local_peak_x_continue import stamp

SOURCES=[
 ('M1p_signed','local_peak_nose_M1p_v10_progress.json'),
 ('M1p_signed','local_peak_nose_M1p_j8_v10_1616_progress.json'),
 ('M1p_signed','local_peak_nose_M1p_ridge_v10_1616_progress.json'),
 ('M1c','local_peak_nose_M1c_v10_1616_progress.json'),
 ('M1c','local_peak_nose_M1c_j8_v10_1616_progress.json'),
 ('M1c','local_peak_nose_M1c_ridge_v10_1616_progress.json'),
 ('M1c','nose_M1c_i11_fine_refine_20261023.json'),
]

def main():
    snapshot();results={};provenance=[]
    for family,name in SOURCES:
        path=ROOT/'outputs'/name;d=json.loads(path.read_text())
        provenance.append(dict(source=family,path=str(path),sha256=sha(path)))
        records=d.get('records')
        if records is None:
            records=[dict(station=x['job']['station'],status=x['result']['status'],result=x['result']) for x in d['results']]
        for record in records:
            if record.get('status')!='candidate_converged':continue
            result=record['result'];station=record['station'];point=result['point']
            assert result['verification']['J'][1][1]<0
            assert abs(point['growth'])<1e-7
            assert abs(result['verification_correction'][0])<1e-4
            assert abs(result['verification_correction'][1])<1e-5
            results[family,station]=dict(job=dict(source=family,station=station,name=name),
                result=result,source_path=str(path),source_sha256=sha(path))
    rows=[]
    for key,record in sorted(results.items()):
        point=record['result']['point'];v=record['result']['verification']
        rows.append(dict(source=key[0],station=key[1],x_mm=35*point['x'],omega=point['omega_r'],
            beta=point['beta_r'],alpha_r=point['alpha_r'],growth=point['growth'],
            h2=v['h'],g_omega_h2=v['J'][0][0],curvature_h2=v['J'][1][1],
            delta_omega_h2=record['result']['verification_correction'][0],
            delta_beta_h2=record['result']['verification_correction'][1],source_path=record['source_path']))
    summary=dict(created=stamp(),qualification='Locally continued interior beta-stationary neutral candidates; no global envelope, missing-branch stability or pole validity claim',
        provenance=provenance,results=list(results.values()),
        unresolved=[dict(source='M1p_signed',stations=list(range(16,0,-1)),reason='Local stationary-root corrector unresolved; finite ridge windows have growing edge samples'),
                    dict(source='M1c',stations=list(range(10,0,-1)),reason='Local stationary-root corrector unresolved after actual iloc transport; finite ridge windows have growing edge samples')],
        counts=dict(accepted=len(rows),new_since_1616=len(rows)-1),
        evidence=dict(M1p='nose_M1p_adaptive_ridge_20261023.json',M1c='nose_M1c_adaptive_ridge_20261023.json'))
    target=ROOT/'outputs/nose_local_candidates_20261023_1616.json'
    if target.exists():raise RuntimeError('Existing export preserved')
    target.write_text(json.dumps(summary,indent=2)+'\n')
    with (ROOT/'outputs/nose_local_candidates_20261023_1616.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    snapshot();print(json.dumps(dict(counts=summary['counts'],endpoints=[v for v in rows if (v['source']=='M1p_signed' and v['station']==17) or (v['source']=='M1c' and v['station']==11)])),flush=True)

if __name__=='__main__':main()
