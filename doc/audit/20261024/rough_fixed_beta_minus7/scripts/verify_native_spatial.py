"""Independent original spatial-PEP checks of rough-case native neutral points."""
import concurrent.futures,datetime,json
from pathlib import Path
import numpy as np
from common import ROOT,sha,snapshot
from profile_interpolation import evaluate,mac

def main():
    manifest=ROOT/'outputs/native_neutral_points.json';d=json.loads(manifest.read_text());assert d['finished']
    stations={r['station'] for r in d['results']};selected={3,4,8,16,17,24,50,60,61,100,133,140,147,200,226,235,242,260,280,294}&stations
    selected|={min(stations),max(stations)}
    jobs=[r for r in d['results'] if r['station'] in selected and r['status']=='converged']
    report=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),native_points_manifest_sha256=sha(manifest),results=[],finished=False)
    out=ROOT/'outputs/native_spatial_verification.json'
    def one(entry):
        source=json.loads(Path(entry['file']).read_text());point=source['point'];station=entry['station']
        r=evaluate('M3',station,point,point['omega_r'],1e-7,f'rough_native_fresh_i{station:03d}_p{entry["ordinal"]}')
        q=np.load(source['vector_file']);delta=abs(complex(r['point']['alpha_r']-point['alpha_r'],r['point']['alpha_i']));overlap=mac(q,r['q'])
        answer=dict(station=station,ordinal=entry['ordinal'],native_point_file=entry['file'],actual_spatial_point=r['point'],
            actual_polynomial_residual=r['polynomial_residual'],spatial_run=r['run'],alpha_difference=delta,raw_MAC=overlap,
            vector_file=r['vector_file'],vector_index=r['vector_index'],side_samples=r['samples'],
            side_growth_sign_crossing=r['samples'][0]['growth']*r['samples'][-1]['growth']<0)
        answer['passed']=delta<1e-5 and abs(r['point']['growth'])<1e-6 and overlap>1-1e-5
        return answer
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        for r in pool.map(one,jobs):
            report['results'].append(r);out.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
            print('SPATIAL_VERIFY',r['station'],r['passed'],r['alpha_difference'],r['raw_MAC'],r['side_growth_sign_crossing'],flush=True)
    report.update(finished=True,all_passed=all(r['passed'] for r in report['results']),original_inputs=snapshot())
    out.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')

if __name__=='__main__':main()
