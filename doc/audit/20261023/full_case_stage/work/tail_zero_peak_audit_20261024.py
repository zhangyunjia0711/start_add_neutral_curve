#!/usr/bin/env python3
"""Finite insert-frequency beta scans audit the large last peak displacement."""
import json
from types import SimpleNamespace
from full_driver import ROOT,sha,snapshot
from local_peak_x_continue import stamp
from tail_lower_ridge_841_fine import recovery_scan,correct_peak

def main():
    snapshot();source=ROOT/'outputs/tail_lower_beta6_841_fine_20261024.json'
    doc=json.loads(source.read_text());assert doc['finished'] and doc['reached_zero']
    origin=min((r['point'] for r in doc['rows'] if r['point']['omega_r']>0),key=lambda p:p['omega_r'])
    args=SimpleNamespace(label='tail841_zero_audit_20261024',solver=ROOT/'outputs/bin/START_LST_wall_local_v10',eigen_jump=4.)
    args.solver_sha=sha(args.solver)
    out=ROOT/'outputs/tail841_zero_peak_audit_20261024.json'
    if out.exists():raise RuntimeError('Existing finite zero audit preserved')
    report=dict(created=stamp(),source_path=str(source),source_sha256=sha(source),
        actual_origin=origin,finite_beta_domain=[-.1,1.5],omega_insertions=[.04,.03,.02,.01],rows=[],
        qualification='Finite beta-window scans and actual local maxima; do not assume an old peak survives when several maxima or edge ascent occur',finished=False)
    for index,w in enumerate(report['omega_insertions']):
        try:
            _,evidence=recovery_scan(args,origin,w,origin['beta_r'],[-.1,1.5],'insert'+str(index))
            row=dict(omega=w,scan=evidence,refined_peaks=[])
            ordered=evidence['samples']
            for j,point in enumerate(ordered):
                if point not in evidence['all_observed_peaks']:continue
                window=[ordered[j-1]['beta_r'],ordered[j+1]['beta_r']]
                peak=correct_peak(args,point,w,point['beta_r'],window,'insert'+str(index)+'_peak'+str(j))
                row['refined_peaks'].append(peak)
            report['rows'].append(row)
        except Exception as error:
            report['rows'].append(dict(omega=w,status='unknown',reason=repr(error)))
        out.write_text(json.dumps(report,indent=2)+'\n')
    report.update(finished=True,original_inputs_unchanged=bool(snapshot()))
    out.write_text(json.dumps(report,indent=2)+'\n')
    print('ZERO AUDIT '+json.dumps([dict(omega=r['omega'],peaks=[dict(status=p['status'],point=p.get('point'),curvature=p.get('curvature'))
        for p in r.get('refined_peaks',[])],reason=r.get('reason')) for r in report['rows']]),flush=True)

if __name__=='__main__':main()
