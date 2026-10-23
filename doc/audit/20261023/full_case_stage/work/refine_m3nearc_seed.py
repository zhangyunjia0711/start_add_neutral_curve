#!/usr/bin/env python3
"""Refine the actual nose M3near conjugate peak and export trace-compatible input."""
import datetime,json
from pathlib import Path
from full_driver import ROOT,sha,snapshot
from fd_envelope import refine

def main():
    snapshot()
    source=ROOT/'outputs/nose119_M3nearc_neutral_initial_20261023.json'
    data=json.loads(source.read_text());point=data['point']
    solver=ROOT/'outputs/bin/START_LST_wall_local_v7'
    expected='7aa6016570f406587156f229107c0e5e4f76c8a85ba89dda267039ff1d688d79'
    assert sha(solver)==expected
    result=refine(119,point,point['omega_r'],point['beta_r'],solver,
        coordinates='wall_local',label='M3nearc_i119_actual_neutral_20261023',
        omega_bounds=(0,80),beta_bounds=(-80,80))
    output=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),
        original_trace_path=str(source),original_trace_sha256=sha(source),original_trace_preserved=True,
        solver=str(solver),solver_sha256=expected,
        results=[dict(job=dict(station=119,source='M3nearc',name='M3nearc_i119_actual_neutral_v7'),result=result)],
        qualification='Single continued local beta-stationary neutral candidate; not a complete beta/family envelope')
    target=ROOT/'outputs/M3nearc_i119_seed_results_20261023.json'
    target.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(path=str(target),status=result['status'],point=result.get('point'),reason=result.get('reason'))),flush=True)
    snapshot()

if __name__=='__main__':main()
