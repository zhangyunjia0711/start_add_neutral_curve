#!/usr/bin/env python3
"""Measure PEP-tolerance sensitivity of finite differences at a stalled point."""
import json
from full_driver import ROOT,sha,snapshot
from fd_envelope import stencil
from local_peak_x_continue import stamp

def main():
    src=ROOT/'outputs/local_peak_nose_M1c_fine_v10_1616_progress.json'
    d=json.loads(src.read_text());point=d['failures'][-1]['result']['trace'][-1]['center']
    solver=ROOT/'outputs/bin/START_LST_wall_local_v10';out=ROOT/'outputs/nose_M1c_i10_fd_tolerance_probe_20261023.json'
    if out.exists():raise RuntimeError('Existing evidence preserved')
    report=dict(created=stamp(),source=str(src),source_sha256=sha(src),actual_center=point,
        qualification='Empirical tolerance and stencil-size sensitivity; not a rigorous error bound or a neutral root',stencils=[],comparisons=[])
    snapshot()
    for h in [1e-6,2.5e-7]:
        pair=[]
        for tight in [False,True]:
            r=stencil(10,point,point['omega_r'],point['beta_r'],solver,'wall_local',h=h,tight=tight,
                label='nose_M1c_i10_noise_'+('tight' if tight else 'normal'))
            report['stencils'].append(r);pair.append(r)
            out.write_text(json.dumps(report,indent=2)+'\n')
        diffs=[abs(complex(a['alpha_r']-b['alpha_r'],a['alpha_i']-b['alpha_i'])) for a,b in zip(pair[0]['samples'],pair[1]['samples'])]
        report['comparisons'].append(dict(h=h,max_alpha_difference=max(diffs),
            center_alpha_difference=abs(complex(pair[0]['center']['alpha_r']-pair[1]['center']['alpha_r'],pair[0]['center']['alpha_i']-pair[1]['center']['alpha_i'])),
            normal_F=pair[0]['F'],tight_F=pair[1]['F'],normal_J=pair[0]['J'],tight_J=pair[1]['J']))
    report['finished']=stamp();out.write_text(json.dumps(report,indent=2)+'\n');snapshot();print(json.dumps(report['comparisons']),flush=True)

if __name__=='__main__':main()
