#!/usr/bin/env python3
"""Finite PEP-tolerance audit of the final accepted auxiliary narrow peak."""
import json
from pathlib import Path
import full_driver as driver
import fd_envelope as fd
from local_peak_x_continue import stamp

def main():
    root=driver.ROOT;source=root/'outputs/aux_neutral_M1p_17_16_s64_20261023.json'
    d=json.loads(source.read_text());record=d['records'][-1];point=record['point'];column=record['auxiliary_column']
    target=root/'outputs/aux_M1p_c28_fd_tolerance_probe_20261023.json'
    if target.exists():raise RuntimeError('Existing audit preserved')
    original_run=driver.run
    def auxiliary_run(spec):
        spec.update(mode=d['template'],auxiliary_input_hashes=d['input_hashes'],station_predictor=False)
        return original_run(spec)
    fd.run=auxiliary_run
    report=dict(created=stamp(),source=str(source),source_sha256=driver.sha(source),auxiliary_column=column,
        theta=record['theta'],actual_point=point,stencils=[],comparisons=[],auxiliary_only=True,
        qualification='Empirical tolerance sensitivity at an auxiliary profile; not a rigorous error bound or physical endpoint')
    driver.snapshot()
    for h in [1.25e-7,6.25e-8]:
        pair=[]
        for tight in [False,True]:
            result=fd.stencil(column,point,point['omega_r'],point['beta_r'],Path(d['solver']),h=h,tight=tight,
                label='aux_M1p_c28_noise_'+('tight' if tight else 'normal'))
            pair.append(result);report['stencils'].append(result)
            target.write_text(json.dumps(report,indent=2)+'\n')
        delta=[abs(complex(a['alpha_r']-b['alpha_r'],a['alpha_i']-b['alpha_i'])) for a,b in zip(pair[0]['samples'],pair[1]['samples'])]
        report['comparisons'].append(dict(h=h,max_alpha_difference=max(delta),
            normal_F=pair[0]['F'],tight_F=pair[1]['F'],normal_J=pair[0]['J'],tight_J=pair[1]['J']))
    report['finished']=stamp();target.write_text(json.dumps(report,indent=2)+'\n');driver.snapshot()
    print(json.dumps(report['comparisons']),flush=True)

if __name__=='__main__':main()
