#!/usr/bin/env python3
"""Compare independently computed temporal-pencil and spatial-PEP primitive vectors."""
import json,sys
from pathlib import Path
import numpy as np
from full_driver import ROOT,sha,snapshot
from local_peak_x_continue import stamp

def main():
    sys.path.insert(0,str(ROOT/'work/pep_matrix_temporal_probe_20261023'))
    import temporal_pencil as tp
    source=ROOT/'outputs/native16_stationary_beta_minimum_20261023.json';minimum=json.loads(source.read_text())
    matrix_path=Path(minimum['source_path']);matrix=json.loads(matrix_path.read_text())
    assert sha(Path(tp.__file__))==matrix['helper_sha256']
    anchor=next(x for x in json.loads(Path(matrix['source_path']).read_text())['cases'] if x['family']=='M1p')['actual_native']
    init=json.loads(Path(matrix['initialization_path']).read_text())['predicted_native16']
    p=tp.Pencil(ROOT/'work/matrix_neutral_M1p16_20261023/native16_matrices',anchor,16,init['beta_r'],init['omega_r'])
    p.build();native=minimum['actual_native_fresh']
    q=tp.vf.vectors(Path(native['native_vector_file']),400)[native['raw_index']].ravel()
    c=matrix['native16_correction']['point'];r,v,_=p.evaluate(c['alpha'],c['beta'],q,c['omega']['real']+1e-3,left=True)
    overlap=tp.mac(q,v);omega_error=abs(complex(r['omega']['real']-c['omega']['real'],r['omega']['imag']))
    passed=bool(minimum['passed'] and overlap>=1-1e-8 and omega_error<1e-8)
    out=ROOT/'outputs/native16_stationary_beta_minimum_raw_vector_20261023.json'
    if out.exists():raise RuntimeError('Existing vector comparison preserved')
    report=dict(created=stamp(),source_path=str(source),source_sha256=sha(source),matrix_source_path=str(matrix_path),
        matrix_source_sha256=sha(matrix_path),helper_sha256=sha(Path(tp.__file__)),
        matrix_evaluation=r,native_spatial_raw_vector=native['native_vector_file'],
        independent_temporal_spatial_raw_primitive_MAC=overlap,temporal_omega_difference=omega_error,passed=passed,
        envelope_candidate=False,qualification='Same native point, independently assembled spatial PEP and reconstructed temporal pencil; beta minimum remains excluded from maximum-growth envelope')
    out.write_text(json.dumps(report,indent=2)+'\n');snapshot();print(json.dumps(dict(passed=passed,raw_MAC=overlap,omega_error=omega_error)),flush=True)

if __name__=='__main__':main()
