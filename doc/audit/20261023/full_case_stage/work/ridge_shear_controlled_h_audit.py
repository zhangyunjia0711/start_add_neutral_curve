#!/usr/bin/env python3
"""Separate same-tolerance h/h2 differences from the tight-PEP comparison."""
import json
from pathlib import Path
from types import SimpleNamespace
from full_driver import ROOT,sha,snapshot
from local_peak_x_continue import stamp
from ridge_shear_corrector import Sampler,accepted_root

def main():
    source_path=ROOT/'outputs/ridge_shear_M1p28_20261023.json';prior=json.loads(source_path.read_text())
    source=json.loads(Path(prior['source_report']).read_text());coarse=prior['validation']['h'];tight=prior['validation']['h2']
    args=SimpleNamespace(label='ridge_shear_controlled_h_M1p28',path_step=5e-6)
    sampler=Sampler(args,source);point=coarse['center'];out=ROOT/'outputs/ridge_shear_controlled_h_M1p28_20261023.json'
    if out.exists():raise RuntimeError('Existing audit preserved')
    snapshot();fine=sampler.stencil(28,point,point['omega_r'],point['beta_r'],coarse['k'],
        coarse['h_omega']/2,coarse['h_eta']/2,'same_tol_h2',False)
    rel=lambda a,b:abs(a-b)/max(abs(b),1e-30)
    report=dict(created=stamp(),source=str(source_path),source_sha256=sha(source_path),auxiliary_only=True,
        fine_same_PEP_tolerance=fine,
        same_tolerance_h_vs_h2=dict(g_u_relative=rel(coarse['J_shear'][0][0],fine['J_shear'][0][0]),
            curvature_relative=rel(coarse['J_shear'][1][1],fine['J_shear'][1][1]),
            original_and_fine_PEP_tol='1e-13',both_roots_pass_engineering_checks=accepted_root(coarse) and accepted_root(fine)),
        fixed_h2_PEP_tol_1e13_vs_1e14=dict(g_u_relative=rel(fine['J_shear'][0][0],tight['J_shear'][0][0]),
            curvature_relative=rel(fine['J_shear'][1][1],tight['J_shear'][1][1])),run_names=sampler.runs,
        qualification='Numerical FD/PEP comparison at an auxiliary accepted root, not a physical endpoint error bound')
    out.write_text(json.dumps(report,indent=2)+'\n');snapshot();print(json.dumps({k:v for k,v in report.items() if k not in ['fine_same_PEP_tolerance','run_names']}),flush=True)

if __name__=='__main__':main()
