#!/usr/bin/env python3
"""Export native profiles projected in the analytic radial basis, not solver FD basis."""
import csv,json
import h5py,numpy as np
from full_driver import ROOT,CASE,sha,snapshot
from local_peak_x_continue import stamp

def main():
    gp=CASE/'M1/in/grid_vertical2d_N.h5';qp=CASE/'M1/in/flow_vertical2d_N.h5'
    snapshot();rows=[]
    with h5py.File(gp,'r') as gf,h5py.File(qp,'r') as qf:
        for i in range(1,22):
            g=gf['grid'][:,i-1,:,0];q=qf['baseflow'][:,i-1,:,0]
            normal=g[0]/np.linalg.norm(g[0]);tangent=np.array([normal[1],-normal[0]])
            n=(g-g[0])@normal;us=q[:,1:3]@tangent;vn=q[:,1:3]@normal
            row=dict(station=i,x_N=float(g[0,0]),x_mm=float(35*g[0,0]),y_N=float(g[0,1]),
                tx=float(tangent[0]),ty=float(tangent[1]),nx=float(normal[0]),ny=float(normal[1]),
                Us_outer=float(us[-1]),Us_min=float(us.min()),Us_max=float(us.max()),Us_max_abs=float(abs(us).max()),
                Us_exactly_zero_all_nodes=bool(np.all(us==0)),Us_exact_zero_nodes=int(np.count_nonzero(us==0)),
                Vn_max_abs=float(abs(vn).max()),W_outer=float(q[-1,3]),n_max_N=float(n[-1]),
                minimum_delta_n=float(np.diff(n).min()),normal_grid_strictly_increasing=bool(np.all(np.diff(n)>0)))
            rows.append(row)
    report=dict(created=stamp(),source_grid=str(gp),source_grid_sha256=sha(gp),
        source_flow=str(qp),source_flow_sha256=sha(qp),projection='native u_s=u*t_x+v*t_y, v_n=u*n_x+v*n_y',
        basis_definition='Analytic radial normal g_wall/|g_wall|, tangent=(n_y,-n_x). Not the solver finite-difference frozen basis.',
        model='Analytic radial projection audit; exact frozen solver basis/Us must be read from separate native metadata',
        qualification='Radial-projection facts only. A zero radial Us at i1 is not a claim of exactly zero Us in the native FD-basis operator, and establishes neither model validity nor root existence.',rows=rows)
    p=ROOT/'outputs/nose_native_profile_metadata_20261023.json'
    if p.exists():raise RuntimeError('Preserve prior metadata')
    p.write_text(json.dumps(report,indent=2)+'\n')
    with (ROOT/'outputs/nose_native_profile_metadata_20261023.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    snapshot();print(json.dumps([v for v in rows if v['station'] in [1,2,10,11,16,17]]),flush=True)

if __name__=='__main__':main()
