"""Audit native data coverage, immutable sources, syntax, and original CFD wall."""
import ast,csv,datetime,json,struct
from pathlib import Path
import numpy as np,h5py
from common import ROOT,sha,snapshot,SOLVER,SOLVER_HASH

def records(p,endian):
    data=p.read_bytes();offset=0;out=[]
    while offset<len(data):
        n=struct.unpack_from(endian+'I',data,offset)[0];offset+=4;assert offset+n+4<=len(data)
        out.append(data[offset:offset+n]);offset+=n;assert struct.unpack_from(endian+'I',data,offset)[0]==n;offset+=4
    return out

def main():
    O=ROOT/'outputs';stem='rough_M3_closed_beta_minus7';s=json.loads((O/(stem+'_summary.json')).read_text())
    info=json.loads((O/'input_conversion.json').read_text());case=Path(info['source_case'])
    with (O/(stem+'_curve.csv')).open() as f:rows=list(csv.DictReader(f))
    assert len(rows)==s['closed_curve_rows']==651
    assert [int(r['order']) for r in rows]==list(range(len(rows)))
    assert all(float(r['beta_r'])==-7 and float(r['alpha_i'])==0 for r in rows)
    assert all(rows[0][k]==rows[-1][k] for k in ['x_mm','omega_r','beta_r','alpha_r','alpha_i'])
    kinds={k:sum(int(r['point_kind'])==k for r in rows) for k in [0,1,2]};assert kinds=={0:644,1:6,2:1}
    assert len({(r['native_station'],r['alpha_r'],r['omega_r']) for r in rows if int(r['point_kind'])==0})==644
    stations={int(r['native_station']) for r in rows if int(r['point_kind'])==0};assert stations==set(range(3,295))
    assert all(Path(r['source_file']).is_file() for r in rows)
    assert s['native_spatial_side_sign_crossings']==s['native_spatial_independent_checks']==48
    for p in (ROOT/'work').glob('*.py'):ast.parse(p.read_text(),filename=str(p))
    rawg=records(case/'plot3dx.unf','>');dims=np.frombuffer(rawg[1],dtype='>i4');assert dims.tolist()==[361,301]
    g=np.frombuffer(rawg[-1],dtype='>f8').reshape((2,301,361)).transpose(1,2,0)
    with h5py.File(ROOT/'M3/in/grid_vertical2d_N.h5') as f:wall=f['grid'][0,:,:,0]*.035
    wall_error=float(np.max(abs(wall-g[0])));assert wall_error<1e-14
    rawwall=[]
    for name in ['R','U','V','W','T','P']:
        data=records(case/(name+'.plot3dq.unf'),'>');assert np.frombuffer(data[1],dtype='>i4').tolist()==[361,301,1]
        rawwall.append(np.frombuffer(data[-1],dtype='>f8').reshape((301,361))[0])
    rawwall=np.column_stack(rawwall)
    with h5py.File(ROOT/'M3/in/flow_vertical2d_N.h5') as f:qwall=f['baseflow'][0,:,:,0]*np.array(info['primitive_scales'])
    wall_flow_error=np.max(abs(rawwall-qwall),axis=0);assert np.max(wall_flow_error)<1e-8
    preprocessing=json.loads((case/'fpp_output/summary.json').read_text())['interpolation']['blocks'][0]
    report=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),rows=651,point_kind_counts=kinds,
        all_beta_minus7=True,first_last_physical_coordinates_exact=True,no_operator_guide_nodes_exported=True,
        native_station_count=292,native_station_range=[3,294],native_neutral_points=644,
        source_files_unchanged=True,original_inputs=snapshot(),solver_sha256=sha(SOLVER),python_syntax_passed=True,
        raw_CFD_shape=[361,301],supplied_FPP_shape=[361,400],raw_CFD_vs_FPP_wall_max_coordinate_error_m=wall_error,
        raw_CFD_vs_FPP_wall_max_primitive_errors=wall_flow_error.tolist(),
        existing_FPP_preprocessing_diagnostics={k:preprocessing[k] for k in ['inverse_map_unconverged_points','inverse_map_extrapolated_points','inverse_map_max_residual','bound_limited_values']},
        model_limit='The supplied FPP interpolated baseflow is used. Its existing interpolation diagnostics are retained; this calculation is not an independent CFD remapping or mesh-convergence study.')
    assert report['solver_sha256']==SOLVER_HASH
    (O/(stem+'_audit.json')).write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps(report))

if __name__=='__main__':main()
