"""Create isolated START inputs from the supplied dimensional rough-case files."""
import configparser,hashlib,json,shutil,sys
from pathlib import Path
import h5py,numpy as np
from inspect_plot3d import CASE,OLD,records

TASK=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/linear_curvature/theta1_0_f_1_over_2/rough/neutral_M3_beta_minus7_20261024')
PREV=OLD.parent.parent/'neutral_M3_beta_minus7_closed_20261024'
FULL=OLD.parent.parent/'neutral_beta_full_20261022'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    assert not TASK.exists(),'Preserve existing task'
    for d in ['work','outputs','runs','inputs','M3/in','reference/outputs','matrix_exports']:(TASK/d).mkdir(parents=True,exist_ok=True)
    g=np.frombuffer(records(CASE/'fpp_output/grid_vertical2d.grd')[-1],dtype='<f8').reshape((2,400,361)).transpose(1,2,0).copy()
    q=np.frombuffer(records(CASE/'fpp_output/flow_vertical2d.dat')[-1],dtype='<f8').reshape((6,400,361)).transpose(1,2,0).copy()
    ref=json.loads((CASE/'fpp_output/summary.json').read_text())['reference']
    scales=np.array([ref['density'],ref['velocity'],ref['velocity'],ref['velocity'],ref['temperature'],ref['pressure']])
    gn=g/.035;qn=q/scales[None,None,:]
    for name,key,data,kind in [('grid_vertical2d_N.h5','grid',gn,'grid'),('flow_vertical2d_N.h5','baseflow',qn,'flow')]:
        split=np.zeros(data.shape+(2,));split[:,:,:,0]=data
        with h5py.File(OLD/name) as old:attrs=dict(old[key].attrs)
        path=TASK/'M3/in'/name
        with h5py.File(path,'w') as f:
            ds=f.create_dataset(key,data=split,compression='gzip',compression_opts=1)
            for k,v in attrs.items():ds.attrs[k]=v
            ds.attrs['plot3d_shape']=np.array([361,400],dtype=np.int64)
            ds.attrs['fortran_logical_shape']=np.array([2,data.shape[2],361,400],dtype=np.int64)
        with h5py.File(path) as f:assert np.array_equal(f[key][:,:,:,0],data) and np.count_nonzero(f[key][:,:,:,1])==0
    start=int(np.argmin(abs(gn[0,:,0]*35-(-25.07056802))))+1
    cfg=configparser.ConfigParser();cfg.read(OLD.parent/'config.ini')
    cfg['domain'].update({'in':'361','jn':'400','istart':'1','iend':'361','iloc':str(start)})
    cfg['filenames']['dir']=str(TASK/'M3/in')+'/'
    cfg['lst']['profile_coordinates']='wall_local'
    cfg['neutral']={'enabled':'F'}
    with (TASK/'M3/config.ini').open('w') as f:cfg.write(f)
    oldstart=json.loads((PREV/'outputs/native_points/i119_p0.json').read_text())
    state={k:oldstart['point'][k] for k in ['omega_r','beta_r','alpha_r']};state['alpha_i']=0.
    hints=[dict(iloc=i,**state) for i in range(1,362)]
    (TASK/'outputs/prior_axis_seeds.json').write_text(json.dumps(dict(records=hints,qualification='These are only fixed coefficient-export targets, not converged rough-case eigenvalues.'),indent=2)+'\n')
    (TASK/'reference/outputs/seed_catalog.json').write_text(json.dumps(dict(grids={'main_841':{'x_nondim':gn[0,:,0].tolist()}}))+'\n')
    original_names=['plot3dx.unf','R.plot3dq.unf','U.plot3dq.unf','V.plot3dq.unf','W.plot3dq.unf','T.plot3dq.unf','P.plot3dq.unf',
        'fpp_output/grid_vertical2d.grd','fpp_output/flow_vertical2d.dat','fpp_output/flow_vertical2d.variables.json','fpp_output/summary.json','fpp_output/fpp_ui_config.yaml']
    files=[dict(path=str(CASE/n),sha256=sha(CASE/n)) for n in original_names]
    files += [dict(path=str(p),sha256=sha(p)) for p in sorted((TASK/'M3').rglob('*')) if p.is_file()]
    (TASK/'outputs/input_manifest.json').write_text(json.dumps(files,indent=2)+'\n')
    metadata=dict(source_case=str(CASE),reference=ref,Lref_m=.035,primitive_scales=scales.tolist(),nominal_Ma=6.,nominal_Re=79100.,
        dimensional_reference_Re=ref['density']*ref['velocity']*.035/ref['viscosity'],
        nominal_Re_note='Use the same nominal Ma=6/Re=79100/Te=226.5 solver configuration as the previous case; dimensional FPP auto-reference gives Re=79091.941226, a 0.0102 percent rounding/model difference.',
        shape=[361,400],start_station=start,start_x_mm=float(g[0,start-1,0]*1000),
        input_grid_sha256=sha(CASE/'fpp_output/grid_vertical2d.grd'),input_flow_sha256=sha(CASE/'fpp_output/flow_vertical2d.dat'),
        H5_roundtrip_exact=True,model='wall-local; local parallel; zero explicit curvature',old_seed_source=str(PREV/'outputs/native_points/i119_p0.json'))
    (TASK/'outputs/input_conversion.json').write_text(json.dumps(metadata,indent=2,allow_nan=False)+'\n')
    shutil.copy2(FULL/'work/full_driver.py',TASK/'work/run_adapter.py')
    common='''"""Isolated rough-case run adapter with immutable source guards."""
import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
PREVIOUS=ROOT/'reference'
import run_adapter as driver
driver.ROOT=ROOT;driver.CASE=ROOT
run,table,sha=driver.run,driver.table,driver.sha
SOLVER=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/bluntplate/neutral_beta_full_20261022/outputs/bin/START_LST_wall_local_v10')
SOLVER_HASH='634b45c5d30e076f5b86ccb136e0734c1905beedde108477dceb8a75adc3853c'
OPTIONS=['-pep_tol','1e-12','-pep_max_it','1000','-pep_nev','1','-pep_ncv','8','-st_pc_factor_mat_solver_type','superlu']
assert sha(SOLVER)==SOLVER_HASH
def snapshot():
    result=json.loads((ROOT/'outputs/input_manifest.json').read_text())
    for e in result:assert sha(e['path'])==e['sha256'],e['path']
    assert sha(SOLVER)==SOLVER_HASH
    return result
'''
    (TASK/'work/common.py').write_text(common)
    for name in ['fixed_beta_pencil.py','profile_interpolation.py','prepare_native_pencils.py','validate_closed_contour.py','native_neutral_points.py','physical_turning_points.py']:
        code=(PREV/'work'/name).read_text()
        if name=='fixed_beta_pencil.py':
            code=code.replace("template=ROOT.parent/str(mode)","template=ROOT/str(mode)").replace('xi<=841','xi<=361').replace(')),840)', ')),360)')
        if name=='prepare_native_pencils.py':code=code.replace('default=735','default=361')
        if name=='profile_interpolation.py':
            code=code.replace('from trace_fixed_neutral import field','')
            code=code.replace("SOURCE=ROOT.parent/'M3/in'","SOURCE=ROOT/'M3/in'")
            a=code.index("            base=AXIS[i];rows,name=field(")
            b=code.index("            with h5py.File(SOURCE/",a)
            code=code[:a]+'''            from fixed_beta_pencil import OperatorLinearFamily
            native=OperatorLinearFamily().native(i)
            directory=Path(native.exports[0]['files'][0]['path']).parent.parent
            with (directory/'wall_local_profile_metadata.csv').open() as f:
                meta=next(q for q in csv.DictReader(line for line in f if not line.startswith('#')) if int(q['iloc'])==i)
'''+code[b:]
            code=code.replace("cfg.read(ROOT.parent/'M3/config.ini')","cfg.read(ROOT/'M3/config.ini')")
        if name=='native_neutral_points.py':code=code.replace("ROOT.parent/'neutral_beta_full_20261022/outputs/seed_catalog.json'","ROOT/'reference/outputs/seed_catalog.json'")
        (TASK/'work'/name).write_text(code)
    shutil.copy2(CASE/'work/inspect_plot3d.py',TASK/'work/inspect_plot3d_source.py')
    shutil.copy2(Path(__file__),TASK/'work/bootstrap_case.py')
    shutil.copy2(CASE/'outputs/input_inspection.json',TASK/'outputs/input_inspection.json')
    print(json.dumps(metadata,indent=2))

if __name__=='__main__':main()
