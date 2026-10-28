"""Full-matrix compatibility and invariance checks on actual FPP profiles."""
import ast,configparser,json,sys
from pathlib import Path
import h5py,numpy as np
from scipy.sparse.linalg import norm
ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/neutral_m3_sweep_consistent_20261028')
sys.path.insert(0,str(ROOT/'work'))
import solver
from slice_input import template
def main():
    previous=ROOT.parent/'75°/neutral_m3_frequency_limits_20261027/work/solver.py'
    fn=next(n for n in ast.parse(previous.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='export')
    compiled=compile(ast.Module(body=[fn],type_ignores=[]),str(previous),'exec')
    scope=dict(solver.__dict__);exec(compiled,scope);newlegacy=scope['export']
    frozen=dict(scope);frozen['SOLVER']=ROOT.parent/'75°/bluntplate/neutral_beta_full_20261022/outputs/bin/START_LST_wall_local_v10';frozen['SOLVER_HASH']='634b45c5d30e076f5b86ccb136e0734c1905beedde108477dceb8a75adc3853c';exec(compiled,frozen);oldlegacy=frozen['export']
    rows=[]
    for case,i in [('s75',119),('s75',680),('s75',681),('s30',241)]:
        a,ra=newlegacy(case,i,-7.5,-6.1,label='new_legacy_regression');b,rb=oldlegacy(case,i,-7.5,-6.1,label='old_v10_regression')
        errors=[float(norm(x-y))/max(1.,float(norm(y))) for x,y in zip(a,b)];assert max(errors)<1e-12
        rows.append(dict(kind='legacy_compatibility',case=case,station=i,errors=errors,new=ra,old=rb));print('LEGACY',case,i,errors,flush=True)
    source,local,receipt=template('s75',681);dest=ROOT/'inputs/consistent_roundoff_regression';(dest/'in').mkdir(parents=True)
    for name,key in [('grid_vertical2d_N.h5','grid'),('flow_vertical2d_N.h5','baseflow')]:
        with h5py.File(source/'in'/name) as h:a=h[key][...]
        if key=='grid':a[:,:,0,0]+=1e-13*(a[:,:,1,0]-a[0,:,1,0])
        with h5py.File(dest/'in'/name,'w') as h:h.create_dataset(key,data=a)
    cfg=configparser.ConfigParser();cfg.read(source/'config.ini')
    with (dest/'config.ini').open('w') as f:cfg.write(f)
    a,ra=solver.export(dest,local,-7.5,-6.1,label='consistent_roundoff');b=solver.Pencil('s75',681,-6.1).spatial(-7.5)
    errors=[float(norm(x-y))/max(1.,float(norm(y))) for x,y in zip(a,b)];assert max(errors)<1e-12
    rows.append(dict(kind='consistent_coordinate_roundoff_invariance',case='s75',station=681,errors=errors,export=ra));print('INVARIANCE',errors,flush=True)
    # Full-field and cropped-native-stencil comparison must use the new mode on both sides.
    text=previous.read_text().replace("'wall_local'","'wall_local_consistent'");fn=next(n for n in ast.parse(text).body if isinstance(n,ast.FunctionDef) and n.name=='export');scope=dict(solver.__dict__);exec(compile(ast.Module(body=[fn],type_ignores=[]),'consistent_full_export','exec'),scope)
    for case,i in [('s30',119),('s45',241),('s60',400),('s60',2),('s75',681)]:
        a,ra=solver.export(case,i,-7.,2.3,label='consistent_stencil_validation');b,rb=scope['export'](case,i,-7.,2.3,label='consistent_full_validation')
        errors=[float(norm(x-y))/max(1.,float(norm(y))) for x,y in zip(a,b)];assert max(errors)<1e-11
        rows.append(dict(kind='consistent_native_stencil',case=case,station=i,errors=errors,cropped=ra,full=rb));print('STENCIL',case,i,errors,flush=True)
    report=dict(passed=True,checks=rows,solver_sha256=solver.SOLVER_HASH)
    (ROOT/'outputs/stencil_validation.json').write_text(json.dumps(report,indent=2)+'\n');(ROOT/'outputs/consistent_validation.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
