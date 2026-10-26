"""Create an isolated two-case fixed-frequency calculation with immutable inputs."""
import configparser,hashlib,json,shutil
from pathlib import Path
import h5py,numpy as np

ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/neutral_xbeta_omega0_20261026')
BASE=ROOT.parent/'bluntplate'
OLD=BASE/'neutral_M3_beta_minus7_closed_20261024'
ROUGH=ROOT.parent/'linear_curvature/theta1_0_f_1_over_2/rough/neutral_M3_beta_minus7_20261024'
SOLVER=BASE/'neutral_beta_full_20261022/outputs/bin/START_LST_wall_local_v10'
HASH='634b45c5d30e076f5b86ccb136e0734c1905beedde108477dceb8a75adc3853c'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    assert not ROOT.exists(),'Preserve any existing calculation'
    assert sha(SOLVER)==HASH
    for d in ['work','outputs','runs','matrix_exports','native_pencils','inputs']:(ROOT/d).mkdir(parents=True)
    cases={};files=[]
    for name,source,station,state,flat in [('bluntplate',BASE/'M3',119,dict(alpha_r=35.93408168,alpha_i=-4.02320790),0.),
        ('rough',ROUGH/'M3',24,json.loads((ROUGH/'outputs/rough_physical_M3_seed.json').read_text())['point'],14.858163)]:
        dest=ROOT/'cases'/name;(dest/'in').mkdir(parents=True)
        cfg=configparser.ConfigParser();cfg.read(source/'config.ini');cfg['filenames']['dir']=str(dest/'in')+'/'
        cfg['lst']['omega']='0, 0';cfg['lst']['beta']='-7, 0';cfg['lst']['profile_coordinates']='wall_local'
        cfg['lst']['alpha']='%.17g, %.17g'%(state['alpha_r'],state['alpha_i']);cfg['neutral']={'enabled':'F'}
        cfg['domain']['iloc']=str(station)
        for item in ['grid_vertical2d_N.h5','flow_vertical2d_N.h5']:
            p=source/'in'/item;shutil.copy2(p,dest/'in'/item)
            assert sha(p)==sha(dest/'in'/item);files.append(dict(path=str(p),sha256=sha(p),role='original_input'))
        files.append(dict(path=str(source/'config.ini'),sha256=sha(source/'config.ini'),role='source_config'))
        with (dest/'config.ini').open('w') as f:cfg.write(f)
        with h5py.File(dest/'in/grid_vertical2d_N.h5') as f:g=f['grid'][:,:,:,0]
        assert g.shape[0]==400 and np.all(np.diff(g[0,:,0])>0)
        cases[name]=dict(name=name,config=str(dest/'config.ini'),in_dir=str(dest/'in'),source_dir=str(source),
            stations=g.shape[1],normal_points=400,x_mm=(g[0,:,0]*35).tolist(),wall_y_mm=(g[0,:,1]*35).tolist(),
            flat_start_mm=flat,seed_station=station,seed=dict(omega_r=0.,omega_i=0.,beta_r=-7.,beta_i=0.,alpha_r=state['alpha_r'],alpha_i=state['alpha_i']))
        for p in sorted(dest.rglob('*')):
            if p.is_file():files.append(dict(path=str(p),sha256=sha(p),role='isolated_input'))
    files.append(dict(path=str(SOLVER),sha256=HASH,role='frozen_solver'))
    (ROOT/'outputs/cases.json').write_text(json.dumps(cases,indent=2)+'\n')
    (ROOT/'outputs/input_manifest.json').write_text(json.dumps(files,indent=2)+'\n')
    shutil.copy2(ROUGH/'work/run_adapter.py',ROOT/'work/run_adapter.py')
    (ROOT/'work/common.py').write_text('''"""Immutable, isolated fixed-omega run context."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
import run_adapter as driver
driver.ROOT=ROOT;driver.CASE=ROOT/'cases'
run,table,sha=driver.run,driver.table,driver.sha
SOLVER=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/bluntplate/neutral_beta_full_20261022/outputs/bin/START_LST_wall_local_v10')
SOLVER_HASH='634b45c5d30e076f5b86ccb136e0734c1905beedde108477dceb8a75adc3853c'
CASES=json.loads((ROOT/'outputs/cases.json').read_text())
def snapshot():
    data=json.loads((ROOT/'outputs/input_manifest.json').read_text())
    for r in data:assert sha(r['path'])==r['sha256'],r['path']
    return data
assert sha(SOLVER)==SOLVER_HASH
''')
    shutil.copy2(Path(__file__),ROOT/'work/bootstrap.py')
    print(json.dumps(dict(root=str(ROOT),cases={k:dict(stations=v['stations'],seed_station=v['seed_station']) for k,v in cases.items()},protected_files=len(files))))

if __name__=='__main__':main()
