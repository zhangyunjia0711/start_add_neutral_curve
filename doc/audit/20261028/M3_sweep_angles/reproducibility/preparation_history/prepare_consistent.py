"""Reuse immutable FPP fields, with a separately built consistent derivative mode."""
import configparser,hashlib,json,re,shutil,sys
from pathlib import Path
import numpy as np
BASE=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate')
OLD=BASE/'neutral_m3_sweep_angles_fpp_20261028'
ROOT=BASE/'neutral_m3_sweep_consistent_20261028'
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    assert not (ROOT/'outputs/cases.json').exists(),'Preserve completed preparation'
    build=json.loads((ROOT/'outputs/build_receipt.json').read_text());assert all(r['returncode']==0 for r in build['commands'])
    for p in (OLD/'work').glob('*.py'):
        s=p.read_text().replace('neutral_m3_sweep_angles_fpp_20261028','neutral_m3_sweep_consistent_20261028')
        s=s.replace("'wall_local'","'wall_local_consistent'")
        if p.name=='common.py':
            s=re.sub(r'^SOLVER=Path\(.*$',"SOLVER=Path("+repr(build['solver'])+")",s,flags=re.M)
            s=re.sub(r"^SOLVER_HASH=.*$","SOLVER_HASH="+repr(build['solver_sha256']),s,flags=re.M)
        (ROOT/'work'/p.name).write_text(s)
    for name in ['frequency_seeds_fpp.py','run_case_fpp.py']:
        s=(HERE/name).read_text().replace('neutral_m3_sweep_angles_fpp_20261028','neutral_m3_sweep_consistent_20261028');(ROOT/'work'/name).write_text(s)
    cases=json.loads((OLD/'outputs/cases.json').read_text());protected=json.loads((OLD/'outputs/input_manifest.json').read_text())
    for case,c in cases.items():
        dest=ROOT/'cases'/case;dest.mkdir();(dest/'in').symlink_to(OLD/'cases'/case/'in',target_is_directory=True)
        cfg=configparser.ConfigParser();cfg.read(OLD/'cases'/case/'config.ini');cfg['filenames']['dir']=str(dest/'in')+'/';cfg['lst']['profile_coordinates']='wall_local_consistent'
        with (dest/'config.ini').open('w') as f:cfg.write(f)
        c['config']=str(dest/'config.ini');c['in_dir']=str(dest/'in');c['normal_derivative_mode']='wall_local_consistent'
        protected.append(dict(path=str(dest/'config.ini'),sha256=sha(dest/'config.ini')))
    protected.append(dict(path=build['solver'],sha256=build['solver_sha256']))
    (ROOT/'outputs/cases.json').write_text(json.dumps(cases,indent=2)+'\n');(ROOT/'outputs/input_manifest.json').write_text(json.dumps(protected,indent=2)+'\n')
    for name in ['input_selection.json','frequency_conversion.json']:shutil.copy2(OLD/'outputs'/name,ROOT/'outputs'/name)
    sys.path.insert(0,str(ROOT/'work'))
    from solver import Pencil
    from core import mac
    pilot=[]
    for case in cases:
        old=json.loads((OLD/'outputs'/f'{case}_neutral_seed_scan.json').read_text());records=[]
        for j,r in enumerate(old['records']):
            if r['f_kHz'] not in [-100.,-50.,0.,50.,100.] or r['station']!=119:continue
            p=Pencil(case,r['station'],r['omega_r']);ref=np.load(r['vector_file']);a,q,check=p.solve_alpha(r['beta_r'],complex(r['alpha_r'],r['alpha_i']),ref,True)
            assert abs(a.imag)<1e-8 and mac(ref,q)>1-1e-7
            vf=ROOT/'outputs'/f'{case}_verified_seed{j}.npy';np.save(vf,q);item=dict(r);item.update(vector_file=str(vf),alpha_r=a.real,alpha_i=a.imag,growth=-a.imag,consistent_mode_check=check,prior_seed_source=str(OLD/'outputs'/f'{case}_neutral_seed_scan.json'));records.append(item)
        (ROOT/'outputs'/f'{case}_neutral_seed_scan.json').write_text(json.dumps(dict(case=case,records=records,failures=[],scans=[],reuse='Every imported curved-station seed re-solved with consistent binary and checked by original polynomial residual and MAC.'),indent=2)+'\n')
        oldpilot=next(r for r in json.loads((OLD/'outputs/angle_mode_pilot.json').read_text()) if r['case']==case)
        p=Pencil(case,119,0.);a,q,check=p.solve_alpha(-7.,complex(oldpilot['alpha_r'],oldpilot['alpha_i']),np.load(oldpilot['vector_file']),True)
        vf=ROOT/'outputs'/f'{case}_mode_seed.npy';np.save(vf,q);item=dict(oldpilot);item.update(alpha_r=a.real,alpha_i=a.imag,growth=-a.imag,vector_file=str(vf),receipt=check);pilot.append(item)
        print('CONSISTENT_CASE',case,'verified_neutral_seeds',len(records),'mode',a,flush=True)
    (ROOT/'outputs/angle_mode_pilot.json').write_text(json.dumps(pilot,indent=2)+'\n');print('PREPARED',ROOT,flush=True)
if __name__=='__main__':main()
