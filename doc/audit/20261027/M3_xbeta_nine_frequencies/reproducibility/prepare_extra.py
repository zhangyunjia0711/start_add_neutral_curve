"""Create an isolated run using proven contour/fold routines and read-only caches."""
import datetime,hashlib,json,shutil
from pathlib import Path
BASE=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°')
SOURCE=BASE/'neutral_xbeta_four_omega_20261026'
ROOT=BASE/'neutral_xbeta_extra_frequencies_20261027'
FS=[-75,-25,25,75]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def replace(s,a,b):
    assert s.count(a)==1,(a[:80],s.count(a))
    return s.replace(a,b)
def main():
    assert not ROOT.exists(),'Preserve any prior attempt'
    for d in ['work','outputs','runs','matrix_exports','inputs']:(ROOT/d).mkdir(parents=True,exist_ok=True)
    for d in ['native_pencils','cases']:(ROOT/d).symlink_to(SOURCE/d,target_is_directory=True)
    names=['run_adapter.py','core.py','solver.py','contour.py','contour3d.py','physical.py','tip_correction.py','finalize.py']
    provenance=[]
    for name in names:
        source=SOURCE/'work'/name;s=source.read_text()
        provenance.append(dict(name=name,source=str(source),sha256=sha(source)))
        if name in ['contour3d.py','finalize.py']:s=replace(s,'choices=[50,100]','choices=[-75,-25,25,75]');s=replace(s,"choices=['A','B']","choices=['S']")
        if name=='contour3d.py':
            s=replace(s,'step=.025\n    state=','step=.025;half=0\n    state=')
            s=replace(s,'for loop in range(4000):','for loop in range(8000):')
            s=replace(s,"accepted=False\n        for retry", "accepted=False;edge=False\n        for retry")
            s=replace(s,"assert 1<=boundary<=family.count,'Neutral curve reaches domain edge'",'if not 1<=boundary<=family.count:\n                    edge=True;break')
            s=replace(s,"        if not accepted:state.update(status='unresolved_minimum_step',finished=True);break",'''        if edge:
            state.setdefault('domain_endpoints',[]).append(dict(index=len(state['records'])-1,point=r))
            if half==0:
                half=1;z=z0.copy();q=q0.copy();p=p0.copy();r=r0.copy();t=-t0.copy();step=.025
                state['reverse_half_start_index']=len(state['records']);store(z,r,q,p,t,step,0);continue
            state.update(status='open_at_domain_boundaries',finished=True);break
        if not accepted:state.update(status='unresolved_minimum_step',finished=True);break''')
            s=replace(s,"if len(state['records'])>30 and arc>.5", "if half==0 and len(state['records'])>30 and arc>.5")
        if name=='contour.py':
            start=s.index('def seed_root(');end=s.index('\ndef main():',start)
            s=s[:start]+'''def seed_root(case,frequency,component):
    source=ROOT/'outputs/neutral_seeds.json';data=json.loads(source.read_text())
    s=next(r for r in data['records'] if r['case']==case and r['f_kHz']==frequency and r['side']=='lower')
    q=np.load(s['vector_file']);i=s['station']
    fresh,v,receipt=original(case,i,s,s['beta_r'],f'{case}_f{frequency}kHz_S_neutral_seed')
    assert abs(fresh['growth'])<1e-8 and mac(v,q)>1-1e-7
    return dict(station=i,point=fresh,receipt=receipt,component=component,frequency_kHz=frequency,source_record=s),v
\n'''+s[end:]
        (ROOT/'work'/name).write_text(s)
    protected=json.loads((SOURCE/'outputs/input_manifest.json').read_text())
    protected += [dict(path=str(p),sha256=sha(p),role='prior_result') for p in sorted((SOURCE/'outputs').glob('*_final.json'))]
    protected += [dict(path=str(p),sha256=sha(p),role='prior_selection') for p in sorted((Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate/75°/outputs/M3_xbeta_frequency_overlay_20261027')).glob('*.csv'))]
    (ROOT/'outputs/input_manifest.json').write_text(json.dumps(protected,indent=2)+'\n')
    shutil.copy2(SOURCE/'outputs/cases.json',ROOT/'outputs/cases.json')
    conversion=json.loads((SOURCE/'outputs/frequency_conversion.json').read_text())
    conversion['frequencies']=[dict(f_kHz=f,omega_nd=6.0727867021443425*f/50) for f in FS]
    (ROOT/'outputs/frequency_conversion.json').write_text(json.dumps(conversion,indent=2)+'\n')
    (ROOT/'work/common.py').write_text('''"""Isolated additional-frequency context; original matrices and CFD stay unchanged."""
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
    shutil.copy2(Path(__file__),ROOT/'work/prepare_extra.py')
    (ROOT/'outputs/bootstrap_receipt.json').write_text(json.dumps(dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),source=str(SOURCE),root=str(ROOT),frequencies_kHz=FS,source_code=provenance,changes=['Signed new frequencies; same representative seeds from direct frequency continuation','Open-domain handling added to existing real-alpha/beta 3D pseudo-arclength guide','Original physical-profile turning point correction retained','Original caches shared for reading; new outputs/runs/physical-profile inputs isolated']),indent=2)+'\n')
    print('PREPARED',ROOT,len(protected),'protected inputs/results',flush=True)
if __name__=='__main__':main()
