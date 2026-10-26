"""Isolated nonzero-frequency task; prior omega=0 calculations stay immutable."""
import hashlib,json,shutil,datetime
from pathlib import Path
BASE=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°')
PREVIOUS=BASE/'neutral_xbeta_omega0_20261026'
ROOT=BASE/'neutral_xbeta_four_omega_20261026'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    assert not ROOT.exists(),'Keep existing results'
    prior=json.loads((PREVIOUS/'outputs/input_manifest.json').read_text())
    assert all(sha(p['path'])==p['sha256'] for p in prior)
    for name in ['work','outputs','runs','matrix_exports','native_pencils','inputs']:(ROOT/name).mkdir(parents=True)
    shutil.copytree(PREVIOUS/'cases',ROOT/'cases')
    cases=json.loads((PREVIOUS/'outputs/cases.json').read_text())
    for case,c in cases.items():
        config=ROOT/'cases'/case/'config.ini'
        config.write_text(config.read_text().replace(str(PREVIOUS),str(ROOT)))
        c['config']=str(config);c['in_dir']=str(config.parent/'in')
    protected=[dict(path=str(p),sha256=sha(p),role='isolated_input') for p in sorted((ROOT/'cases').rglob('*')) if p.is_file()]
    protected += [dict(r,role='prior_protected_input') for r in prior]
    protected += [dict(path=str(PREVIOUS/'outputs'/name),sha256=sha(PREVIOUS/'outputs'/name),role='previous_result') for name in ['cases.json','pilot.json','bluntplate_final.json','rough_final.json']]
    (ROOT/'outputs/cases.json').write_text(json.dumps(cases,indent=2)+'\n')
    (ROOT/'outputs/input_manifest.json').write_text(json.dumps(protected,indent=2)+'\n')
    for name in ['common.py','run_adapter.py']:shutil.copy2(PREVIOUS/'work'/name,ROOT/'work'/name)
    shutil.copy2(PREVIOUS/'work/solver.py',ROOT/'work/core.py')
    shutil.copy2(Path(__file__),ROOT/'work/bootstrap.py')
    (ROOT/'outputs/bootstrap_receipt.json').write_text(json.dumps(dict(created=datetime.datetime.now().isoformat(),previous=str(PREVIOUS),root=str(ROOT),protected_files=len(protected)),indent=2)+'\n')
    print('BOOTSTRAP',ROOT,len(protected),flush=True)
if __name__=='__main__':main()
