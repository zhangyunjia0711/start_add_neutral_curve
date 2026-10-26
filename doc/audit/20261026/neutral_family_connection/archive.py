"""Archive this connectivity correction without modifying earlier results."""
import datetime,hashlib,json,shutil,subprocess,sys
from pathlib import Path
BASE=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°')
sys.path.insert(0,str(BASE/'neutral_xbeta_four_omega_20261026/work'))
from common import snapshot
DEST=Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate/75°/outputs/neutral_frequency_connection_20261026')
REPO=BASE/'bluntplate/neutral_beta_full_20261022/work/solver_source'
TARGET=REPO/'doc/audit/20261026/neutral_family_connection'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    refresh='--refresh' in sys.argv
    if not refresh:
        assert not TARGET.exists(),'Preserve an existing archive'
        assert not subprocess.check_output(['git','-C',str(REPO),'status','--porcelain'],text=True),'Do not mix another task'
    else:
        prefix=str(TARGET.relative_to(REPO))+'/'
        for command in [['diff','--cached','--name-only','-z'],['diff','--name-only','-z'],['ls-files','--others','--exclude-standard','-z']]:
            names=subprocess.check_output(['git','-C',str(REPO)]+command).decode().split('\0')
            assert all(p.startswith(prefix) for p in names if p),'Do not mix another task'
    for name in ['connection_report.json','f100/connection_report.json']:
        report=json.loads((DEST/name).read_text());assert len(report['paths'])==4
        assert all(p['status']=='continuous_neutral_path_verified' for p in report['paths'])
    protected=snapshot();shutil.copy2(__file__,DEST/'archive.py')
    for p in DEST.rglob('*'):
        if p.is_file() and (p.suffix=='.svg' or p.name.startswith('工作日志_')):
            s=p.read_text(encoding='utf-8-sig');p.write_text('\n'.join(v.rstrip() for v in s.splitlines()).rstrip()+'\n',encoding='utf-8')
    (DEST/'.gitattributes').write_text('*.log binary\n*.plt binary\n*.ini binary\n',encoding='ascii')
    files=[dict(path=str(p.relative_to(DEST)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(DEST.rglob('*')) if p.is_file() and p.name not in ['delivery_manifest.json','git_version.json'] and not p.name.startswith('工作日志_')]
    (DEST/'delivery_manifest.json').write_text(json.dumps(dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),files=files,protected_inputs=protected,execution_metadata_excluded=['工作日志_*.md','git_version.json','delivery_manifest.json']),indent=2)+'\n')
    shutil.copytree(DEST,TARGET,dirs_exist_ok=refresh)
    for p in TARGET.rglob('*'):
        if p.is_file():p.chmod(0o644)
    assert all(sha(TARGET/r['path'])==r['sha256'] for r in files)
    print('ARCHIVE',str(TARGET),'files',len(files),'bytes',sum(r['bytes'] for r in files),flush=True)
if __name__=='__main__':main()
