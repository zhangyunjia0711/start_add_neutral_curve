"""Archive this bounded interpretation without modifying any solved result."""
import datetime,hashlib,json,shutil,subprocess
from pathlib import Path
WIN=Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate/75°')
DEST=WIN/'outputs/M3_neutral_interpretation_20261027'
REPO=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/bluntplate/neutral_beta_full_20261022/work/solver_source')
TARGET=REPO/'doc/audit/20261027/M3_neutral_interpretation'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    prefix=str(TARGET.relative_to(REPO))+'/'
    for args in [('diff','--cached','--name-only','-z'),('diff','--name-only','-z'),('ls-files','--others','--exclude-standard','-z')]:
        names=subprocess.check_output(['git','-C',str(REPO),*args]).decode().split('\0')
        assert all(n.startswith(prefix) for n in names if n)
    m=json.loads((DEST/'metrics_provenance.json').read_text());assert all(sha(r['path'])==r['sha256'] for r in m['sources'])
    c=json.loads((DEST/'case_comparison_metadata.json').read_text());assert all(sha(r['config'])==r['config_sha256'] for r in c['cases'])
    log=WIN/'outputs/工作日志_2026-10-27_九频率中性曲线结果解读.md';shutil.copy2(log,DEST/log.name);shutil.copy2(__file__,DEST/'archive.py')
    for p in DEST.iterdir():
        if p.suffix in ['.md','.json','.csv','.svg','.py']:
            s=p.read_text(encoding='utf-8-sig');p.write_text('\n'.join(v.rstrip() for v in s.splitlines()).rstrip()+'\n',encoding='utf-8')
    files=[dict(path=p.name,sha256=sha(p),bytes=p.stat().st_size) for p in sorted(DEST.iterdir()) if p.is_file() and p.name not in ['delivery_manifest.json','git_version.json'] and not p.name.startswith('工作日志_')]
    (DEST/'delivery_manifest.json').write_text(json.dumps(dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),files=files,excluded=['工作日志_*.md','git_version.json','delivery_manifest.json']),ensure_ascii=False,indent=2)+'\n')
    shutil.copytree(DEST,TARGET,dirs_exist_ok=True)
    for p in TARGET.rglob('*'):
        if p.is_file():p.chmod(0o644)
    assert all(sha(TARGET/r['path'])==r['sha256'] for r in files)
    subprocess.run(['git','-C',str(REPO),'add','--',prefix],check=True)
    names=subprocess.check_output(['git','-C',str(REPO),'diff','--cached','--name-only','-z']).decode().split('\0');assert all(n.startswith(prefix) for n in names if n)
    subprocess.run(['git','-C',str(REPO),'diff','--cached','--check'],check=True)
    subprocess.run(['git','-C',str(REPO),'commit','-q','-m','Quantify M3 neutral boundaries and clarify geometric comparison limits'],check=True)
    print('COMMIT',subprocess.check_output(['git','-C',str(REPO),'rev-parse','HEAD'],text=True).strip(),len(files),'science files',sum(r['bytes'] for r in files),'bytes')
if __name__=='__main__':main()
