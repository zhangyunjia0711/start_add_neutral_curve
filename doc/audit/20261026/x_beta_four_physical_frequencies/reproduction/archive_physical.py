"""Archive only this task's final, small deliverables into the existing Git worktree."""
import argparse,hashlib,json,shutil,subprocess,datetime
from pathlib import Path
from common import ROOT,snapshot
DEST=Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate/75°/outputs/M3_xbeta_four_omega_20261026')
REPO=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/bluntplate/neutral_beta_full_20261022/work/solver_source')
TARGET=REPO/'doc/audit/20261026/x_beta_four_physical_frequencies'

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--refresh',action='store_true');args=parser.parse_args()
    snapshot()
    if not args.refresh:
        assert not TARGET.exists(),'Preserve existing archive'
        status=subprocess.check_output(['git','-C',str(REPO),'status','--porcelain'],text=True);assert not status,'Do not mix other changes'
    else:
        prefix=str(TARGET.relative_to(REPO))+'/'
        for command in [['diff','--cached','--name-only','-z'],['diff','--name-only','-z'],['ls-files','--others','--exclude-standard','-z']]:
            names=subprocess.check_output(['git','-C',str(REPO)]+command).decode().split('\0')
            assert all(name.startswith(prefix) for name in names if name),'Do not mix other changes'
    # Normalize generated CSV/layout and execution-log line endings for Git.
    # Preserve binary data and native logs whose hashes certify the Tecplot run.
    for p in DEST.iterdir():
        if p.suffix in ['.csv','.lay','.svg','.json'] or p.name.startswith('工作日志_'):
            text=p.read_text(encoding='utf-8-sig');p.write_text('\n'.join(line.rstrip() for line in text.splitlines())+'\n',encoding='utf-8')
    (DEST/'.gitattributes').write_text('*.plt binary\n*.log binary\n# Preserve byte-identical original solver receipts.\noriginal_checks/** binary\n',encoding='ascii')
    # The calculation log and future Git receipt are intentionally outside this
    # immutable scientific-file hash set, since they record archive execution.
    files=[]
    for p in sorted(DEST.rglob('*')):
        if not p.is_file() or p.name=='delivery_manifest.json' or p.name=='git_version.json' or p.name.startswith('工作日志_'):continue
        files.append(dict(path=str(p.relative_to(DEST)),bytes=p.stat().st_size,sha256=digest(p)))
    manifest=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),files=files,excluded_execution_metadata=['工作日志_*.md','git_version.json','delivery_manifest.json'],calculation_root=str(ROOT),protected_inputs=snapshot())
    (DEST/'delivery_manifest.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    shutil.copytree(DEST,TARGET,dirs_exist_ok=args.refresh)
    for p in TARGET.rglob('*'):
        if p.is_file():p.chmod(0o644)
    assert all(digest(TARGET/q['path'])==q['sha256'] for q in files)
    print(json.dumps(dict(target=str(TARGET),files=len(files),bytes=sum(q['bytes'] for q in files)),ensure_ascii=False))

if __name__=='__main__':main()
