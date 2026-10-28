"""Verify and commit only this task's source changes and immutable delivery."""
import datetime,hashlib,json,shutil,subprocess
from pathlib import Path
BASE=Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate')
DEST=BASE/'outputs/M3_sweep_angles_20261028'
REPO=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/bluntplate/neutral_beta_full_20261022/work/solver_source')
REL=Path('doc/audit/20261028/M3_sweep_angles')
OWNED=['src/public/wall_local_coordinates.f90','tests/wall_local_profiles_test.f90']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(args):
    p=subprocess.run(['git','-c','core.quotepath=false']+args,cwd=REPO,capture_output=True,text=True)
    assert p.returncode==0,(args,p.stdout,p.stderr)
    return p.stdout
def main():
    assert run(['branch','--show-current']).strip()=='codex/full-beta-20261022'
    assert run(['rev-parse','HEAD']).strip()=='b391bde3184540a820e554b153082370e6ef3819'
    if not (REPO/REL).exists():
        assert not run(['diff','--cached','--name-only']).strip(),'Do not include another staged change'
        assert set(run(['diff','--name-only']).splitlines())==set(OWNED)
    else:
        # Resume only an exact, task-owned archive already staged by this script.
        previous=json.loads((REPO/REL/'delivery_manifest.json').read_text())
        assert all(sha(REPO/REL/r['path'])==r['sha256'] for r in previous['files'])
        allowed=set(OWNED+[str(REL/'.gitattributes'),str(REL/'delivery_manifest.json')]+[str(REL/r['path']) for r in previous['files']])
        assert set(run(['diff','--cached','--name-only']).splitlines())==allowed
        assert not run(['diff','--name-only']).strip(),'Do not overwrite unstaged edits during resume'
    receipt=json.loads((DEST/'source_application.json').read_text())
    assert all(sha(REPO/r['path'])==r['after_sha256'] for r in receipt['files'])
    # Keep this orchestration script as part of the reproducible delivery too.
    shutil.copy2(BASE/'work/M3_sweep_angles_20261027/archive_commit.py',DEST/'reproducibility/runtime/archive_commit.py')
    records=[]
    for p in sorted(DEST.rglob('*')):
        if p.is_file() and str(p.relative_to(DEST)) not in ['delivery_manifest.json','git_version.json']:
            records.append(dict(path=str(p.relative_to(DEST)),bytes=p.stat().st_size,sha256=sha(p)))
    manifest=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),files=records,excluded_root_files=['delivery_manifest.json','git_version.json'],scope='Immutable delivery snapshot before commit; final remote SHA is recorded outside this self-referential snapshot.')
    (DEST/'delivery_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    shutil.copytree(DEST,REPO/REL,dirs_exist_ok=True)
    # Preserve data and frozen script bytes exactly, including intentional CRLF.
    (REPO/REL/'.gitattributes').write_text('* -text -whitespace\n')
    assert all(sha(REPO/REL/r['path'])==r['sha256'] for r in records)
    run(['add','--']+OWNED+[str(REL)])
    run(['diff','--cached','--check'])
    staged=run(['diff','--cached','--name-only']).splitlines()
    assert set(staged)==set(OWNED+[str(REL/'.gitattributes'),str(REL/'delivery_manifest.json')]+[str(REL/r['path']) for r in records])
    output=run(['commit','-m','Add consistent normal derivatives and verified M3 sweep-angle curves','-m','Compute 141 curves from eight dense FPP cases and retain 23 original 75-degree curves without recomputation. Archive figures, native Tecplot files, source checks and numerical evidence; record unresolved frequencies and missing inputs explicitly.'])
    commit=run(['rev-parse','HEAD']).strip()
    assert not run(['status','--porcelain']).strip()
    result=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),commit=commit,branch='codex/full-beta-20261022',remote='neutral-backup',remote_url='https://github.com/zhangyunjia0711/start_add_neutral_curve.git',archive_path=str(REL),delivery_manifest_sha256=sha(DEST/'delivery_manifest.json'),committed_paths=len(staged),worktree_clean=True,remote_verified=False)
    (DEST/'git_version.json').write_text(json.dumps(result,indent=2)+'\n')
    print('COMMITTED',commit,len(staged),'owned paths',flush=True)

if __name__=='__main__':main()
