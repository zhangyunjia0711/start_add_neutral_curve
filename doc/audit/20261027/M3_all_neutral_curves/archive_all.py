"""Archive an exact, bounded combination of already verified curve data."""
import datetime
import hashlib
import json
import shutil
import subprocess
from pathlib import Path

ROOT = Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate/75°')
DEST = ROOT / 'outputs/M3_all_neutral_curves_20261027'
LOG = ROOT / 'outputs/工作日志_2026-10-27_全部M3中性曲线合图.md'
REPO = Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/bluntplate/neutral_beta_full_20261022/work/solver_source')
PREFIX = 'doc/audit/20261027/M3_all_neutral_curves/'
PARENT = '1e7abcb7c295a825696fc626ef90528c90f93d1e'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def js(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')


def git(*args):
    result = subprocess.run(['git', '-C', str(REPO), *args], capture_output=True, text=True)
    if result.returncode:
        print('\n'.join((result.stdout + result.stderr).splitlines()[:8]))
        raise RuntimeError(f'git {args[0]} failed {result.returncode}')
    return result.stdout


def main():
    assert git('rev-parse', 'HEAD').strip() == PARENT
    assert git('branch', '--show-current').strip() == 'codex/full-beta-20261022'
    for args in [('diff', '--cached', '--name-only', '-z'), ('diff', '--name-only', '-z'),
                 ('ls-files', '--others', '--exclude-standard', '-z')]:
        assert all(p.startswith(PREFIX) for p in git(*args).split('\0') if p)
    selection = json.loads((DEST / 'selection_manifest.json').read_text())
    validation = json.loads((DEST / 'tecplot_validation.json').read_text())
    assert selection['curve_count'] == 45 and selection['curve_points'] == 14807
    assert selection['frequency_count'] == 25 and selection['critical_points'] == 4
    protected = []
    for c in selection['curves']:
        source = c['source']
        assert sha(source['path']) == sha(DEST / c['output_csv']) == source['sha256']
        protected.append(source)
    for d in selection['duplicates_removed']:
        assert sha(d['source']['path']) == d['source']['sha256']
        protected.append(d['source'])
    source = selection['critical_source']
    assert sha(source['path']) == source['sha256']
    protected.append(source)
    for r in validation['files']:
        assert r['all_values_exact'] and sha(DEST / r['file']) == r['sha256']
        assert sha((DEST / r['file']).with_suffix('.dat')) == r['source_sha256']
    for folder in ['M3_xbeta_extended_frequencies_20261027', 'M3_frequency_limits_20261027']:
        for name in ['selection_manifest.json', 'delivery_manifest.json']:
            source = ROOT / 'outputs' / folder / name
            target = DEST / 'source_manifests' / folder / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
            protected.append(dict(path=str(source), sha256=sha(source)))
    js(DEST / 'source_preservation.json', dict(all_checked_sources_unchanged=True, files=protected))
    shutil.copy2(__file__, DEST / 'archive_all.py')
    shutil.copy2(LOG, DEST / LOG.name)
    (DEST / '.gitattributes').write_text('* -text\n*.md text eol=lf\n*.py text eol=lf\n*.lay -whitespace\n*.log -whitespace\n*.csv -whitespace\n*.json -whitespace\n*.svg -whitespace\n')
    files = [dict(path=str(p.relative_to(DEST)), bytes=p.stat().st_size, sha256=sha(p))
             for p in sorted(DEST.rglob('*')) if p.is_file() and p not in [DEST / 'delivery_manifest.json', DEST / 'git_version.json']]
    js(DEST / 'delivery_manifest.json', dict(time=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),
       files=files, file_count=len(files), total_bytes=sum(p['bytes'] for p in files), git_scope=PREFIX))
    for source in sorted(DEST.rglob('*')):
        if not source.is_file() or source.name == 'git_version.json':
            continue
        target = REPO / PREFIX / source.relative_to(DEST)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        assert sha(source) == sha(target)
    git('add', '-f', '--', PREFIX)
    staged = [p for p in git('diff', '--cached', '--name-only', '-z').split('\0') if p]
    assert len(staged) == len(files) + 1 and all(p.startswith(PREFIX) for p in staged)
    git('diff', '--cached', '--check')
    print(git('commit', '-m', 'Plot all computed M3 neutral curves together').splitlines()[0])
    head = git('rev-parse', 'HEAD').strip()
    assert git('rev-parse', 'HEAD^').strip() == PARENT and not git('status', '--porcelain')
    js(DEST / 'git_version.json', dict(commit=head, parent=PARENT, branch='codex/full-beta-20261022',
       staged_files=len(staged), total_bytes=sum(p['bytes'] for p in files),
       remote='neutral-backup', remote_status='not_yet_verified',
       manifest_sha256=sha(DEST / 'delivery_manifest.json'),
       archived_log_scope='Delivery checks before commit; final push receipt is in the main task log.'))
    print('ARCHIVE_COMMIT_VERIFIED', head, len(staged), 'files', sum(p['bytes'] for p in files), 'bytes; worktree clean')


if __name__ == '__main__':
    main()
