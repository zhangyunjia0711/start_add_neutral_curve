"""Verify the bounded delivery and commit it; remote upload is separate."""
import datetime
import hashlib
import json
import subprocess
from pathlib import Path

WIN = Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate/75°')
DEST = WIN / 'outputs/M3_frequency_limits_20261027'
REPO = Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/bluntplate/neutral_beta_full_20261022/work/solver_source')
PREFIX = 'doc/audit/20261027/M3_frequency_limits/'
EXPECTED_PARENT = '7a78d1441d255414ff964c40e2bcba04e0227b8a'


def git(*args):
    result = subprocess.run(['git', '-C', str(REPO), *args], capture_output=True, text=True)
    if result.returncode:
        print('\n'.join((result.stdout + result.stderr).splitlines()[:8]))
        raise RuntimeError(f'git {args[0]} failed: {result.returncode}')
    return result.stdout


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parent = git('rev-parse', 'HEAD').strip()
    assert parent == EXPECTED_PARENT, 'Recheck changed parent before committing'
    branch = git('branch', '--show-current').strip()
    assert branch == 'codex/full-beta-20261022'
    for args in [('diff', '--cached', '--name-only', '-z'),
                 ('diff', '--name-only', '-z'),
                 ('ls-files', '--others', '--exclude-standard', '-z')]:
        assert all(p.startswith(PREFIX) for p in git(*args).split('\0') if p)
    manifest = json.loads((DEST / 'delivery_manifest.json').read_text())
    for entry in manifest['files']:
        assert sha(DEST / entry['path']) == entry['sha256']
        assert sha(REPO / PREFIX / entry['path']) == entry['sha256']
    assert sha(DEST / 'delivery_manifest.json') == sha(REPO / PREFIX / 'delivery_manifest.json')
    git('add', '-f', '--', PREFIX)
    staged = [p for p in git('diff', '--cached', '--name-only', '-z').split('\0') if p]
    assert len(staged) == manifest['file_count'] + 1
    assert all(p.startswith(PREFIX) for p in staged)
    git('diff', '--cached', '--check')
    output = git('commit', '-m', 'Add verified M3 frequency limits and downstream neutral pockets')
    print(output.splitlines()[0])
    head = git('rev-parse', 'HEAD').strip()
    assert git('rev-parse', 'HEAD^').strip() == parent
    assert not git('status', '--porcelain'), 'Unexpected remaining changes'
    receipt = dict(commit=head, parent=parent, branch=branch, staged_files=len(staged),
                   scientific_manifest_sha256=sha(DEST / 'delivery_manifest.json'),
                   total_bytes=manifest['total_bytes'],
                   remote='neutral-backup', remote_status='not_yet_verified',
                   archived_log_scope='Precommit delivery checks; final push receipt is in the main task log.',
                   time=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat())
    (DEST / 'git_version.json').write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + '\n')
    print('COMMIT_VERIFIED', head, len(staged), 'files; working tree clean')


if __name__ == '__main__':
    main()
