"""Archive verified curves, raw START checks, and bounded diagnostic evidence."""
import datetime
import hashlib
import json
import platform
import shutil
import subprocess
import sys
from pathlib import Path
import numpy as np
import scipy

WIN = Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate/75°')
ROOT = Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/neutral_xomega_multiple_beta_20261027')
DEST = WIN / 'outputs/M3_xomega_multiple_beta_20261027'
LOG = WIN / 'outputs/工作日志_2026-10-27_多β的xω中性曲线.md'
REPO = Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/bluntplate/neutral_beta_full_20261022/work/solver_source')
PREFIX = 'doc/audit/20261027/M3_xomega_multiple_beta/'
PARENT = 'd83049409f89c394f4215ea18595c573af712523'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def js(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')


def git(*args):
    r = subprocess.run(['git', '-C', str(REPO), *args], capture_output=True, text=True)
    if r.returncode:
        print('\n'.join((r.stdout + r.stderr).splitlines()[:8]))
        raise RuntimeError(f'git {args[0]} failed {r.returncode}')
    return r.stdout


def runs_in(value):
    if isinstance(value, dict):
        for key, item in value.items():
            if key == 'run' and isinstance(item, str) and (ROOT / 'runs' / item).is_dir():
                yield item
            else:
                yield from runs_in(item)
    elif isinstance(value, list):
        for item in value:
            yield from runs_in(item)


def main():
    sys.path.insert(0, str(ROOT / 'work'))
    from common import snapshot
    protected = snapshot()
    assert git('rev-parse', 'HEAD').strip() == PARENT
    assert git('branch', '--show-current').strip() == 'codex/full-beta-20261022'
    for args in [('diff', '--cached', '--name-only', '-z'), ('diff', '--name-only', '-z'),
                 ('ls-files', '--others', '--exclude-standard', '-z')]:
        assert all(p.startswith(PREFIX) for p in git(*args).split('\0') if p)
    selection = json.loads((DEST / 'selection_manifest.json').read_text())
    audit = json.loads((DEST / 'numerical_audit.json').read_text())
    validation = json.loads((DEST / 'tecplot_validation.json').read_text())
    assert selection['curve_count'] == 15 and selection['points'] == 6131 and audit['all_corrected_seeds_covered']
    for r in validation['files']:
        assert sha(DEST / r['file']) == r['sha256'] and sha((DEST / r['file']).with_suffix('.dat')) == r['source_sha256']
    evidence, scripts = DEST / 'evidence', DEST / 'reproducibility'
    evidence.mkdir(exist_ok=True)
    scripts.mkdir(exist_ok=True)
    selected_guides, receipts, vectors, provenance = set(), set(), {}, []
    for k, c in enumerate(selection['curves']):
        source = Path(c['source'])
        assert sha(source) == c['source_sha256']
        shutil.copy2(source, evidence / source.name)
        if c['legacy']:
            old = next(r for r in audit['curves'] if r['source'] == c['source'])
            for j, r in enumerate(old['current_joint_matrix_checks']):
                vectors[f'legacy{k:02d}_{j:02d}'] = np.load(r['vector_file'])
            for p in source.parent.glob('*.json'):
                if any(s in p.name for s in ['summary', 'audit', 'verification', 'turning']):
                    target = evidence / 'legacy' / c['case'] / p.name
                    target.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(p, target)
            continue
        d = json.loads(source.read_text())
        selected_guides.add(d['guide'])
        if d.get('dense_source'):
            p = Path(d['dense_source'])
            assert sha(p) == d['dense_source_sha256']
            shutil.copy2(p, evidence / p.name)
        receipts.update(runs_in(d))
        points = d['points'][:-1]
        indices = sorted(set(np.linspace(0, len(points) - 1, 12, dtype=int).tolist() +
                             [j for j, p in enumerate(points) if p['kind'] == 1]))
        for j, p in enumerate(points):
            provenance.append(dict(curve=k, point=j, path=p['vector_file'], sha256=sha(p['vector_file'])))
        for j in indices:
            vectors[f'curve{k:02d}_{j:04d}'] = np.load(points[j]['vector_file'])
    for path in sorted(selected_guides):
        p = Path(path)
        shutil.copy2(p, evidence / p.name)
    diagnostics = []
    for p in sorted((ROOT / 'outputs').glob('*.json')):
        if str(p) in selected_guides:
            continue
        if p.name.endswith('_guide.json'):
            entry = dict(path=str(p), sha256=sha(p), bytes=p.stat().st_size)
            try:
                d = json.loads(p.read_text())
                entry.update(status=d['status'], record_count=len(d['records']), failures=d['failures'],
                             last_point=d['records'][-1]['point'], closure=d.get('closure'),
                             note='Full diagnostic state and vectors retained locally; not a delivered contour')
            except Exception as error:
                entry.update(status='partial_diagnostic', reason=repr(error))
            diagnostics.append(entry)
        else:
            shutil.copy2(p, evidence / p.name)
    js(evidence / 'other_guide_diagnostics.json', dict(guides=diagnostics))
    for p in sorted((ROOT / 'outputs').glob('*.log')):
        shutil.copy2(p, evidence / p.name)
    for name in ['rough_positive_beta_limit_refined.json', 'rough_positive_beta_limit.json']:
        d = json.loads((ROOT / 'outputs' / name).read_text())
        receipts.update(runs_in(d))
        vectors[name.replace('.json', '')] = np.load(d['vector_file'])
    raw = []
    for run in sorted(receipts):
        directory = ROOT / 'runs' / run
        for rel in ['config.ini', 'summary.json', 'solver.log', 'wall_local_profile_metadata.csv',
                    'out/neutral_field.plt', 'out/all_vectors.petsc']:
            p = directory / rel
            if not p.exists():
                continue
            target = DEST / 'evidence_runs' / run / rel
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(p, target)
            raw.append(dict(path=str(target.relative_to(DEST)), source=str(p), bytes=p.stat().st_size, sha256=sha(p)))
    np.savez_compressed(DEST / 'selected_vectors.npz', **vectors)
    js(DEST / 'point_vector_provenance.json', dict(points=provenance, full_vectors_retained_at=str(ROOT), archived_selected_vectors=len(vectors)))
    js(DEST / 'original_run_evidence_manifest.json', dict(runs=len(receipts), files=raw))
    js(DEST / 'runtime_versions.json', dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__, platform=platform.platform()))
    js(DEST / 'source_preservation.json', dict(protected_inputs=protected, selected_sources_unchanged=True))
    for p in sorted((ROOT / 'work').glob('*.py')):
        shutil.copy2(p, scripts / p.name)
    for p in sorted(Path(__file__).parent.glob('*.py')):
        shutil.copy2(p, scripts / p.name)
    shutil.copy2(__file__, DEST / 'archive.py')
    shutil.copy2(LOG, DEST / LOG.name)
    (DEST / '.gitattributes').write_text('* -text\n*.md text eol=lf\n*.py text eol=lf\n*.lay -whitespace\n*.log -whitespace\n*.csv -whitespace\n*.json -whitespace\n*.svg -whitespace\nevidence_runs/** -whitespace\n')
    files = [dict(path=str(p.relative_to(DEST)), sha256=sha(p), bytes=p.stat().st_size)
             for p in sorted(DEST.rglob('*')) if p.is_file() and p not in [DEST / 'delivery_manifest.json', DEST / 'git_version.json']]
    js(DEST / 'delivery_manifest.json', dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),
       files=files, file_count=len(files), total_bytes=sum(r['bytes'] for r in files), git_scope=PREFIX))
    for p in sorted(DEST.rglob('*')):
        if not p.is_file() or p == DEST / 'git_version.json':
            continue
        target = REPO / PREFIX / p.relative_to(DEST)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, target)
        assert sha(p) == sha(target)
    git('add', '-f', '--', PREFIX)
    staged = [p for p in git('diff', '--cached', '--name-only', '-z').split('\0') if p]
    assert len(staged) == len(files) + 1 and all(p.startswith(PREFIX) for p in staged)
    git('diff', '--cached', '--check')
    print(git('commit', '-m', 'Compute verified M3 x-omega contours at multiple beta values').splitlines()[0])
    head = git('rev-parse', 'HEAD').strip()
    assert git('rev-parse', 'HEAD^').strip() == PARENT and not git('status', '--porcelain')
    js(DEST / 'git_version.json', dict(commit=head, parent=PARENT, branch='codex/full-beta-20261022',
       staged_files=len(staged), total_bytes=sum(r['bytes'] for r in files), remote='neutral-backup',
       remote_status='not_yet_verified', manifest_sha256=sha(DEST / 'delivery_manifest.json'),
       archived_log_scope='Delivery checks before commit; final push receipt is in the main task log'))
    print('ARCHIVE_COMMIT_VERIFIED', head, len(staged), 'files', sum(r['bytes'] for r in files), 'bytes', len(receipts), 'frozen runs', flush=True)


if __name__ == '__main__':
    main()
