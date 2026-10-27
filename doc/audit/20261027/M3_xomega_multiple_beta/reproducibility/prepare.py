"""Prepare an isolated, fixed-beta calculation using existing frozen coefficients."""
import datetime
import hashlib
import json
import shutil
from pathlib import Path

BASE = Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°')
SOURCE = BASE / 'neutral_m3_frequency_limits_20261027'
ROOT = BASE / 'neutral_xomega_multiple_beta_20261027'
WIN = Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate/75°')
BETAS = [-30., -20., -10., -7., 0., 10., 20.]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    assert not ROOT.exists(), 'Preserve existing task'
    for name in ['work', 'outputs', 'runs', 'inputs', 'matrix_exports']:
        (ROOT / name).mkdir(parents=True)
    for name in ['cases', 'native_pencils']:
        (ROOT / name).symlink_to(SOURCE / name, target_is_directory=True)
    files = []
    for name in ['common.py', 'core.py', 'solver.py', 'physical.py', 'run_adapter.py']:
        source = SOURCE / 'work' / name
        shutil.copy2(source, ROOT / 'work' / name)
        files.append(dict(path=str(source), sha256=sha(source)))
    protected = json.loads((SOURCE / 'outputs/input_manifest.json').read_text())
    for folder in ['M3_all_neutral_curves_20261027', 'comparison_M3_beta_minus7_20261025', 'M3_beta_minus7_closed_20261024']:
        for name in ['selection_manifest.json', 'delivery_manifest.json']:
            path = WIN / 'outputs' / folder / name
            if path.exists():
                protected.append(dict(path=str(path), sha256=sha(path), role='prior_delivery'))
    for name in ['cases.json', 'frequency_conversion.json']:
        shutil.copy2(SOURCE / 'outputs' / name, ROOT / 'outputs' / name)
    (ROOT / 'outputs/input_manifest.json').write_text(json.dumps(protected, indent=2) + '\n')
    (ROOT / 'outputs/bootstrap_receipt.json').write_text(json.dumps(dict(
        time=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),
        betas=BETAS, root=str(ROOT), source=str(SOURCE), protected_files=len(protected),
        source_scripts=files, model='Unchanged frozen wall-local, zero-explicit-curvature LST'), indent=2) + '\n')
    print('PREPARED', ROOT, len(protected), flush=True)


if __name__ == '__main__':
    main()
