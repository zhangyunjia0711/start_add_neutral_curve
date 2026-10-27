"""Separate run for complete positive/negative frequency limits of tracked M3."""
import datetime,hashlib,json,shutil
from pathlib import Path
BASE=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°')
SRC=BASE/'neutral_xbeta_extended_frequencies_20261027'
ROOT=BASE/'neutral_m3_frequency_limits_20261027'
WIN=Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate/75°')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    assert not ROOT.exists(),'Preserve a prior calculation'
    for name in ['work','outputs','runs','matrix_exports','inputs']:(ROOT/name).mkdir(parents=True,exist_ok=True)
    for name in ['cases','native_pencils']:(ROOT/name).symlink_to(SRC/name,target_is_directory=True)
    sources=[]
    for name in ['common.py','run_adapter.py','core.py','solver.py','physical.py','contour.py','contour3d.py','finalize.py','tip_correction.py','physical_fold_bordered.py']:
        source=SRC/'work'/name;text=source.read_text();sources.append(dict(path=str(source),sha256=sha(source)))
        if name in ['contour3d.py','finalize.py']:text=text.replace("type=int,choices=[-112, -110, 125, 150, 175, 200]","type=float")
        (ROOT/'work'/name).write_text(text)
    protected=json.loads((SRC/'outputs/input_manifest.json').read_text())
    protected += [dict(path=str(p),sha256=sha(p),role='prior_extended_frequency_curve') for p in sorted((WIN/'outputs/M3_xbeta_extended_frequencies_20261027').glob('*.csv'))]
    (ROOT/'outputs/input_manifest.json').write_text(json.dumps(protected,indent=2)+'\n')
    for name in ['cases.json','frequency_conversion.json']:shutil.copy2(SRC/'outputs'/name,ROOT/'outputs'/name)
    for p in Path(__file__).parent.glob('*.py'):shutil.copy2(p,ROOT/'work'/p.name)
    (ROOT/'outputs/bootstrap_receipt.json').write_text(json.dumps(dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),source=str(SRC),root=str(ROOT),protected_files=len(protected),sources=sources,model='Unchanged frozen wall-local zero-explicit-curvature LST'),indent=2)+'\n');print('PREPARED',ROOT,len(protected),flush=True)
if __name__=='__main__':main()
