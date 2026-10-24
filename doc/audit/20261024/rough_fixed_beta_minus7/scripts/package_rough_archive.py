"""Compact audit archive for the computed rough-case fixed-beta curve."""
import datetime,gzip,json,shutil,subprocess
from pathlib import Path
from common import ROOT,sha,snapshot

SOURCE=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/bluntplate/neutral_beta_full_20261022/work/solver_source')

def main():
    target=SOURCE/'doc/audit/20261024/rough_fixed_beta_minus7';assert not target.exists()
    assert not subprocess.check_output(['git','-C',str(SOURCE),'status','--porcelain'])
    snapshot();outputs=target/'outputs';scripts=target/'scripts';outputs.mkdir(parents=True);scripts.mkdir()
    for p in (ROOT/'work').glob('*.py'):shutil.copy2(p,scripts/p.name)
    for p in sorted((ROOT/'outputs').iterdir()):
        if p.is_file() and p.suffix in {'.json','.csv','.dat','.png','.svg','.plt','.lay','.mcr','.md','.log'}:shutil.copy2(p,outputs/p.name)
    points=sorted((ROOT/'outputs/native_points').glob('*.json'));assert len(points)==644
    with (outputs/'native_point_evidence.jsonl.gz').open('wb') as raw,gzip.GzipFile(fileobj=raw,mode='wb',mtime=0) as zipped:
        for p in points:
            record=json.loads(p.read_text());vector=Path(record['vector_file'])
            entry=dict(path=str(p),sha256=sha(p),record=record,vector_sha256=sha(vector))
            zipped.write((json.dumps(entry,separators=(',',':'),allow_nan=False)+'\n').encode())
    manifest=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),
        source_head_before_task=subprocess.check_output(['git','-C',str(SOURCE),'rev-parse','HEAD'],text=True).strip(),
        original_task_root=str(ROOT),archive_relative_path=str(target.relative_to(SOURCE)),native_point_records=644,original_inputs=snapshot(),
        large_evidence='Original CFD inputs, full solver runs, needed raw vectors and PETSc matrices remain at preserved local paths. Superseded j-index seed experiments are excluded from the final M3 result and their discarded files are listed in cleanup_receipt.json.',
        files=[dict(file=str(p.relative_to(target)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(target.rglob('*')) if p.is_file()])
    (target/'archive_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (target/'README.md').write_text('# Rough case: M3 fixed beta=-7\n\n'
        'One verified closed neutral component: 644 native points on 292 stations, 6 re-solved local-profile turning points, '
        'x=-34.932115732 to 623.608951387 mm. See outputs/结果说明.md and the dated Chinese work log. '
        'Native Tecplot binary was read back and all 5859 double values matched exactly. '
        'Scripts are snapshots for the preserved calculation-directory structure; large original data/cache stay local. '
        'This is the declared wall-local parallel model with zero explicit curvature.\n')
    print(json.dumps(dict(archive=str(target),files=sum(p.is_file() for p in target.rglob('*')),bytes=sum(p.stat().st_size for p in target.rglob('*') if p.is_file()))))

if __name__=='__main__':main()
