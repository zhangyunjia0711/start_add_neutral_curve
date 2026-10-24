"""Archive final outputs and compact reproducibility evidence without matrix caches."""
import datetime,gzip,hashlib,json,shutil,subprocess
from pathlib import Path
from common import ROOT,PREVIOUS,sha,snapshot

def main():
    source=PREVIOUS/'work/solver_source'
    target=source/'doc/audit/20261024/M3_fixed_beta_closed'
    assert not target.exists(),'Existing audit archive must not be overwritten'
    assert not subprocess.check_output(['git','-C',str(source),'status','--porcelain'])
    snapshot()
    outputs=target/'outputs';scripts=target/'scripts'
    outputs.mkdir(parents=True);scripts.mkdir()
    for p in (ROOT/'work').glob('*.py'):shutil.copy2(p,scripts/p.name)
    included=[]
    for p in sorted((ROOT/'outputs').iterdir()):
        if p.is_file() and (p.suffix in {'.json','.csv','.dat','.png','.svg','.plt','.lay','.mcr','.md'} or p.name=='batch.log'):
            shutil.copy2(p,outputs/p.name);included.append(dict(file=str(Path('outputs')/p.name),bytes=p.stat().st_size,sha256=sha(p)))
    adapter=scripts/'full_driver_dependency_snapshot.py'
    shutil.copy2(PREVIOUS/'work/full_driver.py',adapter)
    bundle=outputs/'native_point_evidence.jsonl.gz'
    files=sorted((ROOT/'outputs/native_points').glob('*.json'))
    assert len(files)==1565
    with bundle.open('wb') as raw,gzip.GzipFile(fileobj=raw,mode='wb',mtime=0) as zipped:
        for p in files:
            value=json.loads(p.read_text());q=Path(value['vector_file'])
            entry=dict(path=str(p),sha256=sha(p),record=value,vector_sha256=sha(q))
            zipped.write((json.dumps(entry,separators=(',',':'),allow_nan=False)+'\n').encode())
    included.append(dict(file=str(Path('outputs')/bundle.name),bytes=bundle.stat().st_size,sha256=sha(bundle)))
    report=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),
        source_head_before_task=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip(),
        original_task_root=str(ROOT),archive_relative_path=str(target.relative_to(source)),
        output_files=included,native_point_evidence_records=len(files),
        raw_large_evidence='Raw CFD inputs, PETSc matrices, necessary raw vectors, and full solver runs remain at the original task paths; they are not bundled into this Git audit archive.',
        original_inputs=snapshot(),adapter_sha256=sha(adapter),
        script_files=[dict(file=str(p.relative_to(target)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(scripts.glob('*.py'))])
    (target/'archive_manifest.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    (target/'README.md').write_text('# M3 fixed beta=-7 closed neutral curve\n\n'
        'See outputs/结果说明.md and the dated work log for results, model scope, validation, and source paths. '
        'The scripts are snapshots; execute in the preserved task directory described there. '
        'full_driver_dependency_snapshot.py records the run-adapter dependency. '
        'native_point_evidence.jsonl.gz contains all 1565 native-root jobs with raw-vector SHA256 hashes. '
        'Large matrices/vectors and source CFD remain local.\n')
    print(json.dumps(dict(archive=str(target),outputs=len(included),native_point_records=len(files),bytes=sum(p.stat().st_size for p in target.rglob('*') if p.is_file()))))

if __name__=='__main__':main()
