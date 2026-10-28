"""Recover only this task's interrupted generated caches; never touch s75."""
import datetime,hashlib,json,shutil,sys,zipfile
from pathlib import Path
import numpy as np
from scipy.sparse import load_npz
ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/neutral_m3_sweep_consistent_20261028')
def main():
    # All owners must be stopped before moving incomplete generated directories.
    for p in Path('/proc').iterdir():
        if not p.name.isdigit():continue
        try:args=(p/'cmdline').read_bytes().decode().split(chr(0))
        except (OSError,UnicodeDecodeError):continue
        assert not any(a.startswith(str(ROOT/'work')+'/') and Path(a).name in ['precompute.py','run_case_fpp.py','contour3d.py','finalize.py'] for a in args),('Task still active',p,args)
    stamp=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).strftime('%Y%m%dT%H%M%S');dest=ROOT/'work/interrupted_cache'/stamp;rows=[]
    def move(p,reason,category):
        p=p.resolve();assert ROOT.resolve() in p.parents and p.is_dir();target=dest/category/p.name;target.parent.mkdir(parents=True,exist_ok=True);assert not target.exists()
        files=[dict(name=str(f.relative_to(p)),sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in p.rglob('*') if f.is_file()]
        shutil.move(str(p),str(target));rows.append(dict(source=str(p),destination=str(target),reason=reason,files=files))
    for p in list((ROOT/'matrix_exports').iterdir()):
        if not p.is_dir() or not p.name.startswith(('s30','s45','s60')):continue
        report=p/'export.json'
        try:
            if not report.exists():raise RuntimeError('No completed export receipt after interruption')
            r=json.loads(report.read_text())
            for f in r['files']:assert Path(f['path']).exists() and Path(f['path']).stat().st_size>16
        except Exception as e:move(p,repr(e),'exports')
    for case in ['s30','s45','s60']:
        count=0
        for p in list((ROOT/'native_pencils'/case).glob('i*')):
            try:
                r=json.loads((p/'coefficients.json').read_text());assert r['case']==case
                for k in range(3):
                    for j in range(3):assert zipfile.is_zipfile(p/f'C{k}{j}.npz')
                if (p/'frequency_coefficient.json').exists():
                    json.loads((p/'frequency_coefficient.json').read_text());assert zipfile.is_zipfile(p/'G.npz')
                count+=1
            except Exception as e:move(p,repr(e),'native_'+case)
        print('CACHE_AUDIT',case,count,flush=True)
    out=ROOT/'outputs'/f'interruption_recovery_{stamp}.json';out.write_text(json.dumps(dict(moved=rows,preserved_case='s75',policy='Only interrupted generated cache directories moved inside this task; original data untouched.'),indent=2)+'\n');print('RECOVERED',len(rows),out,flush=True)
if __name__=='__main__':main()
