"""Sequential curves within a case, with independently reviewable stage logs."""
import argparse,json,os,subprocess,sys,datetime
from pathlib import Path
ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/neutral_m3_sweep_consistent_20261028')
ENV=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
def main(case,omit0=False):
    assert case in ['s30','s45','s50','s55','s60','s65','s70','s73'], '75 degrees must only reuse existing data'
    source=json.loads((ROOT/'outputs'/f'{case}_neutral_seed_scan.json').read_text())
    jobs=sorted(set((r['f_kHz'],r.get('component','S')) for r in source['records'] if r['side']=='lower'),key=lambda t:(abs(t[0]),t[0],t[1]))
    for f,component in jobs:
        if omit0 and f==0:continue
        tag=f'{case}_f{f}kHz_{component}';guide=ROOT/'outputs'/f'{tag}_v3_guide.json';final=ROOT/'outputs'/f'{tag}_final.json';steps=[]
        if final.exists():continue
        if not guide.exists():steps.append('contour3d.py')
        elif json.loads(guide.read_text())['status']=='running':steps.append('contour3d.py')
        steps.append('finalize.py')
        for script in steps:
            if script=='finalize.py':
                d=json.loads(guide.read_text())
                if d['status'] not in ['closed_verified','open_at_domain_boundaries']:print('UNRESOLVED',tag,d['status'],flush=True);break
            logfile=ROOT/'work'/f'{tag}_{script[:-3]}.log';command=[sys.executable,str(ROOT/'work'/script),case,str(f),component]
            if script=='contour3d.py' and guide.exists():command+=['--resume']
            with logfile.open('a') as out:r=subprocess.run(command,env=ENV,stdout=out,stderr=subprocess.STDOUT)
            print('STEP',tag,script,r.returncode,flush=True)
            if r.returncode:break
    print('CASE_STAGE_FINISHED',case,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case');p.add_argument('--omit-zero',action='store_true');a=p.parse_args();main(a.case,a.omit_zero)
