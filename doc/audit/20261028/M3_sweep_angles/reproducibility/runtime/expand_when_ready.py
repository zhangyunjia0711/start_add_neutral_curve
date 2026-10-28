"""Finish common-frequency curves after the validated stationary guide is ready."""
import argparse,json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/neutral_m3_sweep_consistent_20261028')
ENV=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
def main(case):
    assert case in ['s50','s55','s65','s70','s73']
    guide=ROOT/'outputs'/f'{case}_f0.0kHz_S_v3_guide.json'
    while True:
        if guide.exists():
            try:d=json.loads(guide.read_text())
            except ValueError:time.sleep(2);continue
            if d['status'] in ['closed_verified','open_at_domain_boundaries']:break
            if d['status']!='running':raise RuntimeError('Stationary guide unresolved: '+d['status'])
        time.sleep(2)
    print('STATIONARY_GUIDE_READY',case,flush=True)
    with (ROOT/'work'/f'{case}_expanded_seed_search.log').open('a') as log:
        helper=subprocess.Popen([sys.executable,str(ROOT/'work/frequency_seeds_fpp.py'),case],env=ENV,stdout=log,stderr=subprocess.STDOUT)
        with (ROOT/'work'/f'{case}_expanded_queue.log').open('a') as out:
            worker=subprocess.Popen([sys.executable,str(ROOT/'work/run_batch.py'),case,'expanded'],env=ENV,stdout=out,stderr=subprocess.STDOUT)
            assert helper.wait()==0
            assert worker.wait()==0
    print('EXPANDED_CURVES_FINISHED',case,flush=True)
    with (ROOT/'work'/f'{case}_supplemental_search.log').open('a') as log:
        r=subprocess.run([sys.executable,str(ROOT/'work/additional_seed_probe.py'),case,'--frequencies','-75','-100','-110','-112'],env=ENV,stdout=log,stderr=subprocess.STDOUT)
        assert r.returncode==0
    print('CASE_PIPELINE_FINISHED',case,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case');a=p.parse_args();main(a.case)
