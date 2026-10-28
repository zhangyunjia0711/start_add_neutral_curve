"""Task-owned numerical subprocess queue; no 75-degree solver calls."""
import argparse,datetime,fcntl,json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/neutral_m3_sweep_consistent_20261028')
ENV=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')

def processes():
    result=[]
    for p in Path('/proc').glob('[0-9]*/cmdline'):
        try:result.append(p.read_bytes().split(b'\0')[:-1])
        except OSError:pass
    return result

def active_curve(case,f,component):
    scripts=[str(ROOT/'work'/s).encode() for s in ['contour3d.py','finalize.py']]
    return any(len(a)>=5 and a[1] in scripts and a[2:5]==[case.encode(),str(float(f)).encode(),component.encode()] for a in processes())

def seed_helper(case):
    name=str(ROOT/'work/frequency_seeds_fpp.py').encode()
    return any(len(a)>2 and a[1:3]==[name,case.encode()] for a in processes())

def run_one(case,f,component):
    tag=f'{case}_f{float(f)}kHz_{component}';guide=ROOT/'outputs'/f'{tag}_v3_guide.json';final=ROOT/'outputs'/f'{tag}_final.json'
    if final.exists():return 'already_final'
    steps=[]
    if not guide.exists() or json.loads(guide.read_text())['status']=='running':steps.append('contour3d.py')
    steps.append('finalize.py')
    for script in steps:
        if script=='finalize.py':
            d=json.loads(guide.read_text())
            if d['status'] not in ['closed_verified','open_at_domain_boundaries']:return 'unresolved:'+d['status']
        command=[sys.executable,str(ROOT/'work'/script),case,str(float(f)),component]
        if script=='contour3d.py' and guide.exists():command+=['--resume']
        log=ROOT/'work'/f'{tag}_{script[:-3]}.log'
        with log.open('a') as out:r=subprocess.run(command,env=ENV,stdout=out,stderr=subprocess.STDOUT)
        print('QUEUE_STEP',tag,script,r.returncode,flush=True)
        if r.returncode:return 'failed:'+script
    return 'final'

def main(case,worker):
    assert case in ['s30','s45','s50','s55','s60','s65','s70','s73']
    lockdir=ROOT/'work/job_locks';lockdir.mkdir(exist_ok=True);seen=set();audit=[]
    while True:
        data=json.loads((ROOT/'outputs'/f'{case}_neutral_seed_scan.json').read_text())
        jobs=sorted(set((float(r['f_kHz']),r.get('component','S')) for r in data['records'] if r['side']=='lower'),key=lambda t:(abs(t[0]),t[0],t[1]))
        did=False;busy=False
        for f,component in jobs:
            tag=f'{case}_f{f}kHz_{component}'
            if tag in seen or (ROOT/'outputs'/f'{tag}_final.json').exists():continue
            with (lockdir/(tag+'.lock')).open('a') as lock:
                try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
                except BlockingIOError:busy=True;continue
                if active_curve(case,f,component):busy=True;continue
                status=run_one(case,f,component);seen.add(tag);audit.append(dict(tag=tag,status=status,time=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()));did=True
                (ROOT/'outputs'/f'{case}_queue_worker{worker}.json').write_text(json.dumps(audit,indent=2)+'\n')
        if not did:
            if seed_helper(case) or busy:time.sleep(3);continue
            break
    print('QUEUE_FINISHED',case,worker,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case');p.add_argument('worker');a=p.parse_args();main(a.case,a.worker)
