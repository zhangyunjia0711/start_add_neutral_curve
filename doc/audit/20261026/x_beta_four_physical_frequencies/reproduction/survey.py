"""Broad beta survey at both physical seed locations, with raw vector tracking."""
import concurrent.futures,json,time
import numpy as np
from common import ROOT,CASES,snapshot
from solver import Pencil,original,mac

def track_beta(case,omega,seed_beta,initial,qref,end,step=.5):
    p=Pencil(case,CASES[case]['seed_station'],omega);beta=seed_beta;a=initial;q=qref.copy();records=[];failures=[];roots=[]
    progress=ROOT/'work'/f'{case}_w{int(omega)}_seed{int(seed_beta)}_to{int(end)}_progress.json';started=time.monotonic()
    a,q,r=p.solve_alpha(beta,a,q,True);records.append(r)
    while abs(beta-end)>1e-10:
        db=np.sign(end-beta)*min(step,abs(end-beta));derivative=complex(r['alpha_beta_r'],r['alpha_beta_i']);nb=beta+db
        try:aa,qq,rr=p.solve_alpha(nb,a+db*derivative,q,True)
        except Exception as e:
            step/=2;failures.append(dict(beta=beta,step=step,reason=repr(e)))
            if step<1e-5:break
            continue
        if r['growth']*rr['growth']<0:
            low=(beta,a,q.copy(),r);high=(nb,aa,qq.copy(),rr)
            for it in range(50):
                b0,a0,q0,r0=low;b1,a1,q1,r1=high;theta=-r0['growth']/(r1['growth']-r0['growth']);theta=np.clip(theta,.05,.95);bb=b0+theta*(b1-b0)
                ref,tar=(q0,a0) if theta<.5 else (q1,a1);ar,qr,root=p.solve_alpha(bb,tar,ref,True)
                if abs(root['growth'])<1e-10:break
                if root['growth']*r0['growth']>0:low=(bb,ar,qr,root)
                else:high=(bb,ar,qr,root)
            assert abs(root['growth'])<1e-8
            path=ROOT/'outputs'/f'{case}_w{int(omega)}_seed{int(seed_beta)}_end{int(end)}_root{len(roots)}.npz';np.savez(path,q=qr)
            root.update(vector_file=str(path));roots.append(root)
        beta,a,q,r=nb,aa,qq,rr;records.append(r)
        if rr['raw_MAC']>.999:step=min(.75,step*1.1)
        if len(records)%50==0:
            progress.write_text(json.dumps(dict(case=case,omega=omega,seed_beta=seed_beta,end=end,beta=float(beta),points=len(records),max_growth=max(r['growth'] for r in records),seconds=time.monotonic()-started),indent=2)+'\n')
            print('SCAN',case,omega,seed_beta,end,len(records),float(beta),max(r['growth'] for r in records),flush=True)
        if len(records)>3000 or time.monotonic()-started>300:
            failures.append(dict(beta=float(beta),reason='Survey wall-time/point guard; unresolved range retained'));break
    return dict(case=case,omega=omega,seed_beta=seed_beta,end=end,reached_end=bool(abs(beta-end)<1e-8),records=records,failures=failures,roots=roots)

def one(item):
    case,omega,seed_beta,end=item;pilot=next(p for p in json.loads((ROOT/'outputs/pilot.json').read_text()) if p['case']==case)
    path=ROOT/'outputs'/f'{case}_w{int(omega)}_seed{int(seed_beta)}_to{int(end)}_survey.json'
    if path.exists():
        result=json.loads(path.read_text())
        return dict(path=str(path),case=case,omega=omega,seed_beta=seed_beta,end=end,reached_end=result['reached_end'],maximum_growth=max(r['growth'] for r in result['records']),root_count=len(result['roots']))
    track=next(t for t in pilot['tracks'] if t['omega']==omega and t['beta']==seed_beta);assert track['finished']
    r=track['records'][-1];a=complex(r['alpha_r'],r['alpha_i']);q=np.load(track['vector_file'])['q'];result=track_beta(case,omega,seed_beta,a,q,end)
    path.write_text(json.dumps(result,indent=2)+'\n')
    print('SURVEY',case,omega,seed_beta,end,result['reached_end'],'MAXG',max(r['growth'] for r in result['records']),'ROOTS',[(r['beta_r'],r['alpha_r']) for r in result['roots']],flush=True)
    return dict(path=str(path),case=case,omega=omega,seed_beta=seed_beta,end=end,reached_end=result['reached_end'],maximum_growth=max(r['growth'] for r in result['records']),root_count=len(result['roots']))

if __name__=='__main__':
    snapshot();jobs=[(case,w,b,end) for case in CASES for w in [50.,100.] for b in [-7.,7.] for end in [-200.,200.]]
    with concurrent.futures.ProcessPoolExecutor(max_workers=8) as pool:results=list(pool.map(one,jobs))
    (ROOT/'outputs/survey_summary.json').write_text(json.dumps(results,indent=2)+'\n')
