"""Migrate verified neutral components in frequency, allowing x/alpha/beta to move.

Intermediate operators guide identification only. Every candidate is landed on
an actual native station, then independently checked by the full finalizer.
"""
import argparse,fcntl,json,os,sys,time
from pathlib import Path
import numpy as np
from common import ROOT,CASES,snapshot
from contour import Family
from contour3d import bordered,SCALE
from core import mac
from run_batch import main as recover,run_one,seed_helper
FACTOR=2*np.pi*1000*.035/1810.6307411853716
ALLOWED=['s30','s45','s50','s55','s60','s65','s70','s73']

def correct(family,pred,plane,q,p,maxit=24):
    z=pred.copy()
    for it in range(maxit):
        C,dC,_,_=family.pencil(float(z[0]*SCALE[0]));e,J,t,v,left,res,lr=bordered(C,dC,z,q,p);F=np.r_[e.real,e.imag,np.dot(z-pred,plane)]
        if abs(e)<1e-12 and abs(F[2])<1e-9 and res<1e-8:return z,v,left,t,res,it
        A=np.vstack([J,plane]);weights=1/np.maximum(np.linalg.norm(A,axis=1),1e-30);dz=np.linalg.solve(A*weights[:,None],-F*weights)
        if np.linalg.norm(dz)>.1:dz*=.1/np.linalg.norm(dz)
        z+=dz
    raise RuntimeError('Frequency-surface corrector iteration limit')

def migrate(case,rec,w0,target):
    z=np.array(rec['z']);q=np.load(rec['vector_file']);p=np.load(rec['left_vector_file']);w=float(w0);step=.3;trace=[];failures=[]
    for loop in range(1200):
        if abs(w-target)<1e-11:return dict(finished=True,history=trace,failures=failures),z,q,p
        family=Family(case,w);C,dC,_,_=family.pencil(float(z[0]*SCALE[0]));_,_,plane,q,p,_,_=bordered(C,dC,z,q,p);new=w+np.sign(target-w)*min(step,abs(target-w));nextfamily=Family(case,new)
        try:
            zz,v,left,t,res,it=correct(nextfamily,z,plane,q,p);overlap=mac(q,v);assert overlap>.97 and np.linalg.norm(zz-z)<.25
        except Exception as e:
            step/=2;failures.append(dict(omega_reached=w,step=step,reason=repr(e)))
            if step<1e-5:return dict(finished=False,history=trace,failures=failures),z,q,p
            continue
        w=new;z,q,p=zz,v,left;trace.append(dict(omega_nd=w,z=z.tolist(),raw_MAC=overlap,right_residual=res));
        if it<=5:step=min(.75,step*1.2)
    return dict(finished=False,history=trace,failures=failures,reason='iteration_limit'),z,q,p

def native_land(case,f,z,q,p):
    family=Family(case,f*FACTOR);xi=float(z[0]*SCALE[0]);tries=[]
    for i in sorted(set([int(np.floor(xi)),int(np.ceil(xi))]),key=lambda i:abs(i-xi)):
        if not 1<=i<=family.count:continue
        pred=z.copy();pred[0]=i/SCALE[0]
        try:
            zz,v,left,t,res,it=correct(family,pred,np.array([1.,0.,0.]),q,p,35);assert mac(q,v)>.95
            return dict(station=i,alpha_r=float(zz[1]*SCALE[1]),beta_r=float(zz[2]*SCALE[2]),right_residual=res),v,left
        except Exception as e:tries.append(dict(station=i,reason=repr(e)))
    return dict(unresolved_native_landing=tries),None,None

def matches(case,f,r,q):
    for path in (ROOT/'outputs').glob(f'{case}_f{float(f)}kHz_*_final.json'):
        d=json.loads(path.read_text())
        for v in d['points']:
            if v['kind']!=0:continue
            if abs(v['xi']-r['station'])<1e-8 and abs(v['alpha_r']-r['alpha_r'])<1e-6 and abs(v['beta_r']-r['beta_r'])<1e-6 and mac(np.load(v['vector_file']),q)>1-1e-7:return path.name
    # Do not duplicate a previously unresolved component with an already closed guide.
    for path in (ROOT/'outputs').glob(f'{case}_f{float(f)}kHz_*_v3_guide.json'):
        if path.with_name(path.name.replace('_v3_guide','_final')).exists():continue
        d=json.loads(path.read_text())
        for rec in d['records']:
            v=rec['point']
            if abs(v['xi']-r['station'])<1e-8 and abs(v['alpha_r']-r['alpha_r'])<1e-6 and abs(v['beta_r']-r['beta_r'])<1e-6 and mac(np.load(rec['vector_file']),q)>1-1e-7:return 'unresolved:'+path.name
    return None

def add(case,f,r,q,audit_path):
    path=ROOT/'outputs'/f'{case}_neutral_seed_scan.json';data=json.loads(path.read_text());components=[v.get('component','S') for v in data['records'] if v['side']=='lower' and v['f_kHz']==f];n=1
    while 'F'+str(n) in components:n+=1
    component='F'+str(n);tag=f'{case}_f{float(f)}kHz_{component}';lockdir=ROOT/'work/job_locks';lockdir.mkdir(exist_ok=True)
    with (lockdir/(tag+'.lock')).open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX);vf=ROOT/'outputs'/f'{tag}_surface_seed.npy';np.save(vf,q)
        data['records'].append(dict(case=case,f_kHz=float(f),component=component,station=r['station'],alpha_r=r['alpha_r'],alpha_i=0.,beta_r=r['beta_r'],beta_i=0.,omega_r=f*FACTOR,omega_i=0.,growth=0.,side='lower',kind='native_seed',vector_file=str(vf),frequency_surface_audit=str(audit_path),native_landing_residual=r['right_residual']))
        tmp=path.with_suffix('.surface.tmp');tmp.write_text(json.dumps(data,indent=2)+'\n');tmp.replace(path);status=run_one(case,f,component)
    return component,status

def main(case):
    assert case in ALLOWED
    while seed_helper(case):time.sleep(3)
    recover(case,'surface_recovery');audit_path=ROOT/'outputs'/f'{case}_frequency_surface_search.json';audit=[]
    for sequence in [[-50.,-75.,-100.,-110.,-112.,-150.,-200.,-237.14,-240.2],[50.,75.,100.,125.,150.,175.,200.,230.,250.,264.53,275.,290.,291.6]]:
        for previous,f in zip(sequence,sequence[1:]):
            sources=[]
            for path in (ROOT/'outputs').glob(f'{case}_f{float(previous)}kHz_*_final.json'):
                d=json.loads(path.read_text())
                if d['status']!='closed_verified':continue
                guide=path.with_name(path.name.replace('_final','_v3_guide'));g=json.loads(guide.read_text());sources.append((guide,g))
            added=[]
            for source,guide in sources:
                records=guide['records'][:-1];indices=sorted(set(int(v) for v in np.linspace(0,len(records)-1,5)[:-1]))
                for index in indices:
                    rec=records[index];trial,z,q,p=migrate(case,rec,previous*FACTOR,f*FACTOR);trial.update(case=case,source_guide=str(source),source_index=index,source_f_kHz=previous,target_f_kHz=f)
                    if trial['finished']:
                        landing,v,left=native_land(case,f,z,q,p);trial['native_landing']=landing
                        if v is not None:
                            match=matches(case,f,landing,v);trial['matches_existing_component']=match
                            if match is None:
                                component,status=add(case,f,landing,v,audit_path);trial.update(new_component=component,run_status=status);added.append((component,status));print('SURFACE_COMPONENT',case,f,component,status,flush=True)
                    audit.append(trial);tmp=audit_path.with_suffix('.tmp');tmp.write_text(json.dumps(dict(case=case,attempts=audit,interpretation='Frequency continuation of previously verified components, with free x/alpha/beta and actual-native landing. Failure or inability to land is unresolved, never global stability.'),indent=2)+'\n');tmp.replace(audit_path)
            print('SURFACE_FREQUENCY',case,previous,f,'sources',len(sources),'new',added,flush=True)
    snapshot();print('SURFACE_SEARCH_FINISHED',case,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case');a=p.parse_args();main(a.case)
