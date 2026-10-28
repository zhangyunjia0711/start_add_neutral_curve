"""Identify modes along actual local primitive-profile/normal-grid homotopies.

Synthetic intermediate profiles identify a branch only. Reported end states
solve the native target operator, whose matrix is independently compared.
"""
import configparser,hashlib,json,sys
from pathlib import Path
import h5py,numpy as np
from scipy.sparse.linalg import norm
ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/neutral_m3_sweep_consistent_20261028')
sys.path.insert(0,str(ROOT/'work'))
from common import CASES,sha
from solver import Pencil,export
from core import qep,mac
from pilot_additional import profile

def primitive(case,station):
    n,basis=profile(case,station)
    with h5py.File(ROOT/'cases'/case/'in/flow_vertical2d_N.h5') as h:raw=h['baseflow'][:,station-1,:,0]
    q=np.zeros((400,6));q[:,:raw.shape[1]]=raw;q[:,1:3]=raw[:,1:3]@basis.T
    if raw.shape[1]==5:q[:,5]=q[:,0]*q[:,4]/(1.4*36.)
    return n,q

def matrix_at(source,target,ss,ts,lam,na,qa,nb,qb):
    n=(1-lam)*na+lam*nb;q=(1-lam)*qa+lam*qb;token=hashlib.sha256(float(lam).hex().encode()).hexdigest()[:12]
    d=ROOT/'inputs'/f'mode_profiles_{source}{ss}_{target}{ts}_{token}'
    if not (d/'generation.json').exists():
        (d/'in').mkdir(parents=True,exist_ok=True);grid=np.zeros((400,9,2,2));grid[:,:,0,0]=.0001*np.arange(-4,5)[None,:]+1e-8*n[:,None];grid[:,:,1,0]=n[:,None];flow=np.zeros((400,9,6,2));flow[:,:,:,0]=q[:,None,:]
        for name,key,value in [('grid_vertical2d_N.h5','grid',grid),('flow_vertical2d_N.h5','baseflow',flow)]:
            with h5py.File(d/'in'/name,'w') as h:h.create_dataset(key,data=value)
        cfg=configparser.ConfigParser();cfg.read(ROOT/'cases'/target/'config.ini');cfg['domain'].update({'in':'9','jn':'400','iloc':'5','istart':'1','iend':'9'})
        with (d/'config.ini').open('w') as f:cfg.write(f)
        (d/'generation.json').write_text(json.dumps(dict(source_case=source,target_case=target,source_station=ss,target_station=ts,lambda_value=lam,identification_only=True,definition='Matched-j interpolation of local primitive profiles and monotone physical normal coordinates; consistent local parallel LST; not a new CFD result.'),indent=2)+'\n')
    K,receipt=export(d,5,-7.,0.,label=f'mode_profile_{source}_{target}')
    return n,K,receipt

def transfer(source,target):
    rows=json.loads((ROOT/'outputs/angle_mode_pilot.json').read_text());s=next(r for r in rows if r['case']==source);t=next(r for r in rows if r['case']==target);ss=s['station'];ts=t['station'];na,qa=primitive(source,ss);nb,qb=primitive(target,ts);a=complex(s['alpha_r'],s['alpha_i']);q=np.load(s['vector_file']);lam=0.;step=.05;trace=[];failures=[];oldn=na
    n,K,receipt=matrix_at(source,target,ss,ts,0.,na,qa,nb,qb);native=Pencil(source,ss,0.).spatial(-7.);start_errors=[float(norm(x-y))/max(1.,float(norm(y))) for x,y in zip(K,native)];assert max(start_errors)<2e-12,start_errors
    a,q,_,_=qep(K,a,q,True)
    for iteration in range(500):
        if lam>=1-1e-14:break
        new=min(1.,lam+step);n,K,receipt=matrix_at(source,target,ss,ts,new,na,qa,nb,qb);v=q.reshape((-1,5));ref=np.column_stack([np.interp(n,oldn,v[:,k].real,right=0)+1j*np.interp(n,oldn,v[:,k].imag,right=0) for k in range(5)]).ravel()
        try:
            aa,qq,r,_=qep(K,a,ref,True);overlap=mac(ref,qq);assert overlap>.98
        except Exception as e:
            failures.append(dict(lambda_value=new,step=step,reason=repr(e)));step/=2
            if step<1e-5:raise RuntimeError('Physical mode transfer unresolved: '+repr(failures[-1]))
            continue
        lam=new;a,q,oldn=aa,qq,n;trace.append(dict(lambda_value=lam,alpha_r=a.real,alpha_i=a.imag,physical_adjacent_MAC=overlap,right_residual=r['right_residual'],export=receipt));step=min(.2,step*1.25)
    assert lam>=1-1e-14
    native=Pencil(target,ts,0.).spatial(-7.);end_errors=[float(norm(x-y))/max(1.,float(norm(y))) for x,y in zip(K,native)];assert max(end_errors)<2e-12,end_errors
    aa,qq,r=Pencil(target,ts,0.).solve_alpha(-7.,a,q,True);assert mac(q,qq)>1-1e-7
    path=ROOT/'outputs'/f'{source}_to_{target}_physical_mode_seed.npy';np.save(path,qq);result=dict(source_case=source,target_case=target,source_station=ss,target_station=ts,alpha_r=aa.real,alpha_i=aa.imag,vector_file=str(path),start_matrix_errors=start_errors,end_matrix_errors=end_errors,trace=trace,failures=failures,native_receipt=r,previous_candidate_alpha_difference=abs(aa-complex(t['alpha_r'],t['alpha_i'])),previous_candidate_MAC=mac(qq,np.load(t['vector_file'])),identification_only=True)
    (ROOT/'outputs'/f'{source}_to_{target}_physical_mode_transfer.json').write_text(json.dumps(result,indent=2)+'\n');print('PHYSICAL_MODE',source,target,aa,'steps',len(trace),'previous_difference',result['previous_candidate_alpha_difference'],'MAC',result['previous_candidate_MAC'],flush=True)
    return result
if __name__=='__main__':
    for source,target in [('s50','s55'),('s60','s65'),('s65','s73'),('s65','s70'),('s73','s70')]:transfer(source,target)
