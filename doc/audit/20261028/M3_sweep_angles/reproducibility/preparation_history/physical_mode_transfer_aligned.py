"""First remap the same physical source profile, then blend at common physical n."""
import configparser,hashlib,json,sys
from pathlib import Path
import h5py,numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.sparse.linalg import norm
from physical_mode_transfer import ROOT,primitive
sys.path.insert(0,str(ROOT/'work'))
from solver import Pencil,export
from core import qep,mac

def matrix_profile(source,target,ss,ts,stage,lam,n,q):
    token=hashlib.sha256((stage+float(lam).hex()).encode()).hexdigest()[:12];d=ROOT/'inputs'/f'mode_aligned_{source}{ss}_{target}{ts}_{token}'
    if not (d/'generation.json').exists():
        assert np.all(np.diff(n)>0) and q[:,0].min()>0 and q[:,4].min()>0
        (d/'in').mkdir(parents=True,exist_ok=True);grid=np.zeros((400,9,2,2));grid[:,:,0,0]=.0001*np.arange(-4,5)[None,:]+1e-8*n[:,None];grid[:,:,1,0]=n[:,None];flow=np.zeros((400,9,6,2));flow[:,:,:,0]=q[:,None,:]
        for name,key,value in [('grid_vertical2d_N.h5','grid',grid),('flow_vertical2d_N.h5','baseflow',flow)]:
            with h5py.File(d/'in'/name,'w') as h:h.create_dataset(key,data=value)
        cfg=configparser.ConfigParser();cfg.read(ROOT/'cases'/target/'config.ini');cfg['domain'].update({'in':'9','jn':'400','iloc':'5','istart':'1','iend':'9'})
        with (d/'config.ini').open('w') as f:cfg.write(f)
        (d/'generation.json').write_text(json.dumps(dict(source=source,target=target,stage=stage,lambda_value=lam,identification_only=True),indent=2)+'\n')
    K,receipt=export(d,5,-7.,0.,label=f'mode_aligned_{source}_{target}');return K,receipt

def transfer(source,target):
    rows=json.loads((ROOT/'outputs/angle_mode_pilot.json').read_text());s=next(r for r in rows if r['case']==source);t=next(r for r in rows if r['case']==target);ss=s['station'];ts=t['station'];na,qa=primitive(source,ss);nb,qb=primitive(target,ts);interpolate=PchipInterpolator(na,qa,axis=0,extrapolate=False)
    def sample(n):return interpolate(np.clip(n,na[0],na[-1]))
    a=complex(s['alpha_r'],s['alpha_i']);q=np.load(s['vector_file']);oldn=na;trace=[];failures=[]
    K,receipt=matrix_profile(source,target,ss,ts,'grid',0.,na,qa);native=Pencil(source,ss,0.).spatial(-7.);start_errors=[float(norm(x-y))/max(1.,float(norm(y))) for x,y in zip(K,native)];assert max(start_errors)<2e-12,start_errors
    a,q,_,_=qep(K,a,q,True)
    for stage in ['grid','flow']:
        lam=0.;step=.05
        for it in range(600):
            if lam>=1-1e-14:break
            new=min(1.,lam+step);n=(1-new)*na+new*nb if stage=='grid' else nb
            profile=sample(n) if stage=='grid' else (1-new)*sample(nb)+new*qb
            K,receipt=matrix_profile(source,target,ss,ts,stage,new,n,profile);v=q.reshape((-1,5));ref=np.column_stack([np.interp(n,oldn,v[:,k].real,right=0)+1j*np.interp(n,oldn,v[:,k].imag,right=0) for k in range(5)]).ravel()
            try:aa,qq,r,_=qep(K,a,ref,True);overlap=mac(ref,qq);assert overlap>.98
            except Exception as e:
                failures.append(dict(stage=stage,lambda_value=new,step=step,reason=repr(e)));step/=2
                if step<1e-6:raise RuntimeError('Aligned mode transfer unresolved: '+repr(failures[-1]))
                continue
            lam=new;a,q,oldn=aa,qq,n;trace.append(dict(stage=stage,lambda_value=lam,alpha_r=a.real,alpha_i=a.imag,physical_adjacent_MAC=overlap,right_residual=r['right_residual'],export=receipt));step=min(.1,step*1.2)
        assert lam>=1-1e-14
        print('ALIGNED_STAGE',source,target,stage,a,'steps',len(trace),flush=True)
    native=Pencil(target,ts,0.).spatial(-7.);end_errors=[float(norm(x-y))/max(1.,float(norm(y))) for x,y in zip(K,native)];assert max(end_errors)<2e-12,end_errors
    aa,qq,r=Pencil(target,ts,0.).solve_alpha(-7.,a,q,True);assert mac(q,qq)>1-1e-7;path=ROOT/'outputs'/f'{source}_to_{target}_aligned_mode_seed.npy';np.save(path,qq)
    result=dict(source_case=source,target_case=target,source_station=ss,target_station=ts,alpha_r=aa.real,alpha_i=aa.imag,vector_file=str(path),start_matrix_errors=start_errors,end_matrix_errors=end_errors,trace=trace,failures=failures,native_receipt=r,identification_only=True,definition='Stage 1 remaps the same physical source primitive profiles onto the target normal grid using PCHIP. Stage 2 blends primitive profiles at the same physical normal coordinates. Intermediate profiles are identification aids, not reported CFD or neutral-curve data.')
    (ROOT/'outputs'/f'{source}_to_{target}_aligned_mode_transfer.json').write_text(json.dumps(result,indent=2)+'\n');print('ALIGNED_MODE',source,target,aa,'steps',len(trace),flush=True);return result
if __name__=='__main__':
    for source in ['s65','s73']:transfer(source,'s70')
