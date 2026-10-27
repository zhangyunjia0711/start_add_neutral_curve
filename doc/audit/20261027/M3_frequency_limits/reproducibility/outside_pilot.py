"""Direct signed-frequency/beta mode continuation outside proposed limits.

This independently checks representative stations where a beta-fold projection
cannot be continued. Failed continuations are reported, never counted stable.
"""
import concurrent.futures,json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/neutral_m3_frequency_limits_20261027')
sys.path.insert(0,str(ROOT/'work'))
from common import CASES,sha
from solver import Pencil,mac
from fold_core import W
PREVIOUS=ROOT.parent/'neutral_xbeta_omega0_20261026/outputs'
def one(job):
    case,i=job[:2];sign=job[2] if len(job)>2 else -1;label='positive' if sign>0 else 'negative';limit=json.loads((ROOT/'outputs'/f'{case}_{label}_refined_limit.json').read_text())['critical'];target=limit['f_kHz']+sign*.1;d=json.loads((PREVIOUS/f'{case}_final.json').read_text());r=min([p for p in d['points'] if p['kind']==0 and p['xi']==i],key=lambda p:p['beta_r']);a=complex(r['alpha_r'],r['alpha_i']);b=r['beta_r'];q=np.load(r['vector_file']);f=0.;step=1.;freq=[];failures=[]
    while abs(f-target)>1e-9:
        df=sign*min(step,abs(target-f));p=Pencil(case,i,f*W);a,q,r=p.solve_alpha(b,a,q,True);predict=a+df*W*complex(r['alpha_omega_r'],r['alpha_omega_i'])
        try:aa,v,rr=Pencil(case,i,(f+df)*W).solve_alpha(b,predict,q,True);assert mac(q,v)>.97
        except Exception as e:
            step/=2;failures.append(dict(f_kHz=f,step=step,reason=repr(e)))
            if step<1e-4:break
            continue
        f+=df;a,q,r=aa,v,rr;freq.append(dict(f_kHz=f,alpha_r=a.real,alpha_i=a.imag,MAC=rr['raw_MAC']));step=min(1.,step*1.2)
    state=dict(case=case,station=i,x_mm=CASES[case]['x_mm'][i-1],target_f_kHz=target,reached_f_kHz=f,frequency_path=freq,failures=failures,beta_scans=[],status='running')
    if abs(f-target)>1e-8:state['status']='frequency_unresolved'
    else:
        p=Pencil(case,i,f*W);state['beta_scans'].append(dict(beta=b,alpha_r=a.real,alpha_i=a.imag,growth=-a.imag,MAC=1.))
        for direction,end in [(-1,-64.),(1,64.)]:
            bb=b;aa=a;qq=q.copy();delta=1.;complete=True
            while direction*(end-bb)>1e-9:
                db=direction*min(delta,abs(end-bb))
                try:
                    _,_,rr=p.solve_alpha(bb,aa,qq,True);predict=aa+db*complex(rr['alpha_beta_r'],rr['alpha_beta_i']);at,qt,rt=p.solve_alpha(bb+db,predict,qq,True);assert mac(qq,qt)>.97
                except Exception as e:
                    delta/=2;state['failures'].append(dict(beta=bb,db=db,reason=repr(e)))
                    if delta<1e-4:complete=False;break
                    continue
                if (direction>0 and rr['alpha_beta_i']<0<rt['alpha_beta_i']) or (direction<0 and rt['alpha_beta_i']<0<rr['alpha_beta_i']):
                    anchor_b,anchor_a,anchor_q,anchor_der=bb,aa,qq.copy(),complex(rr['alpha_beta_r'],rr['alpha_beta_i'])
                    def derivative(beta):return p.solve_alpha(float(beta),anchor_a+(beta-anchor_b)*anchor_der,anchor_q,True)[2]['alpha_beta_i']
                    bp=brentq(derivative,min(bb,bb+db),max(bb,bb+db),xtol=1e-9);ap,vp,rp=p.solve_alpha(bp,anchor_a+(bp-anchor_b)*anchor_der,anchor_q,True);state['beta_scans'].append(dict(beta=bp,alpha_r=ap.real,alpha_i=ap.imag,growth=-ap.imag,MAC=rp['raw_MAC'],refined_local_maximum=True))
                bb+=db;aa,qq=at,qt;state['beta_scans'].append(dict(beta=bb,alpha_r=aa.real,alpha_i=aa.imag,growth=-aa.imag,MAC=rt['raw_MAC'],alpha_omega_r=rt['alpha_omega_r'],alpha_beta_i=rt['alpha_beta_i']));delta=min(1.,delta*1.2)
            state.setdefault('directions',[]).append(dict(direction=direction,end_beta=bb,complete=complete))
        state['status']='complete_beta_window' if all(d['complete'] for d in state['directions']) else 'beta_projection_unresolved'
        state['maximum_growth']=max(state['beta_scans'],key=lambda p:p['growth'])
    path=ROOT/'outputs'/f'{case}_outside_{label}_i{i:04d}.json';path.write_text(json.dumps(state,indent=2)+'\n');print('OUTSIDE',label,case,i,state['status'],state.get('maximum_growth'),flush=True);return state
def main():
    jobs=[('bluntplate',i) for i in [360,450,550]]+[('rough',i) for i in [90,140,200]]
    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as pool:results=list(pool.map(one,jobs))
    (ROOT/'outputs/outside_pilot_summary.json').write_text(json.dumps(results,indent=2)+'\n')
if __name__=='__main__':main()
