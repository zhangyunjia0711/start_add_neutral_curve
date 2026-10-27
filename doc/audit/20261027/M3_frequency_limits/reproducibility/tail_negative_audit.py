"""Independent downstream seed for the negative-frequency merger branch.

The first upstream continuation can turn in x or lose a beta projection. A
stationary downstream seed supplies a separate continuation for missing stations.
"""
import concurrent.futures,json,sys
from pathlib import Path
import numpy as np
ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/neutral_m3_frequency_limits_20261027')
sys.path.insert(0,str(ROOT/'work'))
from common import CASES,sha
from solver import Pencil,qep,mac
from fold_core import Equations,S,W
from pilot import fold
from trace_folds import solve,blended
PREVIOUS=ROOT.parent/'neutral_xbeta_omega0_20261026/outputs'
def one(case):
    n=CASES[case]['stations'];d=json.loads((PREVIOUS/f'{case}_final.json').read_text());rows=sorted([p for p in d['points'] if p['kind']==0 and int(p['xi'])==n],key=lambda p:p['beta_r']);assert len(rows)==2
    r=rows[0];q=np.load(r['vector_file']);z=np.array([r['alpha_r'],r['beta_r']])/S[:2];f=0.;step=.05;trace=[];fails=[]
    for k in range(1000):
        p=Pencil(case,n,f*W);_,q,_,left=qep(p.spatial(z[1]*10),z[0]*20,q,True);eq=Equations(p,q,left);r=eq.evaluate(np.r_[z,f/100]);df=-step
        try:
            trial=z+np.linalg.solve(r['J'][:2,:2],-r['J'][:2,2]*df/100);new=Equations(Pencil(case,n,(f+df)*W),q,left)
            for it in range(20):
                rr=new.evaluate(np.r_[trial,(f+df)/100]);F=rr['F'][:2]
                if np.linalg.norm(F)<1e-12 and rr['right_residual']<1e-8:break
                delta=np.linalg.solve(rr['J'][:2,:2],-F)
                if np.linalg.norm(delta)>.1:delta*=.1/np.linalg.norm(delta)
                trial+=delta
            assert np.linalg.norm(F)<1e-12 and rr['right_residual']<1e-8 and rr['MAC']>.97 and np.linalg.norm(trial-z)<.5
        except Exception as e:
            fails.append(dict(f_kHz=f,df=df,reason=repr(e)));step/=2
            if step<.0001:break
            continue
        f+=df;z=trial;q=rr['q'];trace.append(dict(f_kHz=f,alpha_r=float(z[0]*20),beta_r=float(z[1]*10),MAC=rr['MAC']))
        if it<=3:step=min(.5,step*1.2)
    seed=fold(case,n,float(z[0]*20),float(z[1]*10),f,q,case+'_negative_tail');assert seed['status']=='solved'
    state=dict(case=case,seed=seed,frequency_continuation=trace,frequency_failures=fails,records=[],failures=[],status='running');path=ROOT/'outputs'/f'{case}_negative_tail_fold_curve.json';wave=ROOT/'work'/f'{case}_negative_tail_vectors';wave.mkdir(exist_ok=True);q=np.load(seed['vector_file']);z=np.array([seed['alpha_r'],seed['beta_r'],seed['f_kHz']])/S
    existing=json.loads((ROOT/'outputs'/f'{case}_negative_native_fold_curve.json').read_text());front={p['station']:p for p in existing['records']};lower=min(front);previous=n
    for i in range(n,lower-1,-1):
        p=Pencil(case,i,z[2]*100*W)
        try:row,v,zz=solve(p,z,q)
        except Exception as e:
            state['failures'].append(dict(from_station=previous,to_station=i,reason=repr(e)));ok=False
            a=Pencil(case,previous,z[2]*100*W)
            for div in [4,16]:
                try:
                    zt,qt=z.copy(),q.copy()
                    for t in np.linspace(0,1,div+1)[1:]:row,qt,zt=solve(blended(a,p,float(t),zt[2]*100),zt,qt)
                    row,v,zz=solve(p,zt,qt);ok=True;break
                except Exception as ee:state['failures'].append(dict(to_station=i,subdivisions=div,reason=repr(ee)))
            if not ok:state['status']='stopped_with_failure';break
        qp=wave/f'i{i:04d}.npy';np.save(qp,v);row.update(station=i,xi=float(i),x_mm=CASES[case]['x_mm'][i-1],vector_file=str(qp),vector_sha256=sha(qp));state['records'].append(row);z,q,previous=zz,v,i
        if len(state['records'])%10==0:path.write_text(json.dumps(state,indent=2)+'\n')
        if i in front and abs(row['f_kHz']-front[i]['f_kHz'])<1e-5 and abs(row['beta_r']-front[i]['beta_r'])<1e-5:
            state.update(status='matched_upstream_fold_curve',matched_station=i);break
        if len(state['records'])%100==0:print('TAIL_TRACE',case,i,row['f_kHz'],flush=True)
    state['records'].sort(key=lambda p:p['station']);state['native_extreme']=min(state['records'],key=lambda p:p['f_kHz']);path.write_text(json.dumps(state,indent=2)+'\n');print('TAIL_FINAL',case,state['status'],len(state['records']),previous,state['native_extreme']['f_kHz'],flush=True);return state
def main():
    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as pool:results=list(pool.map(one,CASES))
    (ROOT/'outputs/tail_negative_summary.json').write_text(json.dumps(results,indent=2)+'\n')
if __name__=='__main__':main()
