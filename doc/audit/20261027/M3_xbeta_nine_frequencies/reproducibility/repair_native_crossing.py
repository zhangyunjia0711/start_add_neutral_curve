"""Restore native crossings exposed by physical-profile fold correction."""
import json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/neutral_xbeta_extra_frequencies_20261027')
sys.path.insert(0,str(ROOT/'work'))
from common import CASES,sha,snapshot
from solver import Pencil,original,mac
from physical import Pair

def main():
    case='rough';f=-75;tag=f'{case}_f{f}kHz_S';path=ROOT/'outputs'/f'{tag}_final.json';d=json.loads(path.read_text())
    assert not d.get('native_crossing_repairs'),'Do not duplicate a repair'
    tip=d['physical_tips'][0];assert 41<tip['xi']<42 and 40<tip['guide_xi']<41
    i=41;omega=tip['omega_r'];p=Pencil(case,i,omega);qtip=np.load(tip['vector_file']);a0=complex(tip['alpha_r'],tip['alpha_i']);b0=tip['beta_r']
    evals=[]
    def evaluate(beta):
        a,q,r=p.solve_alpha(float(beta),a0,qtip,True);assert mac(qtip,q)>.95
        evals.append(dict(beta=float(beta),alpha_r=a.real,alpha_i=a.imag,MAC=mac(qtip,q)))
        return a,q,r
    edges=np.linspace(b0-2,b0+2,81);vals=[evaluate(b)[0].imag for b in edges]
    brackets=[(l,r) for l,r,gl,gr in zip(edges,edges[1:],vals,vals[1:]) if gl*gr<0];assert len(brackets)==2,brackets
    roots=[]
    for lo,hi in brackets:
        beta=brentq(lambda b:evaluate(b)[0].imag,lo,hi,xtol=1e-11);a,q,r=evaluate(beta);assert abs(a.imag)<1e-8
        roots.append((beta,a,q,r))
    # Follow each native-40 neutral root through the actual local profiles to native 41.
    before=sorted([p for p in d['points'] if p['kind']==0 and p['xi']==40],key=lambda p:p['beta_r']);assert len(before)==2
    pair=Pair(case,40,omega,tag+'_repair41');traces=[];new=[]
    for side,start,end in zip(['lower','upper'],before,roots):
        beta=start['beta_r'];a=complex(start['alpha_r'],start['alpha_i']);q=np.load(start['vector_file']);trace=[]
        for theta in np.linspace(0,1,21):
            pp=pair.pencil(float(theta))
            for iteration in range(20):
                aa,v,r=pp.solve_alpha(beta,a,q,True);assert mac(q,v)>.95;a,q=aa,v
                if abs(a.imag)<1e-11:break
                step=np.clip(-a.imag/r['alpha_beta_i'],-.2,.2);beta+=step
            assert abs(a.imag)<1e-8
            trace.append(dict(theta=float(theta),beta=beta,alpha_r=a.real,alpha_i=a.imag,vector_MAC_to_final=mac(q,end[2])))
        assert abs(beta-end[0])<1e-7 and abs(a-end[1])<1e-7 and mac(q,end[2])>1-1e-7
        beta,a,q,r=end;state=dict(omega_r=omega,beta_r=beta,alpha_r=a.real,alpha_i=a.imag)
        actual,v,receipt=original(case,i,state,beta,f'{tag}_added_native41_{side}')
        assert abs(actual['growth'])<1e-8 and abs(complex(actual['alpha_r'],actual['alpha_i'])-a)<1e-7 and mac(v,q)>1-1e-7
        qpath=ROOT/'outputs'/f'{tag}_native41_{side}.npy';np.save(qpath,q)
        order=(start['guide_order']+tip['guide_order'])/2
        row=dict(case=case,kind=0,xi=41.,x_mm=CASES[case]['x_mm'][40],x_flat_mm=CASES[case]['x_mm'][40]-CASES[case]['flat_start_mm'],
            alpha_r=a.real,alpha_i=a.imag,beta_r=beta,beta_i=0.,omega_r=omega,omega_i=0.,growth=-a.imag,
            right_residual=r['right_residual'],left_residual=r['left_residual'],raw_MAC=r['raw_MAC'],guide_order=order,guide_index=-1,
            alpha_beta_i=r['alpha_beta_i'],vector_file=str(qpath),physical_added_native_station=True,original_receipt=receipt)
        new.append(row);traces.append(dict(side=side,records=trace,original_check=receipt))
    previous=path.with_name(path.stem+'_before_native_repair.json');assert not previous.exists();previous.write_bytes(path.read_bytes())
    points=sorted([r for r in d['points'] if r['kind']!=2]+new,key=lambda p:p['guide_order']);points.append(dict(points[0],kind=2,verified_closure=True))
    d['points']=points;d['native_point_count']+=2;d['native_station_count']+=1
    d['native_crossing_repairs']=[dict(tip_index=0,station=41,reason='Actual-profile x maximum moved from bracket40 to bracket41; add both previously omitted native roots',before_file=str(previous),before_sha256=sha(previous),scanning=evals,traces=traces)]
    d['maximum_growth_error']=max(abs(r['growth']) for r in points);d['maximum_right_residual']=max(r['right_residual'] for r in points)
    path.write_text(json.dumps(d,indent=2)+'\n');snapshot();print('REPAIRED_NATIVE',[(r['xi'],r['beta_r'],r['alpha_r']) for r in new],flush=True)
if __name__=='__main__':main()
