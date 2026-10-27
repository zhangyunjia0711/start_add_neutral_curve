"""Restore an original station exposed by actual-profile turning-point correction."""
import json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/neutral_xbeta_extended_frequencies_20261027')
sys.path.insert(0,str(ROOT/'work'))
from common import CASES,snapshot,sha
from solver import Pencil,original,mac
from physical import Pair
def main():
    case='rough';f=-112;i=24;tag=f'{case}_f{f}kHz_S';path=ROOT/'outputs'/f'{tag}_final.json';d=json.loads(path.read_text());assert not d.get('native_crossing_repairs')
    tip=max(d['physical_tips'],key=lambda p:p['xi']);assert i<tip['xi']<i+1 and i-1<tip['guide_xi']<i
    omega=tip['omega_r'];p=Pencil(case,i,omega);qtip=np.load(tip['vector_file']);a0=complex(tip['alpha_r'],tip['alpha_i']);b0=tip['beta_r'];scan=[]
    def evaluate(b):
        a,q,r=p.solve_alpha(float(b),a0,qtip,True);overlap=mac(qtip,q);assert overlap>.95
        scan.append(dict(beta=float(b),alpha_r=a.real,alpha_i=a.imag,MAC=overlap));return a,q,r
    bs=np.linspace(b0-2,b0+2,161);vs=[evaluate(b)[0].imag for b in bs];brackets=[(l,r) for l,r,gl,gr in zip(bs,bs[1:],vs,vs[1:]) if gl*gr<0];assert len(brackets)==2,brackets
    roots=[]
    for l,r in brackets:
        b=brentq(lambda b:evaluate(b)[0].imag,l,r,xtol=1e-11);a,q,record=evaluate(b);assert abs(a.imag)<1e-8;roots.append((b,a,q,record))
    before=sorted([p for p in d['points'] if p['kind']==0 and p['xi']==i-1],key=lambda p:p['beta_r']);assert len(before)==2
    pair=Pair(case,i-1,omega,tag+f'_repair{i}');traces=[];new=[]
    for side,start,end in zip(['lower','upper'],before,roots):
        b=start['beta_r'];a=complex(start['alpha_r'],start['alpha_i']);q=np.load(start['vector_file']);trace=[]
        for theta in np.linspace(0,1,21):
            pp=pair.pencil(float(theta))
            for iteration in range(25):
                aa,v,r=pp.solve_alpha(b,a,q,True);overlap=mac(q,v);assert overlap>.95;a,q=aa,v
                if abs(a.imag)<1e-11:break
                b+=np.clip(-a.imag/r['alpha_beta_i'],-.15,.15)
            assert abs(a.imag)<1e-8
            trace.append(dict(theta=float(theta),alpha_r=a.real,alpha_i=a.imag,beta=b,adjacent_MAC=overlap))
        assert abs(b-end[0])<1e-7 and abs(a-end[1])<1e-7 and mac(q,end[2])>1-1e-7
        b,a,q,r=end;state=dict(omega_r=omega,beta_r=b,alpha_r=a.real,alpha_i=a.imag);actual,v,receipt=original(case,i,state,b,f'{tag}_added_native{i}_{side}')
        assert abs(actual['growth'])<1e-8 and abs(complex(actual['alpha_r'],actual['alpha_i'])-a)<1e-7 and mac(q,v)>1-1e-7
        qpath=ROOT/'outputs'/f'{tag}_native{i}_{side}.npy';np.save(qpath,q)
        row=dict(case=case,kind=0,xi=float(i),x_mm=CASES[case]['x_mm'][i-1],x_flat_mm=CASES[case]['x_mm'][i-1]-CASES[case]['flat_start_mm'],alpha_r=a.real,alpha_i=a.imag,beta_r=b,beta_i=0.,omega_r=omega,omega_i=0.,growth=-a.imag,right_residual=r['right_residual'],left_residual=r['left_residual'],raw_MAC=r['raw_MAC'],guide_order=(start['guide_order']+tip['guide_order'])/2,guide_index=-1,alpha_beta_i=r['alpha_beta_i'],vector_file=str(qpath),original_receipt=receipt,physical_added_native_station=True)
        new.append(row);traces.append(dict(side=side,records=trace,original_check=receipt))
    previous=path.with_name(path.stem+'_before_native_repair.json');assert not previous.exists();previous.write_bytes(path.read_bytes())
    points=sorted([p for p in d['points'] if p['kind']!=2]+new,key=lambda p:p['guide_order']);points.append(dict(points[0],kind=2,verified_closure=True));d['points']=points;d['native_point_count']+=2;d['native_station_count']+=1
    d['native_crossing_repairs']=[dict(station=i,reason='Actual-profile maximum crosses an additional native station',before_file=str(previous),before_sha256=sha(previous),scanning=scan,traces=traces)]
    d['maximum_growth_error']=max(abs(p['growth']) for p in points);d['maximum_right_residual']=max(p['right_residual'] for p in points);path.write_text(json.dumps(d,indent=2)+'\n');snapshot()
    print('REPAIRED',case,f,[(p['xi'],p['beta_r']) for p in new],flush=True)
if __name__=='__main__':main()
