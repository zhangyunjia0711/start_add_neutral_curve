"""Correct native neutral points and actual profile-interpolated x turning points."""
import argparse,json,csv
import numpy as np
from scipy.sparse.linalg import norm
from common import ROOT,CASES,snapshot
from solver import Pencil,qep,mac,original
OMEGA=None
TAG=None
from physical import Pair

from tip_correction import solve_tip
def tip(case,left,right,index):
    return solve_tip(case,left,right,index,OMEGA,TAG)

def main():
    global OMEGA,TAG
    parser=argparse.ArgumentParser();parser.add_argument('case');parser.add_argument('frequency',type=int,choices=[-75,-25,25,75]);parser.add_argument('component',choices=['S']);args=parser.parse_args();case=args.case;TAG=f'{case}_f{args.frequency}kHz_{args.component}';conversion=json.loads((ROOT/'outputs/frequency_conversion.json').read_text());OMEGA=next(r['omega_nd'] for r in conversion['frequencies'] if r['f_kHz']==args.frequency)
    source=ROOT/'outputs'/f'{TAG}_v3_guide.json';guide=json.loads(source.read_text());assert guide['status'] in ['closed_verified','open_at_domain_boundaries'],guide['status']
    allrecords=guide['records']
    if 'reverse_half_start_index' in guide:
        j=guide['reverse_half_start_index'];records=allrecords[:j][::-1]+allrecords[j+1:]
    else:records=allrecords[:-1] if guide['status']=='closed_verified' else allrecords
    if guide['status']=='closed_verified':
        while len(records)>2 and all(abs(records[0]['point'][key]-records[-1]['point'][key])<1e-7 for key in ['xi','alpha_r','beta_r']):
            assert mac(np.load(records[0]['vector_file']),np.load(records[-1]['vector_file']))>1-1e-7
            records=records[:-1]
    if guide['status']=='closed_verified':
        for j in range(30,len(records)):
            if all(abs(records[0]['point'][key]-records[j]['point'][key])<1e-7 for key in ['xi','alpha_r','beta_r']):
                assert mac(np.load(records[0]['vector_file']),np.load(records[j]['vector_file']))>1-1e-7
                records=records[:j];break
    for order,rec in enumerate(records):rec['order']=float(order)
    expanded=[];extra=0
    for j,rec in enumerate(records):
        expanded.append(rec)
        if j==len(records)-1:break
        nxt=records[j+1];xa,xb=rec['point']['xi'],nxt['point']['xi']
        crossed=[i for i in range(int(np.floor(min(xa,xb)))+1,int(np.ceil(max(xa,xb)))) if min(xa,xb)+1e-9<i<max(xa,xb)-1e-9]
        if xb<xa:crossed=crossed[::-1]
        for i in crossed:
            theta=(i-xa)/(xb-xa);p=dict(rec['point'])
            for key in ['alpha_r','beta_r']:p[key]=(1-theta)*rec['point'][key]+theta*nxt['point'][key]
            p['xi']=float(i);extra+=1
            expanded.append(dict(index=100000+extra,order=j+theta,point=p,vector_file=(rec if theta<.5 else nxt)['vector_file'],crossing_sources=[rec['index'],nxt['index']]))
    wave=ROOT/'work'/f'{TAG}_final';wave.mkdir(exist_ok=True);native=[]
    for rec in expanded:
        r=rec['point'];i=round(r['xi'])
        if abs(i-r['xi'])>1e-9:continue
        if native and abs(native[-1]['xi']-i)<1e-9 and abs(native[-1]['beta_r']-r['beta_r'])<1e-8:continue
        p=Pencil(case,i,OMEGA);q=np.load(rec['vector_file']);beta=r['beta_r'];aa,q,rr=p.solve_alpha(beta,r['alpha_r'],q,True)
        # Refine beta with the original spatial polynomial rather than rounding Im(beta).
        for k in range(6):
            if abs(aa.imag)<1e-11:break
            beta-=aa.imag/rr['alpha_beta_i'];aa,q,rr=p.solve_alpha(beta,aa,q,True)
        assert abs(aa.imag)<1e-8
        row=dict(case=case,kind=0,xi=float(i),x_mm=CASES[case]['x_mm'][i-1],x_flat_mm=CASES[case]['x_mm'][i-1]-CASES[case]['flat_start_mm'],
            alpha_r=aa.real,alpha_i=aa.imag,beta_r=beta,beta_i=0.,omega_r=OMEGA,omega_i=0.,growth=-aa.imag,
            right_residual=rr['right_residual'],left_residual=rr['left_residual'],raw_MAC=rr['raw_MAC'],guide_index=rec['index'],
            alpha_beta_i=rr['alpha_beta_i'],guide_order=rec['order'])
        path=wave/f'q{len(native):04d}.npy';np.save(path,q);row['vector_file']=str(path);native.append(row)
    folds=[]
    for j in range(1,len(records)-1):
        l,m,r=records[j-1],records[j],records[j+1]
        if (m['point']['xi']-l['point']['xi'])*(r['point']['xi']-m['point']['xi']) < -1e-12:
            left=dict(l,turning_guess=m['point']);folds.append(tip(case,left,r,len(folds)))
    points=sorted(native+folds,key=lambda q:q['guide_order'])
    assert len(points)==len(native)+len(folds)
    assert len(folds)>0,'No x turning point found'
    # Every crossed original station needs both boundaries, with no missing stations.
    stations={int(q['xi']) for q in native};expected=set(range(min(stations),max(stations)+1));assert stations==expected
    counts={i:sum(round(q['xi'])==i for q in native) for i in stations};assert all(n>=2 and n%2==0 for n in counts.values()),set(counts.values())
    if guide['status']=='closed_verified':
        repeated=dict(points[0],kind=2,verified_closure=True);points.append(repeated)
    checks=[]
    for k in sorted(set(np.linspace(0,len(points)-1,12,dtype=int).tolist())):
        r=points[k]
        if r['kind']!=0:continue
        fresh,q,receipt=original(case,int(r['xi']),r,r['beta_r'],f'{TAG}_final_check{k:04d}')
        difference=abs(complex(fresh['alpha_r'],fresh['alpha_i'])-complex(r['alpha_r'],r['alpha_i']));overlap=mac(q,np.load(r['vector_file']))
        assert difference<1e-7 and overlap>1-1e-7 and abs(fresh['growth'])<1e-8
        checks.append(dict(point_index=k,alpha_difference=difference,raw_MAC=overlap,receipt=receipt))
    symmetry=[]
    for r in native:
        p=Pencil(case,int(r['xi']),OMEGA);pc=Pencil(case,int(r['xi']),-OMEGA);a=complex(r['alpha_r'],r['alpha_i']);q=np.load(r['vector_file'])
        error=max(float(norm((-1)**(k+j)*pc.C[k][j]-p.C[k][j].conj()))/max(1.,float(norm(p.C[k][j]))) for k in range(3) for j in range(3))
        res=float(np.linalg.norm(sum((-np.conj(a))**k*(-r['beta_r'])**j*pc.C[k][j] for k in range(3) for j in range(3))@q.conj()))
        assert error<1e-12 and res<1e-7;symmetry.append(dict(xi=r['xi'],matrix_symmetry_error=error,conjugate_residual=res))
    result=dict(case=case,f_kHz=args.frequency,component=args.component,omega_r=OMEGA,status=guide['status'],points=points,native_point_count=len(native),native_station_count=len(stations),physical_tips=folds,additional_native_crossings_corrected=extra,
        original_checks=checks,symmetry=symmetry,original_inputs=snapshot(),
        maximum_growth_error=max(abs(q['growth']) for q in points),maximum_right_residual=max(q['right_residual'] for q in points),
        x_range_mm=[min(q['x_mm'] for q in points),max(q['x_mm'] for q in points)],negative_beta_range=[min(q['beta_r'] for q in points),max(q['beta_r'] for q in points)])
    path=ROOT/'outputs'/f'{TAG}_final.json';path.write_text(json.dumps(result,indent=2)+'\n');print('FINALIZE',TAG,len(points),result['x_range_mm'],result['negative_beta_range'],flush=True)

if __name__=='__main__':main()
