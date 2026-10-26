"""Correct native neutral points and actual profile-interpolated x turning points."""
import argparse,json,csv
import numpy as np
from scipy.sparse.linalg import norm
from common import ROOT,CASES,snapshot
from solver import Pencil,qep,mac,original
from physical import Pair

def physical_evaluate(pair,theta,alpha,beta,q):
    p=pair.pencil(float(theta));K=p.transverse(alpha);b,v,r,left=qep(K,beta,q,True)
    Pa=sum(b**j*(p.C[1][j]+2*alpha*p.C[2][j]) for j in range(3));Pb=K[1]+2*b*K[2]
    ba=-np.vdot(left,Pa@v)/np.vdot(left,Pb@v)
    return np.array([b.imag,ba.imag]),b,v,r

def tip(case,left,right,index):
    lp,rp=left['point'],right['point'];f=abs(lp['beta_alpha_i'])/(abs(lp['beta_alpha_i'])+abs(rp['beta_alpha_i']))
    xi=(1-f)*lp['xi']+f*rp['xi'];alpha=(1-f)*lp['alpha_r']+f*rp['alpha_r'];i=int(np.floor(xi));theta=xi-i
    pair=Pair(case,i);beta=(1-f)*lp['beta_r']+f*rp['beta_r'];q=np.load(left['vector_file']);trace=[]
    for it in range(12):
        F,b,v,r=physical_evaluate(pair,theta,alpha,beta,q);beta=b.real;q=v
        trace.append(dict(theta=theta,alpha_r=alpha,beta_r=beta,F=F.tolist(),residual=r))
        if max(abs(F))<1e-10:break
        ht=min(1e-4,max(1e-6,min(theta,1-theta)/4));ha=1e-3
        direction=1 if theta<.999 else -1
        Ft,_,_,_=physical_evaluate(pair,theta+direction*ht,alpha,beta,q)
        Fa,_,_,_=physical_evaluate(pair,theta,alpha+ha,beta,q)
        J=np.column_stack([(Ft-F)/(direction*ht),(Fa-F)/ha]);delta=np.linalg.solve(J,-F)
        damping=1.
        while not 0<theta+damping*delta[0]<1:damping/=2;assert damping>1e-6,('Physical tip outside native bracket',theta,delta)
        theta+=damping*delta[0];alpha+=damping*delta[1]
    assert max(abs(F))<1e-9
    p=pair.pencil(theta);aa,q,rr=p.solve_alpha(beta,alpha,q,True);assert abs(aa.imag)<1e-8
    generation=json.loads((pair.build(theta)/'generation.json').read_text());endpoint_checks=[]
    for t,j in [(0.,i),(1.,i+1)]:
        e=pair.pencil(t);native=Pencil(case,j)
        errors=[float(norm(a-b))/max(1.,float(norm(a))) for a,b in zip(native.spatial(beta),e.spatial(beta))]
        assert max(errors)<1e-9,('Profile endpoint operator mismatch',errors)
        endpoint_checks.append(dict(station=j,relative_matrix_errors=errors))
    point,q,receipt=original(pair.build(theta),5,dict(omega_r=0.,beta_r=beta,alpha_r=alpha,alpha_i=0.),beta,f'{case}_physical_tip{index}')
    assert abs(point['growth'])<1e-8
    out=dict(case=case,kind=1,xi=i+theta,x_mm=generation['x_mm'],x_flat_mm=generation['x_mm']-CASES[case]['flat_start_mm'],
        alpha_r=point['alpha_r'],alpha_i=point['alpha_i'],beta_r=beta,beta_i=0.,omega_r=0.,omega_i=0.,growth=point['growth'],
        right_residual=rr['right_residual'],left_residual=rr['left_residual'],raw_MAC=rr['raw_MAC'],guide_left_index=left['index'],guide_right_index=right['index'],
        guide_order=(1-f)*left['order']+f*right['order'],
        native_bracket=[i,i+1],theta=theta,trace=trace,endpoint_checks=endpoint_checks,original_receipt=receipt,generation=str(pair.build(theta)/'generation.json'))
    path=ROOT/'outputs'/f'{case}_tip{index}.npy';np.save(path,q);out['vector_file']=str(path);return out

def main():
    parser=argparse.ArgumentParser();parser.add_argument('case');args=parser.parse_args();case=args.case
    source=ROOT/'outputs'/f'{case}_guide.json';guide=json.loads(source.read_text());assert guide['status'] in ['closed_verified','open_at_domain_boundaries'],guide['status']
    allrecords=guide['records']
    if 'reverse_half_start_index' in guide:
        j=guide['reverse_half_start_index'];records=allrecords[:j][::-1]+allrecords[j+1:]
    else:records=allrecords
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
    wave=ROOT/'work'/f'{case}_final';wave.mkdir(exist_ok=True);native=[]
    for rec in expanded:
        r=rec['point'];i=round(r['xi'])
        if abs(i-r['xi'])>1e-9:continue
        if native and abs(native[-1]['xi']-i)<1e-9 and abs(native[-1]['beta_r']-r['beta_r'])<1e-8:continue
        p=Pencil(case,i);q=np.load(rec['vector_file']);beta=r['beta_r'];aa,q,rr=p.solve_alpha(beta,r['alpha_r'],q,True)
        # Refine beta with the original spatial polynomial rather than rounding Im(beta).
        for k in range(6):
            if abs(aa.imag)<1e-11:break
            beta-=aa.imag/rr['alpha_beta_i'];aa,q,rr=p.solve_alpha(beta,aa,q,True)
        assert abs(aa.imag)<1e-8
        row=dict(case=case,kind=0,xi=float(i),x_mm=CASES[case]['x_mm'][i-1],x_flat_mm=CASES[case]['x_mm'][i-1]-CASES[case]['flat_start_mm'],
            alpha_r=aa.real,alpha_i=aa.imag,beta_r=beta,beta_i=0.,omega_r=0.,omega_i=0.,growth=-aa.imag,
            right_residual=rr['right_residual'],left_residual=rr['left_residual'],raw_MAC=rr['raw_MAC'],guide_index=rec['index'],
            alpha_beta_i=rr['alpha_beta_i'],guide_order=rec['order'])
        path=wave/f'q{len(native):04d}.npy';np.save(path,q);row['vector_file']=str(path);native.append(row)
    folds=[]
    for l,r in zip(records[:-1],records[1:]):
        if l['point']['beta_alpha_i']*r['point']['beta_alpha_i']<0:
            folds.append(tip(case,l,r,len(folds)))
    points=sorted(native+folds,key=lambda q:q['guide_order'])
    assert len(points)==len(native)+len(folds)
    assert len(folds)>0,'No x turning point found'
    # Every crossed original station needs both boundaries, with no missing stations.
    stations={int(q['xi']) for q in native};expected=set(range(min(stations),CASES[case]['stations']+1));assert stations==expected
    counts={i:sum(round(q['xi'])==i for q in native) for i in stations};assert set(counts.values())=={2},set(counts.values())
    checks=[]
    for k in sorted(set(np.linspace(0,len(points)-1,12,dtype=int).tolist())):
        r=points[k]
        if r['kind']!=0:continue
        fresh,q,receipt=original(case,int(r['xi']),r,r['beta_r'],f'{case}_final_check{k:04d}')
        difference=abs(complex(fresh['alpha_r'],fresh['alpha_i'])-complex(r['alpha_r'],r['alpha_i']));overlap=mac(q,np.load(r['vector_file']))
        assert difference<1e-7 and overlap>1-1e-7 and abs(fresh['growth'])<1e-8
        checks.append(dict(point_index=k,alpha_difference=difference,raw_MAC=overlap,receipt=receipt))
    symmetry=[]
    for r in native:
        p=Pencil(case,int(r['xi']));a=complex(r['alpha_r'],r['alpha_i']);q=np.load(r['vector_file'])
        error=max(float(norm((-1)**(k+j)*p.C[k][j]-p.C[k][j].conj()))/max(1.,float(norm(p.C[k][j]))) for k in range(3) for j in range(3))
        res=float(np.linalg.norm(sum((-np.conj(a))**k*(-r['beta_r'])**j*p.C[k][j] for k in range(3) for j in range(3))@q.conj()))
        assert error<1e-12 and res<1e-7;symmetry.append(dict(xi=r['xi'],matrix_symmetry_error=error,conjugate_residual=res))
    result=dict(case=case,omega_r=0.,status=guide['status'],points=points,native_point_count=len(native),native_station_count=len(stations),physical_tips=folds,additional_native_crossings_corrected=extra,
        original_checks=checks,symmetry=symmetry,original_inputs=snapshot(),
        maximum_growth_error=max(abs(q['growth']) for q in points),maximum_right_residual=max(q['right_residual'] for q in points),
        x_range_mm=[min(q['x_mm'] for q in points),max(q['x_mm'] for q in points)],negative_beta_range=[min(q['beta_r'] for q in points),max(q['beta_r'] for q in points)])
    path=ROOT/'outputs'/f'{case}_final.json';path.write_text(json.dumps(result,indent=2)+'\n');print('FINALIZE',case,len(points),result['x_range_mm'],result['negative_beta_range'],flush=True)

if __name__=='__main__':main()
