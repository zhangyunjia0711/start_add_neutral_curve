"""Follow neutral roots from A(+50) through zero to conjugate B(+50)."""
import concurrent.futures,datetime,json,shutil,sys,time
from pathlib import Path
import numpy as np
FREQUENCY=int(sys.argv[1]) if len(sys.argv)>1 else 50
assert FREQUENCY in [50,100]
BASE=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°')
SOURCE=BASE/'neutral_xbeta_four_omega_20261026'
AUDIT=BASE/('neutral_frequency_connection_20261026'+('_f100' if FREQUENCY==100 else ''))
DEST=Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate/75°/outputs/neutral_frequency_connection_20261026')
if FREQUENCY==100:DEST=DEST/'f100'
sys.path.insert(0,str(SOURCE/'work'))
import common,core,solver
from common import CASES,SOLVER_HASH,sha,snapshot
from solver import Pencil,mac,original,PREVIOUS
core.ROOT=AUDIT;common.driver.ROOT=AUDIT;common.driver.CASE=SOURCE/'cases'
W=6.0727867021443425*FREQUENCY/50

def rows(case,root,filename):
    data=json.loads((root/'outputs'/filename).read_text());i=CASES[case]['seed_station']
    return sorted([p for p in data['points'] if p['kind']==0 and int(p['xi'])==i],key=lambda p:p['beta_r'])

def one(job):
    case,side=job;i=CASES[case]['seed_station'];tag=f'{case}_{side}'
    A=rows(case,SOURCE,f'{case}_f{FREQUENCY}kHz_A_final.json');B=rows(case,SOURCE,f'{case}_f{FREQUENCY}kHz_B_final.json')
    zero=rows(case,PREVIOUS,f'{case}_final.json');assert len(A)==len(B)==len(zero)==2
    k=0 if side=='lower' else 1;start=A[k];expected=B[1-k]
    w=W;b=start['beta_r'];a=complex(start['alpha_r'],start['alpha_i']);q=np.load(start['vector_file'])
    p=Pencil(case,i,w);a,q,r=p.solve_alpha(b,a,q,True)
    wave=AUDIT/'vectors'/tag;wave.mkdir(parents=True,exist_ok=True)
    records=[];attempts=[];checks=[];step=.1;started=time.monotonic()
    checkpoints=[W,W/2,0.,-W/2,-W];next_checkpoint=0
    def save(previous_mac):
        nonlocal next_checkpoint
        path=wave/f'q{len(records):04d}.npy';np.save(path,q)
        row=dict(index=len(records),case=case,side=side,station=i,x_mm=CASES[case]['x_mm'][i-1],omega_nd=float(w),f_kHz=float(FREQUENCY*w/W),beta=float(b),alpha_r=float(a.real),alpha_i=float(a.imag),growth=float(-a.imag),right_residual=r['right_residual'],left_residual=r['left_residual'],previous_vector_MAC=previous_mac,vector_file=str(path),vector_sha256=sha(path))
        records.append(row)
        if next_checkpoint<len(checkpoints) and abs(w-checkpoints[next_checkpoint])<1e-12:
            state=dict(omega_r=w,beta_r=b,alpha_r=a.real,alpha_i=a.imag)
            fresh,v,receipt=original(case,i,state,b,f'{tag}_f{row["f_kHz"]:+g}_neutral_connection')
            error=abs(complex(fresh['alpha_r'],fresh['alpha_i'])-a);overlap=mac(v,q)
            assert error<1e-7 and overlap>1-1e-7 and abs(fresh['growth'])<1e-8
            checks.append(dict(path_index=row['index'],f_kHz=row['f_kHz'],beta=b,alpha_difference=error,raw_MAC=overlap,original_growth=fresh['growth'],receipt=receipt))
            next_checkpoint+=1
        (AUDIT/'outputs'/f'{tag}_progress.json').write_text(json.dumps(dict(records=records,attempts=attempts,original_checks=checks),indent=2)+'\n')
    save(1.)
    while w>-W+1e-12:
        nw=max(-W,w-step)
        if next_checkpoint<len(checkpoints) and nw<checkpoints[next_checkpoint]:nw=checkpoints[next_checkpoint]
        dw=nw-w;ab=complex(r['alpha_beta_r'],r['alpha_beta_i']);aw=complex(r['alpha_omega_r'],r['alpha_omega_i'])
        assert abs(ab.imag)>1e-7,('Frequency projection fold requires pseudoarclength',case,side,w,b)
        dbdw=-aw.imag/ab.imag;nb=b+dbdw*dw;na=a+(aw+ab*dbdw)*dw;qq=q.copy();pn=Pencil(case,i,nw)
        try:
            for iteration in range(10):
                na,qq,rr=pn.solve_alpha(nb,na,qq,True)
                if abs(na.imag)<2e-11:break
                delta=-na.imag/rr['alpha_beta_i'];assert abs(delta)<1.,('Excessive beta correction',delta)
                nb+=delta
            overlap=mac(qq,q)
            assert abs(na.imag)<1e-9 and overlap>.995 and abs(nb-b)<1. and rr['right_residual']<1e-7
        except (AssertionError,RuntimeError) as e:
            step/=2;attempts.append(dict(omega=w,beta=b,step=step,reason=str(e)))
            assert step>1e-6,attempts[-1];continue
        w,b,a,q,r=nw,nb,na,qq,rr;save(overlap)
        if overlap>.999:step=min(.18,step*1.1)
        if len(records)%40==0:print('PATH',tag,len(records),w,b,flush=True)
        assert len(records)<3000 and time.monotonic()-started<300
    target_a=-complex(expected['alpha_r'],expected['alpha_i']).conjugate();target_b=-expected['beta_r'];target_q=np.load(expected['vector_file']).conj()
    endpoint=dict(beta_difference=float(abs(b-target_b)),alpha_difference=float(abs(a-target_a)),conjugate_vector_MAC=mac(q,target_q),expected_beta=target_b,expected_alpha=[target_a.real,target_a.imag],source_B_point=expected)
    assert endpoint['beta_difference']<1e-7 and endpoint['alpha_difference']<1e-7 and endpoint['conjugate_vector_MAC']>1-1e-7,endpoint
    atzero=next(v for v in records if abs(v['omega_nd'])<1e-12);z=zero[k]
    zero_check=dict(beta_difference=abs(atzero['beta']-z['beta_r']),alpha_difference=abs(complex(atzero['alpha_r'],atzero['alpha_i'])-complex(z['alpha_r'],z['alpha_i'])),raw_MAC=mac(np.load(atzero['vector_file']),np.load(z['vector_file'])))
    assert zero_check['beta_difference']<1e-7 and zero_check['alpha_difference']<1e-7 and zero_check['raw_MAC']>1-1e-7
    intermediate=[]
    if FREQUENCY==100:
        for sign,component in [(1,'A'),(-1,'B')]:
            ref=rows(case,SOURCE,f'{case}_f50kHz_{component}_final.json')[k if sign==1 else 1-k]
            node=next(r for r in records if abs(r['f_kHz']-50*sign)<1e-9)
            target_a=complex(ref['alpha_r'],ref['alpha_i']) if sign==1 else -complex(ref['alpha_r'],ref['alpha_i']).conjugate()
            refq=np.load(ref['vector_file']);refq=refq if sign==1 else refq.conj()
            check=dict(f_kHz=50*sign,beta_difference=abs(node['beta']-sign*ref['beta_r']),alpha_difference=abs(complex(node['alpha_r'],node['alpha_i'])-target_a),raw_MAC=mac(np.load(node['vector_file']),refq))
            assert check['beta_difference']<1e-7 and check['alpha_difference']<1e-7 and check['raw_MAC']>1-1e-7
            intermediate.append(check)
    result=dict(case=case,side=side,frequency_limit_kHz=FREQUENCY,status='continuous_neutral_path_verified',station=i,x_mm=CASES[case]['x_mm'][i-1],source_A=start,records=records,attempts=attempts,original_checks=checks,zero_frequency_check=zero_check,intermediate_slice_checks=intermediate,negative_endpoint_matches_conjugate_B=endpoint,maximum_growth_error=max(abs(t['growth']) for t in records),minimum_neighbor_MAC=min(t['previous_vector_MAC'] for t in records))
    (AUDIT/'outputs'/f'{tag}_path.json').write_text(json.dumps(result,indent=2)+'\n')
    print('CONNECTED',tag,len(records),endpoint['beta_difference'],endpoint['alpha_difference'],endpoint['conjugate_vector_MAC'],flush=True)
    return result

def main():
    assert not (AUDIT/'outputs/connection_report.json').exists(),'Preserve completed verification'
    for d in [AUDIT/'outputs',AUDIT/'runs',DEST]:d.mkdir(parents=True,exist_ok=True)
    snapshot();jobs=[(case,side) for case in ['bluntplate','rough'] for side in ['lower','upper']]
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:results=list(pool.map(one,jobs))
    report=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),method='Fixed native x; real frequency continuation and real beta neutral correction; original spatial PEP',source=str(SOURCE),audit_root=str(AUDIT),solver_sha256=SOLVER_HASH,paths=results,protected_inputs=snapshot())
    path=AUDIT/'outputs/connection_report.json';path.write_text(json.dumps(report,indent=2)+'\n');shutil.copy2(path,DEST/path.name)
    shutil.copy2(__file__,DEST/'verify.py')
    for r in results:
        for c in r['original_checks']:
            source=AUDIT/'runs'/c['receipt']['run'];target=DEST/'original_checks'/source.name;target.mkdir(parents=True,exist_ok=True)
            for relative in ['config.ini','solver.log','summary.json','out/neutral_field.plt','out/neutral_run.log']:
                p=source/relative
                if p.exists():
                    dst=target/relative;dst.parent.mkdir(exist_ok=True);shutil.copy2(p,dst)
    print('ALL_PATHS_CONNECTED',len(results),flush=True)
if __name__=='__main__':main()
