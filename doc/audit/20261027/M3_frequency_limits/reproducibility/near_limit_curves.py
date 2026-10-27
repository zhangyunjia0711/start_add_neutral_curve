"""Compute closed neutral loops using actual critical branches and profiles.

Every delivered root solves the frozen original spatial polynomial. Native
stations are all retained; small loops receive extra actual-profile samples.
"""
import concurrent.futures,json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/neutral_m3_frequency_limits_20261027')
sys.path.insert(0,str(ROOT/'work'))
from common import CASES,sha,snapshot
from solver import Pencil,mac,original
from physical import Pair
from trace_folds import solve
from fold_core import S,W
JOBS=[('bluntplate',1,f) for f in [230.,250.,275.,290.,291.6]]+[('rough',1,f) for f in [230.,250.,260.,264.5]]+[('bluntplate',-1,f) for f in [-113.,-114.,-114.12]]+[('rough',-1,f) for f in [-114.,-115.,-115.95]]
def setfreq(p,f):
    p.omega=f*W;p.C=[list(row) for row in p.baseC];p.C[0][0]=p.C[0][0]+p.omega*p.G
def root_pair(p,cap,qcap,f):
    a,b,fc=cap['alpha_r'],cap['beta_r'],cap['f_kHz'];q=qcap.copy();a=complex(a)
    for frequency in np.linspace(fc,f,max(2,int(abs(fc-f)/3)+2))[1:]:
        setfreq(p,float(frequency));a,q,r=p.solve_alpha(b,a,q,True)
    setfreq(p,f);trace=[]
    # Follow the local maximum growth, then bracket its two neutral boundaries.
    for it in range(60):
        a,q,r=p.solve_alpha(b,a,q,True);h=.02;am,_,_=p.solve_alpha(b-h,a,q);ap,_,_=p.solve_alpha(b+h,a,q);curvature=(ap.imag-2*a.imag+am.imag)/(h*h);der=r['alpha_beta_i'];trace.append(dict(beta=b,alpha_r=a.real,alpha_i=a.imag,beta_derivative=der,beta_curvature=curvature))
        if abs(der)<1e-9:break
        step=float(np.clip(-der/curvature,-1.,1.)) if curvature>0 else float(-np.sign(der)*.4);accepted=False
        for retry in range(12):
            aa,v,rr=p.solve_alpha(b+step,a,q,True)
            if aa.imag<a.imag+1e-12:b+=step;a,q=aa,v;accepted=True;break
            step/=2
        assert accepted,('Growth maximum iteration stalled',f,trace[-1])
    assert a.imag<0 and abs(r['alpha_beta_i'])<1e-6,('No interior instability',f,trace[-1])
    centre=dict(beta=b,alpha_r=a.real,alpha_i=a.imag);out=[]
    for direction in [-1,1]:
        lastb=b;lasta=a;lastq=q.copy();h=.05
        for k in range(100):
            bb=lastb+direction*h;aa,v,rr=p.solve_alpha(bb,lasta,lastq,True)
            if aa.imag>0:break
            lastb,lasta,lastq=bb,aa,v;h=min(1.,h*1.35)
        else:raise RuntimeError(('Neutral boundary not bracketed',f,direction,b))
        lo,hi=sorted([lastb,bb]);anchor_b=(lastb+bb)/2;anchor_a,anchor_q,_=p.solve_alpha(anchor_b,(lasta+aa)/2,lastq,True)
        def fun(beta):return p.solve_alpha(float(beta),anchor_a,anchor_q)[0].imag
        beta=brentq(fun,lo,hi,xtol=2e-11);alpha,vector,record=p.solve_alpha(beta,anchor_a,anchor_q,True);assert abs(alpha.imag)<1e-8
        out.append((dict(alpha_r=alpha.real,alpha_i=alpha.imag,beta_r=float(beta),growth=-alpha.imag,right_residual=record['right_residual'],raw_MAC=record['raw_MAC'],alpha_beta_i=record['alpha_beta_i']),vector))
    assert out[0][0]['beta_r']<out[1][0]['beta_r']
    return out,dict(centre=centre,peak_trace=trace)

def one(job):
    case,sign,f=job[:3];variant=job[3] if len(job)>3 else '';label=f'{case}_{"positive" if sign>0 else "negative"}'+('_'+variant if variant else '');tag=f'{case}'+('_'+variant if variant else '')+f'_f{f:g}kHz';source=ROOT/'outputs'/f'{label}_native_fold_curve.json';d=json.loads(source.read_text());records={int(p['station']):p for p in d['records']};inside=sorted(i for i,p in records.items() if sign*(p['f_kHz']-f)>0);assert inside and inside==list(range(min(inside),max(inside)+1))
    pairs={};cap_cache={};wave=ROOT/'work'/f'{tag}_neutral_vectors';wave.mkdir(exist_ok=True);audit=[]
    def cap_at(x):
        key=float(x).hex()
        if key in cap_cache:return cap_cache[key]
        i=min(int(np.floor(x)),CASES[case]['stations']-1);t=x-i
        if i not in pairs:pairs[i]=Pair(case,i,0.,tag+'_actual_profile')
        p=pairs[i].pencil(float(t));a=records.get(i,records[min(records)]);b=records.get(i+1,a);z=np.array([(1-t)*a[k]+t*b[k] for k in ['alpha_r','beta_r','f_kHz']])/S;qref=np.load((a if t<.5 else b)['vector_file']);cap,q,z=solve(p,z,qref,max_nfev=100);cap.update(xi=float(x),x_mm=(1-t)*CASES[case]['x_mm'][i-1]+t*CASES[case]['x_mm'][i],generation=str(pairs[i].build(float(t))/'generation.json'));cap_cache[key]=(cap,q,p);return cap_cache[key]
    tips=[]
    for lo,hi in [(min(inside)-1,min(inside)),(max(inside),max(inside)+1)]:
        if lo not in records and variant=='secondary':
            walk=json.loads((ROOT/'outputs'/f'{case}_negative_secondary_left_physical.json').read_text())['records'];outside=[r['xi'] for r in walk if r['xi']<hi and sign*(r['f_kHz']-f)<0];assert outside;lo=max(outside)
        else:assert lo in records and hi in records
        x=brentq(lambda x:cap_at(float(x))[0]['f_kHz']-f,float(lo),float(hi),xtol=1e-8);cap,q,p=cap_at(x);assert abs(cap['f_kHz']-f)<1e-6
        cap=dict(cap,kind=1,case=case,omega_r=f*W,f_kHz=f,alpha_i=0.,growth=0.);tips.append((cap,q))
    xleft,xright=tips[0][0]['xi'],tips[1][0]['xi'];samples=[float(i) for i in range(int(np.ceil(xleft)),int(np.floor(xright))+1)]
    extra=[]
    if len(samples)<14:
        extra=[float(x) for x in np.linspace(xleft+.01*(xright-xleft),xright-.01*(xright-xleft),25) if all(abs(x-i)>1e-7 for i in samples)]
    samples=sorted(samples+extra);low=[];high=[]
    def store(row,q,name):
        file=wave/(name+'.npy');np.save(file,q);row.update(vector_file=str(file),vector_sha256=sha(file),x_flat_mm=row['x_mm']-CASES[case]['flat_start_mm'],beta_i=0.,omega_i=0.);return row
    for k,x in enumerate(samples):
        if abs(x-round(x))<1e-9:
            i=round(x);cap=records[i];qcap=np.load(cap['vector_file']);p=Pencil(case,i,f*W);kind=0;generation=None;xmm=CASES[case]['x_mm'][i-1]
        else:cap,qcap,p=cap_at(x);kind=3;generation=cap['generation'];xmm=cap['x_mm']
        roots,checks=root_pair(p,cap,qcap,f);audit.append(dict(xi=x,kind=kind,checks=checks))
        for side,collection,(row,q) in zip(['lower','upper'],[low,high],roots):
            row.update(case=case,f_kHz=f,omega_r=f*W,kind=kind,xi=x,x_mm=xmm,generation=generation);collection.append(store(row,q,f'{k:04d}_{side}'))
    left=store(tips[0][0],tips[0][1],'left_tip');right=store(tips[1][0],tips[1][1],'right_tip');points=[left]+low+[right]+high[::-1]+[dict(left,kind=2,verified_closure=True)];checks=[]
    for k in sorted(set([0,len(points)-1,len(low)+1]+np.linspace(1,len(points)-2,6,dtype=int).tolist())):
        row=points[k]
        if row['kind']==2:continue
        domain=case if row['kind']==0 else Path(row['generation']).parent;i=int(row['xi']) if row['kind']==0 else 5
        fresh,v,receipt=original(domain,i,row,row['beta_r'],f'{tag}_actual_check{k}');difference=abs(complex(fresh['alpha_r'],fresh['alpha_i'])-complex(row['alpha_r'],row['alpha_i']));overlap=mac(v,np.load(row['vector_file']));assert abs(fresh['growth'])<1e-8 and difference<1e-7 and overlap>1-1e-7;checks.append(dict(point_index=k,alpha_difference=difference,MAC=overlap,receipt=receipt))
    native={round(p['xi']) for p in points if p['kind']==0};expected=set(range(int(np.ceil(xleft)),int(np.floor(xright))+1));assert native==expected
    result=dict(case=case,sign=sign,f_kHz=f,omega_r=f*W,status='closed_verified',points=points,native_stations=len(native),physical_extra_samples=len(extra),root_checks=audit,original_checks=checks,critical_source=str(source),critical_source_sha256=sha(source),x_range_mm=[left['x_mm'],right['x_mm']],beta_range=[min(p['beta_r'] for p in points),max(p['beta_r'] for p in points)],maximum_growth_error=max(abs(p['growth']) for p in points),maximum_right_residual=max(p['right_residual'] for p in points));path=ROOT/'outputs'/f'{tag}_neutral.json';path.write_text(json.dumps(result,indent=2)+'\n');print('LIMIT_CURVE',case,f,len(points),result['x_range_mm'],flush=True);return {k:result[k] for k in ['case','sign','f_kHz','status','x_range_mm','beta_range','maximum_growth_error','maximum_right_residual']}
def safe_one(job):
    case,sign,f=job[:3];variant=job[3] if len(job)>3 else '';tag=f'{case}'+('_'+variant if variant else '')+f'_f{f:g}kHz';file=ROOT/'outputs'/f'{tag}_neutral.json'
    if file.exists():
        d=json.loads(file.read_text());assert d['status']=='closed_verified';return {k:d[k] for k in ['case','sign','f_kHz','status','x_range_mm','beta_range','maximum_growth_error','maximum_right_residual']}
    try:return one(job)
    except Exception as error:
        import traceback
        d=dict(case=case,sign=sign,f_kHz=f,status='failed',reason=repr(error),traceback=traceback.format_exc());(ROOT/'outputs'/f'{tag}_attempt_failure.json').write_text(json.dumps(d,indent=2)+'\n');print('CURVE_FAILED',case,f,repr(error),flush=True);return d
def main():
    snapshot()
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:results=list(pool.map(safe_one,JOBS))
    (ROOT/'outputs/near_limit_curves_summary.json').write_text(json.dumps(results,indent=2)+'\n');snapshot()
if __name__=='__main__':main()
