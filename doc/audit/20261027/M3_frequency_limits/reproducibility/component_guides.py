"""Whole neutral loops for the newly found connected downstream pocket."""
import concurrent.futures,json,subprocess,sys
from pathlib import Path
import numpy as np
ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/neutral_m3_frequency_limits_20261027')
sys.path.insert(0,str(ROOT/'work'))
from common import CASES,driver,sha
from contour import Family
from contour3d import bordered,SCALE
from solver import mac
W=6.0727867021443425/50
FS=[-112.,-150.,-200.]
def prepare():
    path=ROOT/'work/contour.py';s=path.read_text();anchor="    fresh,v,receipt=original(case,i,s,s['beta_r'],f'{case}_f{frequency}kHz_S_neutral_seed')\n";assert s.count(anchor)==1
    if "s.get('kind')=='operator_guide'" not in s:s=s.replace(anchor,"    if s.get('kind')=='operator_guide':\n        return dict(station=i,point=s,receipt=dict(model='Original-operator guide seed; independent physical checks performed after finalization'),component=component,frequency_kHz=frequency,source_record=s),q\n"+anchor)
    path.write_text(s)
    path=ROOT/'work/contour3d.py';s=path.read_text();anchor="    family=Family(case,omega);seed,q=seed_root(case,f,component);s=seed['point'];pencil=Pencil(case,seed['station'],omega)\n    _,q,_,p=qep(pencil.spatial(s['beta_r']),complex(s['alpha_r'],s['alpha_i']),q,True)\n"
    replacement="    family=Family(case,omega);seed,q=seed_root(case,f,component);s=seed['point']\n    C,dC,_,_=family.pencil(float(seed['station']));K=[sum(s['beta_r']**j*C[k][j] for j in range(3)) for k in range(3)]\n    _,q,_,p=qep(K,complex(s['alpha_r'],s['alpha_i']),q,True)\n"
    if anchor in s:s=s.replace(anchor,replacement)
    else:assert replacement in s
    path.write_text(s)
    conversion=ROOT/'outputs/frequency_conversion.json';d=json.loads(conversion.read_text());d['frequencies']=[dict(f_kHz=f,omega_nd=W*f) for f in FS];conversion.write_text(json.dumps(d,indent=2)+'\n');seeds=[]
    for case in CASES:
        arc=json.loads((ROOT/'outputs'/f'{case}_negative_secondary_arc.json').read_text())
        for f in FS:
            ref=min(arc['records'],key=lambda r:abs(r['point']['f_kHz']-f));r=ref['point'];q=np.load(ref['vector_file']);p=np.load(ref['left_vector_file']);z=np.array([r['xi'],r['alpha_r'],r['beta_r']])/SCALE;family=Family(case,W*f)
            for it in range(30):
                C,dC,i,t=family.pencil(float(z[0]*SCALE[0]));e,J,_,v,w,res,_=bordered(C,dC,z,q,p)
                if abs(e)<1e-12 and res<1e-8:break
                dz=np.linalg.solve(J[:,:2],-np.array([e.real,e.imag]));size=np.linalg.norm(dz)
                if size>.05:dz*=.05/size
                z[:2]+=dz
            assert abs(e)<1e-12 and res<1e-8 and mac(q,v)>.95
            file=ROOT/'outputs'/f'{case}_secondary_f{f:g}_guide_seed.npy';np.save(file,v);xi,alpha,beta=z*SCALE
            row=dict(case=case,f_kHz=f,station=float(xi),xi=float(xi),side='lower',kind='operator_guide',alpha_r=float(alpha),alpha_i=0.,beta_r=float(beta),beta_i=0.,omega_r=float(f*W),omega_i=0.,growth=0.,right_residual=res,vector_file=str(file),vector_sha256=sha(file),source_arc=str(ROOT/'outputs'/f'{case}_negative_secondary_arc.json'));seeds.append(row);print('COMPONENT_SEED',case,f,xi,alpha,beta,flush=True)
    (ROOT/'outputs/neutral_seeds.json').write_text(json.dumps(dict(records=seeds),indent=2)+'\n')
def one(job):
    case,f=job;tag=f'{case}_f{float(f)}kHz_S';result=dict(case=case,f_kHz=f,stages=[])
    for script,suffix in [('contour3d.py','v3_guide'),('finalize.py','final')]:
        path=ROOT/'outputs'/f'{tag}_{suffix}.json';log=ROOT/'outputs'/f'{tag}_{suffix}.log'
        if not path.exists():
            with log.open('w') as stream:p=subprocess.run([sys.executable,str(ROOT/'work'/script),case,str(f),'S'],env=driver.ENV,stdout=stream,stderr=subprocess.STDOUT)
            code=p.returncode
        else:code=0
        d=json.loads(path.read_text()) if path.exists() else {};result['stages'].append(dict(script=script,exit=code,status=d.get('status','NO_RESULT'),records=len(d.get('points',d.get('records',[]))),log=str(log)));print('COMPONENT_STAGE',case,f,result['stages'][-1],flush=True)
        if code or d.get('status')!='closed_verified':break
    return result
def main():
    if not (ROOT/'outputs/neutral_seeds.json').exists() or len(json.loads((ROOT/'outputs/neutral_seeds.json').read_text())['records'])<6:prepare()
    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as pool:results=list(pool.map(one,[(case,f) for case in CASES for f in FS]))
    (ROOT/'outputs/component_curves_summary.json').write_text(json.dumps(results,indent=2)+'\n')
if __name__=='__main__':main()
