"""Independent nontrivial section return and turning-point corrections."""
import json,datetime
import numpy as np
from common import ROOT,PREVIOUS,sha
from fixed_beta_pencil import OperatorLinearFamily
from profile_interpolation import mac
SCALE=np.array([100.,20.])

def main():
    source=ROOT/'outputs/M3b7_closed_arc_v2.json';d=json.loads(source.read_text());assert d['status']=='closed_contour_verified'
    family=OperatorLinearFamily();records=d['records'];grid=np.array(json.loads((PREVIOUS/'outputs/seed_catalog.json').read_text())['grids']['main_841']['x_nondim'])*35
    def ev(z,q,w):
        p=family.pencil(float(z[0]*100));r,v=p.evaluate(float(z[1]*20),q,w)
        i=p.native_bracket[0];r.update(xi=float(z[0]*100),x_mm=float((1-p.theta)*grid[i-1]+p.theta*grid[i]),
            native_bracket=p.native_bracket,theta=p.theta,alpha_i=0.)
        return r,v
    start=records[0];t=np.array(start['tangent']);zstart=np.array(start['z'])
    a,b=records[-3],records[-2];za,zb=np.array(a['z']),np.array(b['z']);sa,sb=np.dot(za-zstart,t),np.dot(zb-zstart,t)
    assert sa<0<=sb
    z=za+(-sa)/(sb-sa)*(zb-za);initial_z=z.copy();q=np.load(b['vector_file']);w=b['point']['omega_r'];trace=[]
    for k in range(10):
        r,q=ev(z,q,w);w=r['omega_r'];F=np.array([r['omega_i'],np.dot(z-zstart,t)])
        trace.append(dict(z=z.tolist(),point=r,F=F.tolist()))
        if abs(F[0])<1e-11 and abs(F[1])<1e-10:break
        J=np.array([[r['omega_xi_i']*100,r['omega_alpha_i']*20],t]);z-=np.linalg.solve(J,F)
    assert abs(r['omega_i'])<1e-10 and np.linalg.norm(z-zstart)<1e-7
    overlap=mac(q,np.load(start['vector_file']));assert overlap>1-1e-7
    closure=dict(initial_section_interpolation_z=initial_z.tolist(),initial_distance_to_start=float(np.linalg.norm(initial_z-zstart)),
        actual_corrected_point=r,scaled_return_error=float(np.linalg.norm(z-zstart)),raw_MAC_to_start=overlap,
        omega_difference=abs(r['omega_r']-start['point']['omega_r']),trace=trace,passed=True)
    folds=[]
    for a,b in zip(records,records[1:]):
        if a['tangent'][0]*b['tangent'][0]>=0:continue
        z=(np.array(a['z'])+np.array(b['z']))/2;q=np.load(a['vector_file']);w=a['point']['omega_r'];history=[]
        for it in range(18):
            r,q=ev(z,q,w);w=r['omega_r'];F=np.array([r['omega_i'],r['omega_alpha_i']]);row=dict(z=z.tolist(),point=r,F=F.tolist());history.append(row)
            if max(abs(F))<1e-10:break
            J=np.empty((2,2));J[0]=[r['omega_xi_i']*100,r['omega_alpha_i']*20]
            for j,h in enumerate([1e-5,1e-5]):
                zz=z.copy();zz[j]+=h;rp,_=ev(zz,q,w);zz[j]-=2*h;rm,_=ev(zz,q,w)
                J[1,j]=(rp['omega_alpha_i']-rm['omega_alpha_i'])/(2*h)
            delta=np.linalg.solve(J,-F);scale=min(1,.03/max(abs(delta[0]),1e-100),.1/max(abs(delta[1]),1e-100))
            row['J']=J.tolist();z+=scale*delta
        assert max(abs(F))<1e-8,('Unresolved actual turning point',r,F)
        path=ROOT/'work'/('M3b7_fold_q%d.npy'%len(folds));np.save(path,q)
        folds.append(dict(arc_segment=[a['index'],b['index']],z=z.tolist(),point=r,vector_file=str(path),
            vector_sha256=sha(path),trace=history,passed=True))
    report=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),
        arc_path=str(source),arc_sha256=sha(source),closure=closure,folds=folds,
        qualification='Closure of one continued M3 neutral component under stated adjacent-native-operator interpolation.')
    (ROOT/'outputs/closed_contour_validation.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print('CLOSURE',closure['initial_distance_to_start'],closure['scaled_return_error'],overlap)
    print('FOLDS',[(r['point']['xi'],r['point']['x_mm'],r['point']['omega_r']) for r in folds])

if __name__=='__main__':main()
