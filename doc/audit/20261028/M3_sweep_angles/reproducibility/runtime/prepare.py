"""Isolated angle calculations; protect all original input files.

Priority is the containing 1.1m directory, as explicitly requested by user.
Missing fine-grid full fields are kept pending, never replaced by coarse fields.
"""
import configparser,datetime,hashlib,json,shutil,sys
from pathlib import Path
import h5py,numpy as np
from scipy.interpolate import PchipInterpolator
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from inventory import read
WIN=HERE.parents[1]
BASE=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate')
ROOT=BASE/'neutral_m3_sweep_angles_20261027'
PREV=BASE/'75°/neutral_m3_frequency_limits_20261027'
REF=BASE/'75°/neutral_xbeta_four_omega_20261026/cases/bluntplate'
U=1810.6307411853716;L=.035;T=226.5;RHO=.018410817069110273
SOURCES=[(0,'0°sweptangle/rough'),(10,'10°swptangle/baseflow'),(30,'30°/1.1m'),(45,'45°swptangle/1.1m'),(60,'60°/1.1m'),(73,'73°/1.1m'),(74,'74°/rough'),(74.2,'74.2'),(74.5,'74.5°/rough')]
FREQ=[-240.2,-237.14,-200,-150,-112,-110,-100,-75,-50,-25,0,25,50,75,100,125,150,175,200,230,250,264.53,275,290,291.6]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    assert not (ROOT/'outputs/input_manifest.json').exists(),'Preserve completed preparation'
    for d in ['work','outputs','cases','runs','matrix_exports','native_pencils','inputs']:(ROOT/d).mkdir(parents=True,exist_ok=True)
    for name in ['common.py','run_adapter.py','core.py','solver.py','contour.py','contour3d.py','physical.py','tip_correction.py','finalize.py']:
        s=(PREV/'work'/name).read_text()
        if name=='solver.py':
            s=s.replace("PREVIOUS=ROOT.parent/'neutral_xbeta_omega0_20261026'","PREVIOUS=ROOT")
            s=s.replace("prior=json.loads((directory/'coefficients.json').read_text())","from core import Pencil as BetaPencil\n    BetaPencil(case,i)\n    prior=json.loads((directory/'coefficients.json').read_text())")
        if name=='physical.py':s=s.replace("p.baseC=[list(row) for row in p.C]","p.baseC=[list(row) for row in p.C]")
        if name=='contour.py':
            start=s.index('        path=ROOT/');end=s.index('\n    def pencil',start)
            s=s[:start]+"        return Pencil(self.case,i,self.omega)\n"+s[end:]
        if name=='contour3d.py':
            s=s.replace("and mac(q0,q)>1-1e-7:","and mac(q0,q)>1-1e-7 and np.dot(t,t0)>.7:")
        (ROOT/'work'/name).write_text(s)
    with h5py.File(REF/'in/grid_vertical2d_N.h5') as h:rg=h['grid'][...][:,:,:,0]
    rx=rg[0,:,0]*L;ry=rg[0,:,1]*L
    rt=np.column_stack([ry/.035,-rx/.035]);rt[rx>=-1e-12]=[1,0];rn=np.column_stack([-rt[:,1],rt[:,0]])
    refn=np.einsum('jik,ik->ji',rg-rg[0],rn)*L
    protected=[];cases={};selection=[]
    for angle,relative in SOURCES:
        src=WIN/relative;tag='s'+str(angle).replace('.','p');dest=ROOT/'cases'/tag
        (dest/'in').mkdir(parents=True,exist_ok=True)
        paths=[src/'plot3dx.unf']+[src/(k+'.plot3dq.unf') for k in 'RUVWTP']
        assert all(p.exists() for p in paths),src
        protected.extend(dict(path=str(p),sha256=sha(p),role='original_angle_input') for p in paths)
        g=read(paths[0],True);q=np.stack([read(p) for p in paths[1:]],axis=-1)
        ni,nj=g.shape[:2];wall=g[:,0];assert np.all(np.diff(wall[:,0])>0),relative
        normal=np.column_stack([wall[:,0]/.035,wall[:,1]/.035]);normal[wall[:,0]>=-1e-12]=[0,1]
        assert np.max(abs(np.linalg.norm(normal,axis=1)-1))<1e-7
        nn=np.einsum('ijk,ik->ij',g-wall[:,None,:],normal);assert np.all(np.diff(nn,axis=1)>0)
        target=np.empty((400,ni));newq=np.empty((400,ni,6))
        for i in range(ni):
            target[:,i]=np.array([np.interp(wall[i,0],rx,refn[j]) for j in range(400)])
            # Clip the outer domain by uniform scale only when source extent is shorter.
            target[:,i]*=min(1.,nn[i,-1]/target[-1,i]);target[0,i]=0
            target[:,i]=np.minimum(target[:,i],nn[i,-1])
            newq[:,i]=PchipInterpolator(nn[i],q[i],axis=0,extrapolate=False)(target[:,i])
        assert np.isfinite(newq).all() and np.all(newq[:,:,0]>0) and np.all(newq[:,:,4]>0),(relative,int(np.count_nonzero(~np.isfinite(newq))),np.nanmin(newq[:,:,0]),np.nanmin(newq[:,:,4]))
        ng=wall[None,:,:]+target[:,:,None]*normal[None,:,:]
        gd=np.zeros((400,ni,2,2));gd[:,:,:,0]=ng/L
        bd=np.zeros((400,ni,6,2));bd[:,:,:,0]=newq/np.array([RHO,U,U,U,T,RHO*U*U])
        files=[]
        for name,key,data in [('grid_vertical2d_N.h5','grid',gd),('flow_vertical2d_N.h5','baseflow',bd)]:
            p=dest/'in'/name
            with h5py.File(p,'w') as h:h.create_dataset(key,data=data,compression='gzip',compression_opts=1)
            with h5py.File(p) as h:assert np.array_equal(h[key][...],data)
            files.append(dict(path=str(p),sha256=sha(p)))
        cfg=configparser.ConfigParser();cfg.read(REF/'config.ini');seed=int(np.argmin(abs(wall[:,0]-rx[118])))+1
        cfg['filenames']['dir']=str(dest/'in')+'/';cfg['domain'].update({'in':str(ni),'jn':'400','iloc':str(seed),'istart':'1','iend':str(ni)})
        with (dest/'config.ini').open('w') as f:cfg.write(f)
        cfgfile=dest/'config.ini';protected.extend(files+[dict(path=str(cfgfile),sha256=sha(cfgfile),role='isolated_angle_config')])
        cases[tag]=dict(name=tag,angle=angle,config=str(cfgfile),in_dir=str(dest/'in'),source_dir=str(src),stations=ni,normal_points=400,x_mm=(wall[:,0]*1000).tolist(),wall_y_mm=(wall[:,1]*1000).tolist(),flat_start_mm=0.,seed_station=seed)
        row=dict(angle=angle,case=tag,source=str(src),original_shape=[ni,nj],calculation_shape=[ni,400],source_priority='1.1m dense' if '/1.1m' in str(src) else 'No 1.1m directory: available coarse input',normal_resampling='PCHIP of actual native primitive profiles versus projected wall distance, 75-degree normal distribution mapped by wall x; not CFD refinement',files=files)
        selection.append(row);print('PREPARED',angle,relative,ni,nj,seed,flush=True)
    # Exact 75-degree anchor, linked to immutable previous inputs and coefficient caches.
    (ROOT/'cases/s75').symlink_to(REF,target_is_directory=True)
    priorcases=json.loads((PREV/'outputs/cases.json').read_text());c=dict(priorcases['bluntplate']);c.update(name='s75',angle=75);cases['s75']=c
    selection.append(dict(angle=75,case='s75',source=str(REF),source_priority='Existing verified 75-degree 1.1m baseline',reuse=True))
    for p in [REF/'config.ini',REF/'in/grid_vertical2d_N.h5',REF/'in/flow_vertical2d_N.h5']:protected.append(dict(path=str(p),sha256=sha(p),role='verified75_anchor'))
    (ROOT/'outputs/cases.json').write_text(json.dumps(cases,indent=2)+'\n')
    (ROOT/'outputs/input_manifest.json').write_text(json.dumps(protected,indent=2)+'\n')
    (ROOT/'outputs/input_selection.json').write_text(json.dumps(dict(selected=selection,pending_dense_fullfields=[50,55,63,65,68,70,80],empty_or_incomplete=[15,78],rule='User: all results within 1.1m are dense-grid; never override by file name.'),indent=2)+'\n')
    (ROOT/'outputs/frequency_conversion.json').write_text(json.dumps(dict(Lref_m=L,Uref_m_s=U,definition='omega_nd=2*pi*f_kHz*1000*Lref/Uref; signed cyclic physical frequency',frequencies=[dict(f_kHz=f,omega_nd=float(2*np.pi*1000*f*L/U)) for f in FREQ]),indent=2)+'\n')
    for p in HERE.glob('*.py'):shutil.copy2(p,ROOT/'work'/p.name)
    print('ROOT',ROOT,'protected',len(protected),flush=True)
if __name__=='__main__':main()
