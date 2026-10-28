"""Use exact existing 1.1m/fpp_output profiles, with explicit unit conversion."""
import configparser,hashlib,json,shutil,sys,datetime
from pathlib import Path
import h5py,numpy as np
HERE=Path(__file__).resolve().parent
WIN=HERE.parents[1]
sys.path.insert(0,str(WIN/'work'))
from convert_plot3d_2d_to_h5 import read_grid,read_flow
BASE=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate')
TRIAL=BASE/'neutral_m3_sweep_angles_20261027'
ROOT=BASE/'neutral_m3_sweep_angles_fpp_20261028'
OLD=BASE/'75°/neutral_m3_frequency_limits_20261027'
REF=BASE/'75°/neutral_xbeta_four_omega_20261026/cases/bluntplate'
L=.035;U=1810.6307411853716;RHO=.018410817069110273;T=226.5
FREQ=[-240.2,-237.14,-200,-150,-112,-110,-100,-75,-50,-25,0,25,50,75,100,125,150,175,200,230,250,264.53,275,290,291.6]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    assert not (ROOT/'outputs/input_manifest.json').exists(),'Preserve completed preparation'
    for n in ['work','outputs','cases','runs','matrix_exports','native_pencils','inputs']:(ROOT/n).mkdir(parents=True,exist_ok=True)
    for p in (TRIAL/'work').glob('*.py'):shutil.copy2(p,ROOT/'work'/p.name)
    for name in ['seeds.py','precompute.py','slice_input.py']:shutil.copy2(HERE/name,ROOT/'work'/name)
    # Runtime modules resolve their own isolated directory, independent of date.
    for p in (ROOT/'work').glob('*.py'):
        s=p.read_text();s=s.replace("ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/neutral_m3_sweep_angles_20261027')","ROOT=Path(__file__).resolve().parent.parent")
        if p.name=='physical.py':s=s.replace('flow=np.zeros((400,9,6,2))','flow=np.zeros((400,9,q.shape[1],2))')
        if p.name=='precompute.py':s=s.replace("    assert (ROOT/'outputs/stencil_validation.json').exists()","    assert (ROOT/'outputs/stencil_validation.json').exists()")
        p.write_text(s)
    cases={};protected=[];selection=[]
    for angle,rel in [(30,'30°'),(45,'45°swptangle'),(60,'60°'),(75,'75°')]:
        src=WIN/rel/'1.1m/fpp_output';tag=f's{angle}';dest=ROOT/'cases'/tag;(dest/'in').mkdir(parents=True,exist_ok=True)
        summary=json.loads((src/'summary.json').read_text());refs=summary['reference'];meta=json.loads((src/'flow_vertical2d.variables.json').read_text())
        shape,x,y,ge=read_grid(src/'grid_vertical2d.grd');names,arrays,fe=read_flow(src/'flow_vertical2d.dat',shape)
        assert list(names)==meta['variables'];assert shape[1]==400
        # Physical intermediate, then canonical nondimensionalization matching 75-degree LST.
        nd=summary['nondimensional'];grid=np.zeros((shape[1],shape[0],2,2));flow=np.zeros((shape[1],shape[0],len(names),2))
        grid[:,:,0,0]=x.T*(refs['length'] if nd else 1.)/L;grid[:,:,1,0]=y.T*(refs['length'] if nd else 1.)/L
        scales=[refs[k] if nd else 1. for k in ['density','velocity','velocity','velocity','temperature','pressure']][:len(names)]
        canonical=[RHO,U,U,U,T,RHO*U*U][:len(names)]
        for k,a in enumerate(arrays):flow[:,:,k,0]=a.T*scales[k]/canonical[k]
        assert np.isfinite(grid).all() and np.isfinite(flow).all() and np.min(flow[:,:,0,0])>0 and np.min(flow[:,:,4,0])>0
        files=[]
        for name,key,a in [('grid_vertical2d_N.h5','grid',grid),('flow_vertical2d_N.h5','baseflow',flow)]:
            p=dest/'in'/name
            with h5py.File(p,'w') as h:h.create_dataset(key,data=a)
            with h5py.File(p) as h:assert np.array_equal(a,h[key][...])
            files.append(dict(path=str(p),sha256=sha(p)))
        cfg=configparser.ConfigParser();cfg.read(REF/'config.ini');cfg['filenames']['dir']=str(dest/'in')+'/';cfg['domain'].update({'in':str(shape[0]),'jn':'400','iloc':'119','istart':'1','iend':str(shape[0])})
        with (dest/'config.ini').open('w') as f:cfg.write(f)
        for p in [src/'grid_vertical2d.grd',src/'flow_vertical2d.dat',src/'flow_vertical2d.variables.json',src/'summary.json']+[Path(r['path']) for r in files]+[dest/'config.ini']:protected.append(dict(path=str(p),sha256=sha(p)))
        cases[tag]=dict(name=tag,angle=angle,config=str(dest/'config.ini'),in_dir=str(dest/'in'),source_dir=str(src),stations=shape[0],normal_points=400,x_mm=(grid[0,:,0,0]*35).tolist(),wall_y_mm=(grid[0,:,1,0]*35).tolist(),flat_start_mm=0.,seed_station=119)
        selection.append(dict(angle=angle,path=str(src),source_nondimensional=nd,source_reference=refs,canonical_reference=dict(length=L,density=RHO,velocity=U,temperature=T),variables=list(names),shape=list(shape),conversion='Only reference-unit conversion and transposition; no profile interpolation',files=files))
        print('FPP',angle,shape,'source_nondimensional',nd,'variables',len(names),flush=True)
        if angle==75:
            errors={}
            for name,key,a in [('grid_vertical2d_N.h5','grid',grid),('flow_vertical2d_N.h5','baseflow',flow)]:
                with h5py.File(REF/'in'/name) as h:b=h[key][...]
                errors[key]=float(np.max(abs(a-b)))
            selection[-1]['previous75_anchor_max_absolute_errors']=errors
            selection[-1]['previous75_curves_reusable']=max(errors.values())<2e-12
    protected.append(dict(path=str(OLD/'outputs/input_manifest.json'),sha256=sha(OLD/'outputs/input_manifest.json')))
    (ROOT/'outputs/cases.json').write_text(json.dumps(cases,indent=2)+'\n');(ROOT/'outputs/input_manifest.json').write_text(json.dumps(protected,indent=2)+'\n')
    (ROOT/'outputs/input_selection.json').write_text(json.dumps(dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),selected=selection,pending_angles=[50,55,63,65,68,70,73,80],rule='User supplied fpp directory location and will supplement missing files; dense-grid status by 1.1m directory.'),indent=2)+'\n')
    (ROOT/'outputs/frequency_conversion.json').write_text(json.dumps(dict(Lref_m=L,Uref_m_s=U,definition='omega_nd=2*pi*f_kHz*1000*Lref/Uref; all fields converted to same canonical reference',frequencies=[dict(f_kHz=f,omega_nd=float(2*np.pi*1000*f*L/U)) for f in FREQ]),indent=2)+'\n')
    # Previous native-stencil verification remains applicable to unchanged exporter;
    # additionally verify this stage against its full FPP fields before preparation.
    shutil.copy2(HERE/'verify_stencils.py',ROOT/'work/verify_stencils.py')
    p=ROOT/'work/verify_stencils.py';s=p.read_text().replace("ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/neutral_m3_sweep_angles_20261027')","ROOT=Path(__file__).resolve().parent.parent").replace("[('s73',119),('s73',241),('s73',400),('s73',2),('s60',119)]","[('s30',119),('s45',241),('s60',400),('s60',2),('s75',119)]");p.write_text(s)
    print('PREPARED',ROOT,len(protected),flush=True)
if __name__=='__main__':main()
