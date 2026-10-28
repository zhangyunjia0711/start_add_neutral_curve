"""Add only newly available exact 1.1m FPP input fields; preserve existing cases."""
import configparser,datetime,hashlib,json,sys
from pathlib import Path
import h5py,numpy as np
WIN=Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate')
ROOT=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/neutral_m3_sweep_consistent_20261028')
sys.path.insert(0,str(WIN/'work'))
from convert_plot3d_2d_to_h5 import read_grid,read_flow
L=.035;U=1810.6307411853716;RHO=.018410817069110273;T=226.5
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):
    q=p.with_suffix('.add.tmp');q.write_text(json.dumps(d,indent=2)+'\n');q.replace(p)
def main():
    cases=json.loads((ROOT/'outputs/cases.json').read_text());protected=json.loads((ROOT/'outputs/input_manifest.json').read_text());selection=json.loads((ROOT/'outputs/input_selection.json').read_text());added=[]
    for angle in [50,55,65,70,73]:
        tag=f's{angle}'
        if tag in cases:continue
        src=WIN/(str(angle)+'°')/'1.1m/fpp_output';dest=ROOT/'cases'/tag;assert not dest.exists();(dest/'in').mkdir(parents=True)
        sources=[src/n for n in ['grid_vertical2d.grd','flow_vertical2d.dat','flow_vertical2d.variables.json','summary.json','fpp_ui_config.yaml']];before={str(p):sha(p) for p in sources}
        summary=json.loads((src/'summary.json').read_text());refs=summary['reference'];meta=json.loads((src/'flow_vertical2d.variables.json').read_text());shape,x,y,_=read_grid(src/'grid_vertical2d.grd');names,arrays,_=read_flow(src/'flow_vertical2d.dat',shape)
        assert list(names)==meta['variables'] and list(names)[:5]==['rho','u','v','w','T'] and shape[1]==400
        nd=summary['nondimensional'];grid=np.zeros((shape[1],shape[0],2,2));flow=np.zeros((shape[1],shape[0],len(names),2))
        grid[:,:,0,0]=x.T*(refs['length'] if nd else 1.)/L;grid[:,:,1,0]=y.T*(refs['length'] if nd else 1.)/L
        scales=[refs[k] if nd else 1. for k in ['density','velocity','velocity','velocity','temperature','pressure']][:len(names)];canonical=[RHO,U,U,U,T,RHO*U*U][:len(names)]
        for k,a in enumerate(arrays):flow[:,:,k,0]=a.T*scales[k]/canonical[k]
        assert np.isfinite(grid).all() and np.isfinite(flow).all() and flow[:,:,0,0].min()>0 and flow[:,:,4,0].min()>0
        assert before=={str(p):sha(p) for p in sources},'Input changed while being read'
        files=[]
        for name,key,a in [('grid_vertical2d_N.h5','grid',grid),('flow_vertical2d_N.h5','baseflow',flow)]:
            path=dest/'in'/name
            with h5py.File(path,'w') as h:h.create_dataset(key,data=a)
            with h5py.File(path) as h:assert np.array_equal(a,h[key][...])
            files.append(dict(path=str(path),sha256=sha(path)))
        wallx=grid[0,:,0,0]*35;seed=int(np.argmin(abs(wallx+25.070568019862893)))+1
        cfg=configparser.ConfigParser();cfg.read(ROOT/'cases/s60/config.ini');cfg['filenames']['dir']=str(dest/'in')+'/';cfg['domain'].update({'in':str(shape[0]),'jn':'400','iloc':str(seed),'istart':'1','iend':str(shape[0])})
        with (dest/'config.ini').open('w') as f:cfg.write(f)
        protected.extend(dict(path=path,sha256=h) for path,h in before.items());protected.extend(files);protected.append(dict(path=str(dest/'config.ini'),sha256=sha(dest/'config.ini')))
        cases[tag]=dict(name=tag,angle=angle,config=str(dest/'config.ini'),in_dir=str(dest/'in'),source_dir=str(src),stations=shape[0],normal_points=400,x_mm=wallx.tolist(),wall_y_mm=(grid[0,:,1,0]*35).tolist(),flat_start_mm=0.,seed_station=seed,normal_derivative_mode='wall_local_consistent')
        rec=dict(angle=angle,path=str(src),source_nondimensional=nd,source_reference=refs,canonical_reference=dict(length=L,density=RHO,velocity=U,temperature=T),variables=list(names),shape=list(shape),x_range_mm=[float(wallx.min()),float(wallx.max())],conversion='Only reference-unit conversion and transposition; no profile interpolation',files=files,source_hashes=before,source_files_stable_during_read=True)
        selection['selected'].append(rec);added.append(rec);print('ADDED_FPP',angle,shape,'seed',seed,'x_range_mm',rec['x_range_mm'],flush=True)
    selection['pending_angles']=[63,68,80];selection['additional_input_check_time']=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
    write(ROOT/'outputs/cases.json',cases);write(ROOT/'outputs/input_manifest.json',protected);write(ROOT/'outputs/input_selection.json',selection);write(ROOT/'outputs/additional_fpp_inputs.json',added)
if __name__=='__main__':main()
