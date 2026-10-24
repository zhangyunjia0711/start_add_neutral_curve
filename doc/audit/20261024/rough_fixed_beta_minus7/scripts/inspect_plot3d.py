"""Read-only inspection of the supplied rough-case postprocessed Plot3D inputs."""
import json,struct,hashlib
from pathlib import Path
import numpy as np,h5py

CASE=Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate/75°/linear_curvature/theta1_0_f_1_over_2/rough')
OLD=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/bluntplate/M3/in')

def records(path,endian='<'):
    raw=path.read_bytes();offset=0;rows=[]
    while offset<len(raw):
        n=struct.unpack_from(endian+'I',raw,offset)[0];offset+=4
        assert offset+n+4<=len(raw),(path,offset,n,len(raw))
        rows.append(raw[offset:offset+n]);offset+=n
        assert struct.unpack_from(endian+'I',raw,offset)[0]==n;offset+=4
    return rows

def gridflow():
    g=records(CASE/'fpp_output/grid_vertical2d.grd');f=records(CASE/'fpp_output/flow_vertical2d.dat')
    print('RECORD_SIZES',list(map(len,g)),list(map(len,f)))
    for r in g[:2]:print('GRID_HEADER',np.frombuffer(r,dtype='<i4').tolist())
    for r in f[:2]:print('FLOW_HEADER',np.frombuffer(r,dtype='<i4').tolist())
    n=361*400
    ga=np.frombuffer(g[-1],dtype='<f8').reshape((2,400,361)).transpose(1,2,0)
    fa=np.frombuffer(f[-1],dtype='<f8').reshape((6,400,361)).transpose(1,2,0)
    assert ga.shape==(400,361,2) and fa.shape==(400,361,6)
    normal=np.linalg.norm(ga-ga[0][None,:,:],axis=2)
    wall=ga[0];t=np.gradient(wall,axis=0);t/=np.linalg.norm(t,axis=1)[:,None]
    localn=np.stack([-t[:,1],t[:,0]],axis=1)
    if np.median(np.sum((ga[1]-wall)*localn,axis=1))<0:localn=-localn
    proj=np.einsum('jik,ik->ji',ga-wall[None,:,:],localn)
    report=dict(grid_shape=list(ga.shape),flow_shape=list(fa.shape),x_wall_range_m=[float(wall[:,0].min()),float(wall[:,0].max())],
        x_monotone=bool(np.all(np.diff(wall[:,0])>0)),minimum_normal_increment_m=float(np.diff(proj,axis=0).min()),
        primitive_min=fa.min(axis=(0,1)).tolist(),primitive_max=fa.max(axis=(0,1)).tolist(),
        farfield_first_middle_last=fa[-1,[0,180,360]].tolist(),wall_first_middle_last=wall[[0,180,360]].tolist(),
        source_files=[dict(path=str(CASE/'fpp_output'/name),sha256=hashlib.sha256((CASE/'fpp_output'/name).read_bytes()).hexdigest()) for name in ['grid_vertical2d.grd','flow_vertical2d.dat','flow_vertical2d.variables.json','summary.json','fpp_ui_config.yaml']])
    (CASE/'outputs/input_inspection.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(report,indent=2))
    for name,key in [('grid_vertical2d_N.h5','grid'),('flow_vertical2d_N.h5','baseflow')]:
        with h5py.File(OLD/name) as h:
            a=h[key];print('OLD_H5',name,a.shape,{k:np.array(v).tolist() for k,v in a.attrs.items()})
            print('OLD_VALUES',a[-1,[0,118,840],:,0].tolist())

if __name__=='__main__':gridflow()
