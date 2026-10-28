"""Read-only, bounded inventory of original angle inputs (Python 3.8)."""
import hashlib,json,math,re,struct
from pathlib import Path
import numpy as np
from scipy.io import FortranFile
ROOT=Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate')
OUT=ROOT/'work/M3_sweep_angles_20261027'

def read(p,grid=False):
    with FortranFile(str(p),'r',header_dtype='>u4') as f:
        nb=f.read_ints('>i4');dims=f.read_ints('>i4');v=f.read_reals('>f8')
    assert nb.tolist()==[1],(p,nb)
    ni,nj=map(int,dims[:2]);assert v.size==ni*nj*(2 if grid else 1),(p,dims,v.size)
    return v.reshape((ni,nj,2),order='F') if grid else v.reshape((ni,nj),order='F')

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    rows=[];missing=[]
    for angle in sorted(ROOT.iterdir()):
        if not angle.is_dir() or not re.match(r'^\d',angle.name):continue
        if angle.name=='75°':
            missing.append(dict(directory=angle.name,note='Use existing verified baseline/optimized data separately.'));continue
        stack=[(angle,0)]
        while stack:
            d,depth=stack.pop()
            children=list(d.iterdir())
            if (d/'plot3dx.unf').is_file() and all((d/(v+'.plot3dq.unf')).is_file() for v in 'RUVWT'):
                try:
                    grid=read(d/'plot3dx.unf',True);flow={v:read(d/(v+'.plot3dq.unf')) for v in 'RUVWT'}
                    assert all(a.shape==grid.shape[:2] for a in flow.values())
                    x,y=grid[:,:,0],grid[:,:,1];i=int(np.argmax(x[:,0]));u,v,w=[float(flow[k][i,-1]) for k in 'UVW'];speed=math.sqrt(u*u+v*v+w*w)
                    wall=grid[:,0,:];t=np.gradient(wall,axis=0);t/=np.linalg.norm(t,axis=1)[:,None];curvature=np.gradient(np.unwrap(np.arctan2(t[:,1],t[:,0])))/np.linalg.norm(np.gradient(wall,axis=0),axis=1)
                    row=dict(path=str(d),relative=str(d.relative_to(ROOT)),ni=x.shape[0],nj=x.shape[1],wall_x_mm=[float(x[:,0].min()*1000),float(x[:,0].max()*1000)],wall_y_mm=[float(y[:,0].min()*1000),float(y[:,0].max()*1000)],farfield_uvw=[u,v,w],speed=speed,angle_from_farfield=math.degrees(math.atan2(w,math.hypot(u,v))),T=float(flow['T'][i,-1]),rho=float(flow['R'][i,-1]),max_wall_velocity=max(float(np.max(np.abs(flow[k][:,0]))) for k in 'UVW'),curvature_range=[float(curvature.min()),float(curvature.max())],source_files=[dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in [d/'plot3dx.unf']+[d/(k+'.plot3dq.unf') for k in 'RUVWT']])
                    rows.append(row)
                except Exception as e:rows.append(dict(path=str(d),error=repr(e)))
            if depth<4:
                stack.extend((c,depth+1) for c in children if c.is_dir() and c.name not in ['work','outputs','tmp','DNS','fpp_output','__pycache__'] and not c.name.startswith('neutral'))
    (OUT/'inventory.json').write_text(json.dumps(dict(cases=rows,notes=missing),indent=2)+'\n')
    for r in rows:
        print(r.get('relative',r['path']),r.get('ni'),r.get('nj'),'x=',r.get('wall_x_mm'),'angle=',round(r.get('angle_from_farfield',0),5),'speed=',round(r.get('speed',0),3),'error=',r.get('error',''),flush=True)
if __name__=='__main__':main()
