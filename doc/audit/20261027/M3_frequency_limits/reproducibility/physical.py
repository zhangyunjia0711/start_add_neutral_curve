"""Actual local-profile interpolation for sub-grid x extrema (not CFD refinement)."""
import configparser,csv,hashlib,json
from pathlib import Path
import h5py,numpy as np
from solver import Pencil,export,PREVIOUS
from common import ROOT,CASES,sha

class Pair:
    def __init__(self,case,i,omega,tag):
        self.case,self.i,self.omega,self.tag=case,i,omega,tag;self.source=ROOT/'cases'/case/'in';self.profiles=[]
        for station in [i,i+1]:
            prior=json.loads((PREVIOUS/'native_pencils'/case/f'i{station:03d}'/'coefficients.json').read_text());directory=Path(prior['exports'][0]['files'][0]['path']).parent.parent
            with (directory/'wall_local_profile_metadata.csv').open() as f:
                meta=next(q for q in csv.DictReader(line for line in f if not line.startswith('#')) if int(q['iloc'])==station)
            with h5py.File(self.source/'grid_vertical2d_N.h5') as f:g=f['grid'][:,station-1,:,0]
            with h5py.File(self.source/'flow_vertical2d_N.h5') as f:q=f['baseflow'][:,station-1,:,0]
            t=np.array([float(meta['t_x']),float(meta['t_y'])]);normal=np.array([float(meta['n_x']),float(meta['n_y'])])
            n=(g-g[0])@normal;n[0]=0.;local=q.copy();local[:,1]=q[:,1:3]@t;local[:,2]=q[:,1:3]@normal
            assert np.all(np.diff(n)>0)
            self.profiles.append(dict(n=n,q=local,x=float(g[0,0]),metadata=meta))
        flags=[q['metadata']['identity_flat'].strip() for q in self.profiles];assert flags[0]==flags[1]
        self.flat=flags[0]=='T';self.cache={}
    def build(self,theta):
        assert 0<=theta<=1
        token=hashlib.sha256(float(theta).hex().encode()).hexdigest()[:12]
        directory=ROOT/'inputs'/f'{self.tag}_pair{self.i:03d}_{token}'
        if (directory/'generation.json').exists():return directory
        (directory/'in').mkdir(parents=True,exist_ok=True)
        a,b=self.profiles;n=(1-theta)*a['n']+theta*b['n'];q=(1-theta)*a['q']+theta*b['q'];x=(1-theta)*a['x']+theta*b['x']
        width=max(1e-4,abs(b['x']-a['x'])/4);shear=0 if self.flat else 1e-8
        grid=np.zeros((400,9,2,2));grid[:,:,0,0]=x+width*np.arange(-4,5)[None,:]+shear*n[:,None];grid[:,:,1,0]=n[:,None]
        flow=np.zeros((400,9,6,2));flow[:,:,:,0]=q[:,None,:]
        cfg=configparser.ConfigParser();cfg.read(ROOT/'cases'/self.case/'config.ini')
        cfg['domain'].update({'in':'9','jn':'400','iloc':'5','istart':'1','iend':'9'})
        with (directory/'config.ini').open('w') as f:cfg.write(f)
        files=[]
        for name,key,data in [('grid_vertical2d_N.h5','grid',grid),('flow_vertical2d_N.h5','baseflow',flow)]:
            with h5py.File(self.source/name) as f:attrs=dict(f[key].attrs)
            path=directory/'in'/name
            with h5py.File(path,'w') as f:
                ds=f.create_dataset(key,data=data,compression='gzip',compression_opts=1)
                for k,v in attrs.items():ds.attrs[k]=v
                ds.attrs['plot3d_shape']=np.array([9,400],dtype=np.int64);ds.attrs['fortran_logical_shape']=np.array([2,data.shape[2],9,400],dtype=np.int64)
            files.append(dict(path=str(path),sha256=sha(path)))
        (directory/'generation.json').write_text(json.dumps(dict(case=self.case,left=self.i,right=self.i+1,theta=theta,x_mm=35*x,model='Adjacent native local n and primitive q interpolation at matched j; not a new CFD station',files=files),indent=2)+'\n')
        return directory
    def pencil(self,theta):
        key=float(theta).hex()
        if key in self.cache:return self.cache[key]
        directory=self.build(theta);H=16.
        minus,rm=export(directory,5,-H,0.,label=f'{self.case}_physical');zero,r0=export(directory,5,0.,0.,label=f'{self.case}_physical');plus,rp=export(directory,5,H,0.,label=f'{self.case}_physical')
        frequency,rf=export(directory,5,0.,1.,label=f'{self.tag}_physical_frequency');p=object.__new__(Pencil);p.case=str(directory);p.i=5;p.omega=self.omega;p.G=frequency[0]-zero[0]
        p.C=[[zero[k],(plus[k]-minus[k])/(2*H),(plus[k]+minus[k]-2*zero[k])/(2*H*H)] for k in range(3)]
        p.baseC=[list(row) for row in p.C];p.C[0][0]=p.C[0][0]+self.omega*p.G;p.report=dict(exports=[rm,r0,rp,rf]);self.cache[key]=p;return p
