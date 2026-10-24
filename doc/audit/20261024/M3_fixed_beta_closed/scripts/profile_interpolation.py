"""Documented, adjacent-CFD local-profile interpolation for neutral turning points.

Interpolated physical x is a model coordinate, never an extra CFD grid station.
Every endpoint operator must reproduce the actual native operator and vector.
"""
import configparser,csv,hashlib,json,re,struct
from pathlib import Path
import h5py,numpy as np
from common import ROOT,PREVIOUS,SOLVER,SOLVER_HASH,OPTIONS,run,table,sha
from trace_fixed_neutral import field

FLOAT=r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?'
SOURCE=ROOT.parent/'M3/in'
AXIS={int(q['iloc']):q for q in json.loads((ROOT/'outputs/prior_axis_seeds.json').read_text())['records']}

def mac(a,b):return float(abs(np.vdot(a,b))**2/(np.vdot(a,a).real*np.vdot(b,b).real))

def vector(path):
    data=path.read_bytes();offset=0;out=[]
    while offset<len(data):
        cls,n=struct.unpack_from('>ii',data,offset);offset+=8;assert cls==1211214 and n==2000
        end=offset+16*n;out.append(np.frombuffer(data[offset:end],dtype='>c16').astype(complex));offset=end
    assert offset==len(data);return out

def evaluate(mode,station,state,omega,h,label):
    spec=dict(station=station,mode=str(mode),seed={k:float(state[k]) for k in ['omega_r','beta_r','alpha_r','alpha_i']},
        seed_shift=1e-5,solver=str(SOLVER),solver_sha256=SOLVER_HASH,coordinates='wall_local',task=33,
        outer_parameter='iloc',outer_min=station,outer_max=station,outer_points=1,
        inner_parameter='omega',inner_min=omega-h,inner_max=omega+h,inner_points=3,
        sample_only=True,anchor_row_seed=True,step=.05,jump=4,adaptive=True,
        vector_tracking=True,minimum_overlap=.95,minimum_step=1e-6,retry_eigenpairs=8,
        options=['-pep_tol','1e-13','-pep_max_it','1000','-pep_nev','1','-pep_ncv','8',
                 '-st_pc_factor_mat_solver_type','superlu','-pep_view_vectors','binary:out/all_vectors.petsc:default:append'],timeout=120)
    token=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()[:12];spec['name']=label+'_'+token
    summary=run(spec);directory=ROOT/'runs'/spec['name']
    assert summary['exit']!=124 and not any(summary[k] for k in ['petsc_errors','unused_options','failed_field_records','coordinate_metadata_missing'])
    values=sorted(table(directory/'out/neutral_field.plt'),key=lambda q:q['omega_r']);assert len(values)==3
    assert all(q['status']==0 and q['beta_r']==-7 for q in values)
    log=(directory/'solver.log').read_text();pairs=[]
    for m in re.finditer(r'^\s*\d+\s+\(\s*('+FLOAT+r')\s*,\s*('+FLOAT+r')\s*\)\s+The error is\s+('+FLOAT+r')',log,re.M):
        pairs.append((complex(float(m[1].replace('D','E')),float(m[2].replace('D','E'))),float(m[3].replace('D','E'))))
    rawfile=directory/'out/all_vectors.petsc';qs=vector(rawfile);assert len(qs)==len(pairs)
    center=values[1];target=complex(center['alpha_r'],center['alpha_i'])
    k=min(range(len(pairs)),key=lambda k:abs(pairs[k][0]-target));assert abs(pairs[k][0]-target)<1e-8 and pairs[k][1]<1e-7
    matching=[j for j in range(len(pairs)) if abs(pairs[j][0]-target)<1e-8]
    assert all(mac(qs[k],qs[j])>1-1e-8 for j in matching),'Ambiguous eigenvector association'
    return dict(point=center,samples=values,g=center['growth'],g_omega=(values[2]['growth']-values[0]['growth'])/(2*h),
        g_omega_omega=(values[2]['growth']-2*center['growth']+values[0]['growth'])/(h*h),h=h,
        run=spec['name'],vector_file=str(rawfile),vector_index=k,polynomial_residual=pairs[k][1],q=qs[k])

class Pair:
    def __init__(self,left,right):
        assert right==left+1
        self.left,self.right=left,right;self.profiles={}
        for i in [left,right]:
            base=AXIS[i];rows,name=field(i,base,0.,1e-9,2,f'M3b7_basis_i{i:03d}')
            with (ROOT/'runs'/name/'wall_local_profile_metadata.csv').open() as f:
                meta=next(q for q in csv.DictReader(line for line in f if not line.startswith('#')) if int(q['iloc'])==i)
            with h5py.File(SOURCE/'grid_vertical2d_N.h5','r') as f:g=f['grid'][:,i-1,:,0]
            with h5py.File(SOURCE/'flow_vertical2d_N.h5','r') as f:q=f['baseflow'][:,i-1,:,0]
            t=np.array([float(meta['t_x']),float(meta['t_y'])]);normal=np.array([float(meta['n_x']),float(meta['n_y'])])
            n=(g-g[0])@normal;n[0]=0.;local=q.copy();local[:,1]=q[:,1:3]@t;local[:,2]=q[:,1:3]@normal
            assert np.all(np.diff(n)>0) and np.all(local[:,0]>0) and np.all(local[:,4]>0)
            self.profiles[i]=dict(n=n,q=local,x=float(g[0,0]),metadata=meta)
        flags=[self.profiles[i]['metadata']['identity_flat'].strip() for i in [left,right]]
        assert flags[0]==flags[1],'No interpolation across a discrete derivative-path change'
        self.flat=flags[0]=='T'

    def build(self,theta):
        assert 0<=theta<=1
        token=hashlib.sha256(float(theta).hex().encode()).hexdigest()[:12]
        directory=ROOT/'inputs'/f'pair{self.left:03d}_{self.right:03d}_n9_{token}'
        if (directory/'generation.json').exists():return directory
        (directory/'in').mkdir(parents=True,exist_ok=True)
        a,b=self.profiles[self.left],self.profiles[self.right]
        n=(1-theta)*a['n']+theta*b['n'];q=(1-theta)*a['q']+theta*b['q'];x=(1-theta)*a['x']+theta*b['x']
        width=max(1e-4,abs(b['x']-a['x'])/4);shear=0 if self.flat else 1e-8
        grid=np.zeros((400,9,2,2));grid[:,:,0,0]=x+width*np.arange(-4,5)[None,:]+shear*n[:,None]
        grid[:,:,1,0]=n[:,None];flow=np.zeros((400,9,6,2));flow[:,:,:,0]=q[:,None,:]
        cfg=configparser.ConfigParser();cfg.read(ROOT.parent/'M3/config.ini')
        cfg['domain']['in']='9';cfg['domain']['jn']='400';cfg['domain']['iloc']='5';cfg['domain']['istart']='1';cfg['domain']['iend']='9'
        with (directory/'config.ini').open('w') as f:cfg.write(f)
        files=[]
        for name,key,data in [('grid_vertical2d_N.h5','grid',grid),('flow_vertical2d_N.h5','baseflow',flow)]:
            with h5py.File(SOURCE/name,'r') as f:attrs=dict(f[key].attrs)
            path=directory/'in'/name
            with h5py.File(path,'w') as f:
                ds=f.create_dataset(key,data=data,compression='gzip',compression_opts=1)
                for k,v in attrs.items():ds.attrs[k]=v
                ds.attrs['plot3d_shape']=np.array([9,400],dtype=np.int64)
                ds.attrs['fortran_logical_shape']=np.array([2,data.shape[2],9,400],dtype=np.int64)
            files.append(dict(path=str(path),sha256=sha(path)))
        metadata=dict(left=self.left,right=self.right,theta=theta,x_N=x,x_mm=35*x,flat_derivative_path=self.flat,
            model='Linear interpolation of adjacent native frozen-local n and primitive q at matched j; not a CFD refinement',
            physical_x_formula='(1-theta)*x_native_left+theta*x_native_right',files=files,
            endpoint_metadata=[a['metadata'],b['metadata']])
        (directory/'generation.json').write_text(json.dumps(metadata,indent=2)+'\n');return directory

    def eval(self,theta,state,omega,h,label):
        directory=self.build(theta);r=evaluate(directory,5,state,omega,h,label)
        metadata=json.loads((directory/'generation.json').read_text());r.update(theta=theta,x_mm=metadata['x_mm'],generation=str(directory/'generation.json'))
        return r

    def endpoint_checks(self,state,omega,label):
        checks=[]
        for theta,i in [(0.,self.left),(1.,self.right)]:
            native=evaluate('M3',i,state,omega,1e-6,label+f'_native{i}')
            aux=self.eval(theta,native['point'],omega,1e-6,label+f'_endpoint{i}')
            difference=abs(complex(native['point']['alpha_r']-aux['point']['alpha_r'],native['point']['alpha_i']-aux['point']['alpha_i']))
            overlap=mac(native['q'],aux['q'])
            check=dict(native_station=i,theta=theta,alpha_difference=difference,raw_MAC=overlap,
                native_run=native['run'],interpolated_endpoint_run=aux['run'],passed=difference<1e-7 and overlap>1-1e-8)
            checks.append(check);assert check['passed'],check
        return checks
