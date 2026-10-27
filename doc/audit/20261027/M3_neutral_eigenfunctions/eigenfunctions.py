"""Inspect archived spatial M3 eigenvectors on both neutral boundaries.

Run under WSL Python. Input files are read only; all outputs are task scoped.
"""
import csv, datetime, hashlib, json, re, shutil, sys
sys.dont_write_bytecode=True
from pathlib import Path
import h5py
import numpy as np
from scipy.sparse.linalg import norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

BASE=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°')
WIN=Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate/75°')
SOURCE=BASE/'neutral_xbeta_four_omega_20261026'
CONNECTION=BASE/'neutral_frequency_connection_20261026_f100'
DEST=WIN/'outputs/M3_neutral_eigenfunctions_20261027'
OVERLAY=WIN/'outputs/M3_xbeta_frequency_overlay_20261027'
sys.path.insert(0,str(SOURCE/'work'))
from solver import Pencil, PREVIOUS
from common import snapshot, SOLVER, SOLVER_HASH

FS=[-100,-50,0,50,100]
CASES=['bluntplate','rough']
SIDES=['lower','upper']
COMPS=['rho','u','v','w','T']
LABELS=[r'$|\hat\rho|/S$',r'$|\hat u_t|/S$',r'$|\hat v_n|/S$',r'$|\hat w|/S$',r'$|\hat T|/S$']
COLORS={-100:'#762a83',-50:'#2166ac',0:'#222222',50:'#1b9e77',100:'#d73027'}
TEC_COLORS={-100:'PURPLE',-50:'BLUE',0:'BLACK',50:'GREEN',100:'RED'}
NAMES=['M3_neutral_eigenfunctions','bluntplate_M3_neutral_eigenfunctions','rough_M3_neutral_eigenfunctions']
VARS=['n_mm','f_kHz','beta','alpha_r','alpha_i','x_mm','case_id','boundary_id']
VARS += [c+'_'+part for c in COMPS for part in ['abs','real','imag']]

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def js(p,obj):p.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n')
def complex_interp(x,n,q):return np.interp(x,n,q.real)+1j*np.interp(x,n,q.imag)
def integral(a,n):return np.trapz(a,x=n,axis=0)
def vmac(q,r,n):
    z=integral(np.sum(q[:,1:4].conj()*r[:,1:4],axis=1),n)
    return float(abs(z)**2/(integral(np.sum(abs(q[:,1:4])**2,axis=1),n)*integral(np.sum(abs(r[:,1:4])**2,axis=1),n)))

def profile(case,station):
    prior=PREVIOUS/'native_pencils'/case/f'i{station:03d}'/'coefficients.json'
    a=json.loads(prior.read_text());run=Path(a['exports'][0]['files'][0]['path']).parent.parent
    path=run/'wall_local_profile_metadata.csv'
    with path.open() as f:meta=next(q for q in csv.DictReader(s for s in f if not s.startswith('#')) if int(q['iloc'])==station)
    grid=SOURCE/'cases'/case/'in/grid_vertical2d_N.h5'
    flow=SOURCE/'cases'/case/'in/flow_vertical2d_N.h5'
    with h5py.File(grid) as f:g=f['grid'][:,station-1,:,0]
    with h5py.File(flow) as f:base=f['baseflow'][:,station-1,:,0]
    normal=np.array([float(meta['n_x']),float(meta['n_y'])])
    tangent=np.array([float(meta['t_x']),float(meta['t_y'])])
    n=35*((g-g[0])@normal);n[0]=0
    local=base.copy();local[:,1]=base[:,1:3]@tangent;local[:,2]=base[:,1:3]@normal
    assert np.all(np.diff(n)>0) and len(n)==400
    return n,local,dict(metadata=meta,metadata_file=str(path),metadata_sha256=sha(path),grid_path=str(grid),grid_sha256=sha(grid),flow_path=str(flow),flow_sha256=sha(flow))

def collect():
    DEST.mkdir(parents=True,exist_ok=True)
    reportfile=CONNECTION/'outputs/connection_report.json'
    original=json.loads(reportfile.read_text());assert sha(SOLVER)==SOLVER_HASH==original['solver_sha256']
    protected=snapshot()
    modes=[];profiles={};sources=[];continuous=[]
    for p in original['paths']:
        case,side=p['case'],p['side'];station=p['station']
        if case not in profiles:profiles[case]=profile(case,station)
        n,base,prov=profiles[case]
        chosen=[]
        for f in FS:
            record=next(r for r in p['records'] if abs(r['f_kHz']-f)<1e-9)
            raw=np.load(record['vector_file']).reshape(-1,5)
            assert sha(record['vector_file'])==record['vector_sha256'] and raw.shape==(400,5) and np.isfinite(raw).all()
            scale=float(np.max(np.linalg.norm(raw[:,1:4],axis=1)));assert scale>0
            raw_scaled=raw/scale
            csvfile=OVERLAY/f'{case}_f{f:+d}kHz.csv'
            with csvfile.open() as stream:pts=[r for r in csv.DictReader(stream) if float(r['point_kind'])==0 and float(r['station_coordinate'])==station]
            pts.sort(key=lambda r:float(r['beta']));assert len(pts)==2
            point=pts[SIDES.index(side)]
            betaerr=abs(float(point['beta'])-record['beta'])
            alpha=complex(record['alpha_r'],record['alpha_i'])
            alphaerr=abs(alpha-complex(float(point['alpha_r']),float(point['alpha_i'])))
            assert betaerr<1e-7 and alphaerr<1e-7
            pencil=Pencil(case,station,record['omega_nd']);K=pencil.spatial(record['beta'])
            q=raw.ravel();res=(K[0]@q)+alpha*(K[1]@q)+alpha**2*(K[2]@q)
            den=sum(abs(alpha)**k*norm(K[k]) for k in range(3))*np.linalg.norm(q)
            residual=float(np.linalg.norm(res)/den)
            assert residual<1e-9 and abs(alpha.imag)<1e-8
            check=next(r for r in p['original_checks'] if abs(r['f_kHz']-f)<1e-9)
            assert check['raw_MAC']>1-1e-9 and check['alpha_difference']<1e-7
            r=dict(case=case,side=side,station=station,f_kHz=f,x_mm=p['x_mm'],beta=record['beta'],omega_nd=record['omega_nd'],alpha_r=alpha.real,alpha_i=alpha.imag,
                   source_vector=record['vector_file'],source_sha256=record['vector_sha256'],velocity_scale_S=scale,original_frozen_check=check,
                   overlay_path=str(csvfile),overlay_sha256=sha(csvfile),overlay_beta_error=betaerr,overlay_alpha_error=alphaerr,
                   matrix_scaled_residual=residual,source_right_residual=record['right_residual'],raw=raw,unphased=raw_scaled,n=n)
            modes.append(r);chosen.append(r)
        q0=next(r for r in chosen if r['f_kHz']==0)['unphased']
        for r in chosen:r['velocity_weighted_MAC_to_f0']=vmac(q0,r['unphased'],n)
        # Continuous metrics use the same archived neutral path, not interpolation between frequencies.
        for r in p['records']:
            raw=np.load(r['vector_file']).reshape(400,5)
            assert sha(r['vector_file'])==r['vector_sha256']
            S=float(np.max(np.linalg.norm(raw[:,1:4],axis=1)));q=raw/S
            continuous.append(metrics(dict(case=case,side=side,f_kHz=r['f_kHz'],x_mm=p['x_mm'],beta=r['beta'],alpha_r=r['alpha_r'],alpha_i=r['alpha_i'],n=n,q=q,velocity_weighted_MAC_to_f0=vmac(q0,q,n))))
        sources.append(dict(case=case,side=side,station=station,x_mm=p['x_mm'],profile=prov,continuous_points=len(p['records']),minimum_neighbor_raw_MAC=p['minimum_neighbor_MAC']))
    assert len(modes)==20
    # A single physical w probe fixes the arbitrary overall phase for all modes.
    probe_grid=np.linspace(0,min(r['n'][-1] for r in modes),10001)
    min_w=np.min([abs(complex_interp(probe_grid,r['n'],r['unphased'][:,3])) for r in modes],axis=0)
    probe=float(probe_grid[np.argmax(min_w)]);strength=float(np.max(min_w));assert strength>0.02
    for r in modes:
        w=complex_interp(probe,r['n'],r['unphased'][:,3]);angle=float(np.angle(w))
        r.update(phase_reference_n_mm=probe,phase_rotation_rad=-angle,q=r['unphased']*np.exp(-1j*angle))
        z=complex_interp(probe,r['n'],r['q'][:,3]);assert z.real>0 and abs(z.imag)<1e-12
        assert np.max(abs(abs(r['q'])-abs(r['unphased'])))<1e-12
        # Directly check invariance under an arbitrary complex rescaling.
        a=r['raw']*(2.7*np.exp(1.234j));S=np.max(np.linalg.norm(a[:,1:4],axis=1));a/=S
        a*=np.exp(-1j*np.angle(complex_interp(probe,r['n'],a[:,3])))
        assert np.max(abs(a-r['q']))<1e-11
        r.update(metrics(r))
    matrixfiles={}
    for case in CASES:
        i=next(r['station'] for r in modes if r['case']==case)
        for path in sorted((PREVIOUS/'native_pencils'/case/f'i{i:03d}').glob('C*.npz')):matrixfiles[str(path)]=sha(path)
        path=SOURCE/'native_pencils'/case/f'i{i:03d}'/'G.npz';matrixfiles[str(path)]=sha(path)
    src=BASE/'bluntplate/neutral_beta_full_20261022/work/solver_source/src/lst'
    refs=[]
    for name in ['module_lst.f90','lst_eqn.f90','lst_solver_abstract.f90']:
        file=src/name;target=DEST/'source_evidence'/name;target.parent.mkdir(exist_ok=True);shutil.copy2(file,target)
        refs.append(dict(path=str(file),sha256=sha(file)))
    boundary_checks=dict(wall_velocity_and_temperature_max_over_S=float(max(abs(r['q'][0,1:]).max() for r in modes)),farfield_velocity_and_temperature_max_over_S=float(max(abs(r['q'][-1,1:]).max() for r in modes)))
    assert boundary_checks['wall_velocity_and_temperature_max_over_S']<1e-12 and boundary_checks['farfield_velocity_and_temperature_max_over_S']<1e-12
    evidence=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),source_report=str(reportfile),source_report_sha256=sha(reportfile),boundary_checks=boundary_checks,
       solver_sha256=SOLVER_HASH,vector_order=['rho','u_t','v_n','w','T'],length_reference_m=.035,
       amplitude_normalization='Divide the entire complex vector by S=max_n sqrt(|u_t|^2+|v_n|^2+|w|^2); primitive variables retain solver nondimensional reference units.',
       phase_normalization='Rotate the entire complex vector so interpolated w at one common physical n is positive real; choose n by maximin |w|/S across all 20 modes.',
       phase_reference_n_mm=probe,minimum_probe_w_over_S=strength,source_code=refs,coefficient_hashes=matrixfiles,
       source_paths=sources,protected_inputs_before=protected,
       records=[{k:v for k,v in r.items() if k not in ['raw','unphased','q','n']} for r in modes],
       limitations=['Each case uses one common native streamwise station, not an entire streamwise eigenfunction sweep.','The two cases have slightly different x stations.','Frequency, neutral beta and alpha change together.','Linear eigenvectors give relative shapes, not absolute physical disturbance amplitudes.','Continuous 139-point paths supply shape metrics; original matrix residuals are rechecked for the 20 plotted frequency checkpoints.'])
    js(DEST/'eigenfunction_verification.json',evidence)
    data={}
    for r in modes:
        key=f"{r['case']}_{r['side']}_f{r['f_kHz']:+d}"
        data[key+'_raw']=r['raw'];data[key+'_normalized']=r['q']
    for case,(n,base,_) in profiles.items():data[case+'_n_mm']=n;data[case+'_local_baseflow']=base
    np.savez_compressed(DEST/'eigenfunctions_original_and_normalized.npz',**data)
    writecsv(DEST/'frequency_evolution_metrics.csv',continuous)
    writecsv(DEST/'selected_mode_metrics.csv',[metrics(r) for r in modes])
    print('COLLECTED',len(modes),'vectors; phase n_mm',probe,'minimum |w|/S',strength,'max residual',max(r['matrix_scaled_residual'] for r in modes),flush=True)
    return modes,continuous,evidence

def metrics(r):
    q=r['q'];n=r['n'];energy=np.sum(abs(q[:,1:4])**2,axis=1)
    m={k:r[k] for k in ['case','side','f_kHz','x_mm','beta','alpha_r','alpha_i','velocity_weighted_MAC_to_f0']}
    m['velocity_squared_centroid_n_mm']=float(integral(n*energy,n)/integral(energy,n))
    for j,c in enumerate(COMPS):
        peak=int(np.argmax(abs(q[:,j])));m[c+'_peak_over_S']=float(abs(q[peak,j]));m[c+'_peak_n_mm']=float(n[peak])
    return m

def writecsv(path,records):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)

def extent(modes):
    limits=[]
    for r in modes:
        q=abs(r['q']);mask=np.any(q>q.max(axis=0)[None,:]*1e-3,axis=1)
        limits.append(r['n'][np.where(mask)[0][-1]])
    return min(max(r['n'][-1] for r in modes),max(limits)*1.1)

def save(fig,name):
    for ext in ['png','svg']:fig.savefig(DEST/(name+'.'+ext),dpi=200,facecolor='white')
    plt.close(fig)

def plot(modes,continuous,evidence):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
    handles=[Line2D([0],[0],color=COLORS[f],lw=2,label=f'{f:+d} kHz' if f else '0 kHz') for f in FS]
    full_lim=extent(modes);lim=5.;evidence['plot_detail_n_max_mm']=lim;evidence['plot_full_n_max_mm']=full_lim
    for case in CASES:
        subset=[r for r in modes if r['case']==case];x=subset[0]['x_mm']
        fig,axes=plt.subplots(2,5,figsize=(16.2,8));fig.subplots_adjust(left=.06,right=.985,top=.83,bottom=.14,hspace=.55,wspace=.34)
        for row,side in enumerate(SIDES):
            group=sorted([r for r in subset if r['side']==side],key=lambda r:r['f_kHz'])
            for col in range(5):
                ax=axes[row,col]
                for r in group:ax.plot(r['n'],abs(r['q'][:,col]),color=COLORS[r['f_kHz']],lw=1.55)
                ax.set_xlim(0,lim);ax.set_xlabel('Wall-normal n (mm)');ax.set_ylabel(LABELS[col]);ax.grid(alpha=.2)
                if col==0:ax.text(0,1.18,('Lower' if row==0 else 'Upper')+' beta neutral boundary',transform=ax.transAxes,fontweight='bold')
            values=';  '.join(f"{r['f_kHz']:+d}: {r['beta']:.4f}" for r in group)
            axes[row,2].text(.5,1.18,'f (kHz): beta = '+values,transform=axes[row,2].transAxes,ha='center',fontsize=9)
        fig.suptitle(f'M3 eigenfunction magnitudes | {case} | x = {x:.5f} mm',y=.965,fontsize=15)
        fig.text(.5,.90,'Same local station and representative family; each full vector uses one velocity normalization S.',ha='center')
        fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.045),ncol=5,frameon=False)
        fig.text(.5,.014,'u_t: local wall tangent; v_n: local wall normal; w: spanwise. Density and temperature retain their nondimensional reference units.',ha='center',fontsize=9)
        for ext in ['png','svg']:fig.savefig(DEST/(case+'_eigenfunction_magnitudes.'+ext),dpi=200,facecolor='white')
        for ax in axes.flat:ax.set_xlim(0,full_lim)
        save(fig,case+'_eigenfunction_magnitudes_full_tail')
        # Full complex components, same phase gauge for all five frequencies.
        for side in SIDES:
            fig,axes=plt.subplots(2,5,figsize=(16.2,7.1));fig.subplots_adjust(left=.06,right=.985,top=.83,bottom=.16,hspace=.35,wspace=.34)
            for col,c in enumerate(COMPS):
                for row,part in enumerate(['real','imag']):
                    ax=axes[row,col]
                    for r in subset:
                        if r['side']==side:ax.plot(r['n'],getattr(r['q'][:,col],part),color=COLORS[r['f_kHz']],lw=1.5)
                    ax.set_xlim(0,lim);ax.grid(alpha=.2);ax.axhline(0,color='.5',lw=.5)
                    ax.set_xlabel('Wall-normal n (mm)');ax.set_ylabel(f'{part}({c}) / S')
            fig.suptitle(f'M3 complex eigenfunctions | {case}, {side} beta boundary | x = {x:.5f} mm',y=.965,fontsize=14)
            fig.text(.5,.90,f"Common phase: w(n = {evidence['phase_reference_n_mm']:.5f} mm) is positive real for every mode.",ha='center')
            fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.045),ncol=5,frameon=False)
            save(fig,case+'_'+side+'_complex_eigenfunctions')
    fig,axes=plt.subplots(2,2,figsize=(12.5,8.2));fig.subplots_adjust(left=.08,right=.98,bottom=.16,top=.86,hspace=.32,wspace=.23)
    for row,side in enumerate(SIDES):
        for col,case in enumerate(CASES):
            ax=axes[row,col]
            for r in modes:
                if (r['case'],r['side'])==(case,side):ax.plot(r['n'],abs(r['q'][:,3]),color=COLORS[r['f_kHz']],lw=1.8)
            ax.set(xlim=(0,lim),xlabel='Wall-normal n (mm)',ylabel=r'$|\hat w|/S$',title=f'{case}: {side} beta boundary');ax.grid(alpha=.2)
    fig.suptitle('M3 spanwise eigenfunction evolution on the neutral boundaries',y=.96,fontsize=14)
    fig.text(.5,.90,'bluntplate: x = -25.07057 mm; rough: x = -25.02701 mm',ha='center')
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.045),ncol=5,frameon=False)
    for ext in ['png','svg']:fig.savefig(DEST/('M3_eigenfunction_w_comparison.'+ext),dpi=200,facecolor='white')
    for ax in axes.flat:ax.set_xlim(0,full_lim)
    save(fig,'M3_eigenfunction_w_comparison_full_tail')
    fig,axes=plt.subplots(2,3,figsize=(14,8));fig.subplots_adjust(left=.065,right=.98,bottom=.17,top=.86,hspace=.37,wspace=.25)
    fields=[('w_peak_n_mm','Native-grid peak location of |w| (mm)'),('velocity_squared_centroid_n_mm','Velocity-squared centroid (mm)'),('velocity_weighted_MAC_to_f0','Velocity shape MAC relative to f = 0')]
    for row,case in enumerate(CASES):
        for col,(key,label) in enumerate(fields):
            ax=axes[row,col]
            for side,style in [('lower','-'),('upper','--')]:
                records=sorted([r for r in continuous if (r['case'],r['side'])==(case,side)],key=lambda r:r['f_kHz'])
                ax.plot([r['f_kHz'] for r in records],[r[key] for r in records],style,color='#3b6c9e' if side=='lower' else '#d47b24',label=side+' beta boundary')
                for r in records:
                    if r['f_kHz'] in FS:ax.plot(r['f_kHz'],r[key],'o',ms=4,color=COLORS[int(r['f_kHz'])])
            ax.set(xlabel='Physical frequency f (kHz)',ylabel=label,title=case,xlim=(-100,100));ax.grid(alpha=.2)
            if col==2:ax.set_ylim(0,1.03)
    fig.suptitle('Shape changes along the continuously tracked neutral family',y=.96,fontsize=14)
    fig.text(.5,.90,'139 archived frequency points per boundary; beta and alpha vary. Peak steps reflect the native n grid.',ha='center')
    axes[0,0].legend(frameon=False,fontsize=9)
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.045),ncol=5,frameon=False)
    save(fig,'M3_eigenfunction_frequency_metrics')
    # Locate all selected eigenfunctions directly on the preceding x-beta overlay.
    fig,axes=plt.subplots(1,2,figsize=(12.5,5.5));fig.subplots_adjust(left=.07,right=.98,bottom=.22,top=.83,wspace=.24)
    for ax,case in zip(axes,CASES):
        for f in FS:
            a=np.genfromtxt(OVERLAY/f'{case}_f{f:+d}kHz.csv',delimiter=',',names=True)
            ax.plot(a['x_mm'],a['beta'],color=COLORS[f],lw=1.4)
            for r in modes:
                if r['case']==case and r['f_kHz']==f:ax.plot(r['x_mm'],r['beta'],'o' if r['side']=='lower' else 's',color=COLORS[f],ms=5)
        x=next(r['x_mm'] for r in modes if r['case']==case);ax.axvline(x,color='.5',lw=.8,ls=':')
        ax.set(xlim=(-35,60),ylim=(-26,10),xlabel='Global x (mm)',ylabel='Neutral beta',title=f'{case}: sampled x = {x:.5f} mm');ax.grid(alpha=.2)
    fig.suptitle('Locations of the 20 eigenfunctions on the x-beta neutral curves',y=.95,fontsize=14)
    fig.text(.5,.88,'Circle: lower-beta root. Square: upper-beta root. Same representative as the five-frequency overlay.',ha='center',fontsize=10)
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.045),ncol=5,frameon=False)
    save(fig,'M3_eigenfunction_sampling_points')
    js(DEST/'eigenfunction_verification.json',evidence)

def tecplot(modes):
    ordered=sorted(modes,key=lambda r:(CASES.index(r['case']),SIDES.index(r['side']),FS.index(r['f_kHz'])))
    datasets={NAMES[0]:ordered,NAMES[1]:[r for r in ordered if r['case']=='bluntplate'],NAMES[2]:[r for r in ordered if r['case']=='rough']}
    zones=[]
    for name,group in datasets.items():
        with (DEST/(name+'.dat')).open('w') as f:
            f.write('TITLE="Normalized spatial M3 neutral eigenfunctions"\nVARIABLES='+','.join('"'+v+'"' for v in VARS)+'\n')
            for r in group:
                zone=f"{r['case']}_{r['side']}_f{r['f_kHz']:+d}kHz"
                arr=np.column_stack([r['n']]+[np.full(len(r['n']),r[k]) for k in ['f_kHz','beta','alpha_r','alpha_i','x_mm']]+[np.full(len(r['n']),CASES.index(r['case'])+1),np.full(len(r['n']),SIDES.index(r['side'])+1)]+[a for j in range(5) for a in [abs(r['q'][:,j]),r['q'][:,j].real,r['q'][:,j].imag]])
                assert arr.shape==(400,len(VARS)) and np.isfinite(arr).all()
                f.write(f'ZONE T="{zone}", I=400, DATAPACKING=POINT\nDT=('+','.join(['DOUBLE']*len(VARS))+')\n')
                np.savetxt(f,arr,fmt='%.17g')
                csvname=zone+'.csv'
                if name==NAMES[0]:np.savetxt(DEST/csvname,arr,delimiter=',',header=','.join(VARS),comments='',fmt='%.17g')
                zones.append(dict(dataset=name,zone=zone,points=400))
    js(DEST/'tecplot_zone_manifest.json',dict(variables=VARS,zones=zones,normalization='See eigenfunction_verification.json'))
    macro=['#!MC 1410']
    for name,group in datasets.items():
        macro+=['$!READDATASET \''+'"'+name+'.plt"'+'\'',' READDATAOPTION=NEW',' RESETSTYLE=YES',' INITIALPLOTTYPE=XYLINE','$!PLOTTYPE=XYLINE','$!DELETELINEMAPS']
        for j,r in enumerate(group,1):
            zone=f"{r['case']}_{r['side']}_f{r['f_kHz']:+d}kHz"
            line_style='SOLID' if r['side']=='lower' else 'DASHED'
            line_color=TEC_COLORS[r['f_kHz']]
            macro+=['$!CREATELINEMAP',f'$!LINEMAP [{j}]',f" NAME='{zone}'",f' ASSIGN{{ZONE={j} XAXISVAR=1 YAXISVAR=18}}',' LINES{SHOW=YES COLOR='+line_color+' LINEPATTERN='+line_style+'}',' SYMBOLS{SHOW=NO}']
        macro+=['$!ACTIVELINEMAPS=['+','.join(map(str,range(1,len(group)+1)))+']','$!VIEW FIT',f'$!SAVELAYOUT "{name}.lay"',f'$!WRITEDATASET "{name}_roundtrip.dat"',' BINARY=NO',' USEPOINTFORMAT=YES',' PRECISION=17',' INCLUDETEXT=NO',' INCLUDEGEOM=NO',' INCLUDECUSTOMLABELS=NO']
    macro+=['$!QUIT'];(DEST/'tecplot_roundtrip.mcr').write_text('\n'.join(macro)+'\n')

def heatmap():
    report=json.loads((CONNECTION/'outputs/connection_report.json').read_text());arrays={};provenance=[]
    fig,axes=plt.subplots(2,2,figsize=(12.5,8.4));fig.subplots_adjust(left=.075,right=.88,bottom=.12,top=.85,hspace=.31,wspace=.22)
    for p in report['paths']:
        case,side=p['case'],p['side'];n,_,_=profile(case,p['station'])
        recs=sorted(p['records'],key=lambda r:r['f_kHz']);freq=np.array([r['f_kHz'] for r in recs]);values=[]
        for r in recs:
            q=np.load(r['vector_file']).reshape(400,5);assert sha(r['vector_file'])==r['vector_sha256']
            S=np.max(np.linalg.norm(q[:,1:4],axis=1));values.append(abs(q[:,3])/S)
            provenance.append(dict(case=case,side=side,f_kHz=r['f_kHz'],path=r['vector_file'],sha256=r['vector_sha256']))
        key=case+'_'+side;arrays[key+'_f_kHz']=freq;arrays[key+'_n_mm']=n;arrays[key+'_w_abs_over_S']=np.array(values)
        ax=axes[SIDES.index(side),CASES.index(case)]
        im=ax.pcolormesh(freq,n,np.array(values).T,shading='auto',vmin=0,vmax=1,cmap='viridis',rasterized=True)
        ax.set(ylim=(0,5),xlim=(-100,100),xlabel='Physical frequency f (kHz)',ylabel='Wall-normal n (mm)',title=case+': '+side+' beta boundary')
        for f in FS:ax.axvline(f,color='white',lw=.6,alpha=.45,ls=':')
    cax=fig.add_axes([.91,.15,.018,.65]);fig.colorbar(im,cax=cax,label=r'$|\hat w|/S$')
    fig.suptitle('Continuous evolution of the spanwise eigenfunction on the neutral family',y=.96,fontsize=14)
    fig.text(.5,.90,'139 archived eigenvectors per boundary; same local x, changing neutral beta and alpha.',ha='center')
    save(fig,'M3_eigenfunction_frequency_heatmap')
    np.savez_compressed(DEST/'continuous_w_amplitude_maps.npz',**arrays)
    js(DEST/'continuous_vector_sources.json',dict(records=provenance))

def parse(path):
    variables=[];zones=[];inside=False
    for line in path.read_text(encoding='ascii').splitlines():
        s=line.strip()
        if s.startswith('VARIABLES'):inside=True;variables+=re.findall(r'"([^"]+)"',s);continue
        if inside and s.startswith('"'):variables+=re.findall(r'"([^"]+)"',s);continue
        inside=False
        if s.startswith('ZONE '):zones.append(dict(name=re.search(r'T\s*=\s*"([^"]+)"',s).group(1),rows=[],double=False))
        elif s.startswith('DT=') or s.startswith('DT ='):
            assert len(re.findall('DOUBLE',s))==len(variables);zones[-1]['double']=True
        elif zones and re.match(r'^[+-]?(?:[0-9]|\.[0-9])',s):
            row=[float(t) for t in s.split()];assert len(row)==len(variables) and np.isfinite(row).all();zones[-1]['rows'].append(row)
    assert all(z['double'] for z in zones);return variables,zones

def verify_native():
    assert 'ERROR' not in (DEST/'batch.log').read_text().upper();checks=[]
    for name in NAMES:
        source=DEST/(name+'.dat');back=DEST/(name+'_roundtrip.dat');binary=DEST/(name+'.plt')
        va,za=parse(source);vb,zb=parse(back);assert va==vb==VARS and za==zb
        assert binary.read_bytes()[:8]==b'#!TDV112' and (DEST/(name+'.lay')).exists()
        for z in za:
            a=np.array(z['rows']);assert a.shape==(400,len(VARS)) and np.all(np.diff(a[:,0])>0)
            for k in range(1,8):assert np.all(a[:,k]==a[0,k])
            for j in range(5):assert np.allclose(a[:,8+3*j],abs(a[:,9+3*j]+1j*a[:,10+3*j]),rtol=1e-14,atol=1e-14)
        checks.append(dict(file=binary.name,sha256=sha(binary),source_sha256=sha(source),roundtrip_sha256=sha(back),zones=len(za),points=sum(len(z['rows']) for z in za),double_values=sum(len(z['rows'])*len(va) for z in za),all_values_exact=True))
    js(DEST/'tecplot_validation.json',dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),tool='Tecplot 360 EX 2024 R1 / Preplot',binary_magic='#!TDV112',batch_log_sha256=sha(DEST/'batch.log'),files=checks))
    print('NATIVE_VERIFIED',[(r['file'],r['zones'],r['points'],r['double_values']) for r in checks],flush=True)

def main():
    if '--verify-native' in sys.argv:verify_native();return
    modes,continuous,evidence=collect();plot(modes,continuous,evidence);tecplot(modes);heatmap()
    assert snapshot()==evidence['protected_inputs_before']
    shutil.copy2(__file__,DEST/Path(__file__).name)
    print('PLOTS_AND_ASCII_READY',DEST,flush=True)

if __name__=='__main__':main()
