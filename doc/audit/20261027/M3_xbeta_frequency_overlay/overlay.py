"""Overlay one continuously identified M3 representative at five frequencies."""
import argparse,csv,datetime,hashlib,json,math,re,shutil
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

BASE=Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate/75°')
DEST=BASE/'outputs/M3_xbeta_frequency_overlay_20261027'
NONZERO=BASE/'outputs/M3_xbeta_four_omega_20261026'
ZERO=BASE/'outputs/M3_omega0_xbeta_20261026'
CONNECTION=BASE/'outputs/neutral_frequency_connection_20261026/f100/connection_report.json'
FS=[-100,-50,0,50,100]
CASES=['bluntplate','rough']
COLORS={-100:'#762a83',-50:'#2166ac',0:'#222222',50:'#1b9e77',100:'#d73027'}
TEC_COLORS={-100:'PURPLE',-50:'BLUE',0:'BLACK',50:'GREEN',100:'RED'}
VARS=['x_mm','beta','f_kHz','omega_nd','alpha_r','alpha_i','growth','x_flat_mm','station_coordinate','point_kind','case_id','right_residual','raw_MAC']
NAME='M3_xbeta_five_frequencies_same_branch'
NAMES=[NAME]+[case+'_M3_xbeta_five_frequencies' for case in CASES]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def load():
    datasets={};sources=[];checks=[];report=json.loads(CONNECTION.read_text())
    for case in CASES:
        for f in FS:
            directory=ZERO if f==0 else NONZERO
            component='negative_beta_zero' if f==0 else ('A' if f>0 else 'B')
            filename=f'{case}_M3_omega0_xbeta.csv' if f==0 else f'{case}_f{f:+d}kHz_{component}.csv'
            path=directory/filename;manifest=json.loads((directory/'delivery_manifest.json').read_text())
            previous=next(r for r in manifest['files'] if r['path']==filename);assert sha(path)==previous['sha256'],path
            with path.open(newline='') as stream:raw=list(csv.DictReader(stream))
            a=[]
            for r in raw:
                row={key:float(r[key]) for key in VARS if key not in ['f_kHz','omega_nd']}
                row.update(f_kHz=float(f),omega_nd=float(r['omega'] if f==0 else r['omega_nd']))
                if f!=0:assert float(r['f_kHz'])==f
                a.append([row[k] for k in VARS])
            a=np.array(a);assert np.isfinite(a).all() and np.max(abs(a[:,5]))<1e-8
            expected_omega=6.0727867021443425*f/50;assert np.all(a[:,3]==expected_omega)
            if f:
                assert np.array_equal(a[0,:4],a[-1,:4]) and a[-1,9]==2
            else:
                assert not np.array_equal(a[0,:2],a[-1,:2]) and abs(a[0,0]-1100)<1e-8 and abs(a[-1,0]-1100)<1e-8
            i=119 if case=='bluntplate' else 24
            seed=sorted([r for r in a if r[9]==0 and r[8]==i],key=lambda r:r[1]);assert len(seed)==2
            refs=[]
            for p in report['paths']:
                if p['case']!=case:continue
                refs.append(next(r for r in p['records'] if abs(r['f_kHz']-f)<1e-9))
            refs.sort(key=lambda r:r['beta']);assert len(refs)==2
            errors=[dict(beta_error=abs(r[1]-q['beta']),alpha_error=abs(complex(r[4],r[5])-complex(q['alpha_r'],q['alpha_i']))) for r,q in zip(seed,refs)]
            assert all(r['beta_error']<1e-7 and r['alpha_error']<1e-7 for r in errors)
            checks.append(dict(case=case,f_kHz=f,seed_station=i,comparison_to_continuous_frequency_path=errors))
            sources.append(dict(case=case,f_kHz=f,source_component=component,path=str(path),sha256=sha(path),points=len(a),closed=bool(f),x_range_mm=[float(a[:,0].min()),float(a[:,0].max())],beta_range=[float(a[:,1].min()),float(a[:,1].max())]))
            datasets[case,f]=a
    return datasets,dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),selection='One representative continued from stationary beta=-7 seed; no companion branch added',mapping='positive f: prior A; negative f: prior B; zero: original negative-beta branch',source_files=sources,continuity_source=dict(path=str(CONNECTION),sha256=sha(CONNECTION)),continuity_checks=checks)

def plot(datasets):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
    fig,axes=plt.subplots(2,2,figsize=(13.1,9.3));fig.subplots_adjust(left=.07,right=.985,bottom=.15,top=.89,wspace=.19,hspace=.33)
    for col,case in enumerate(CASES):
        for row in [0,1]:
            ax=axes[row,col]
            for f in FS:
                a=datasets[case,f];ax.plot(a[:,0],a[:,1],color=COLORS[f],lw=1.9 if f==0 else 1.6)
            ax.set_xlim(-35,1120 if row==0 else 160);ax.set_ylim(-25,10);ax.grid(alpha=.2);ax.set_axisbelow(True)
            ax.set_xlabel('Global x (mm)');ax.set_ylabel(r'$\beta$ (nondimensional)')
            ax.set_title(f'{case}: '+('full computed extent' if row==0 else 'upstream detail'))
    handles=[Line2D([0],[0],color=COLORS[f],lw=2,label=f'f = {f:+d} kHz' if f else 'f = 0') for f in FS]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.055),ncol=5,frameon=False)
    fig.suptitle('M3 x-beta neutral curves at five fixed physical frequencies',y=.97,fontsize=14)
    fig.text(.5,.02,'One representative branch continued from the stationary beta = -7 seed. The f = 0 boundary remains open at x = 1100 mm.',ha='center',fontsize=9)
    for ext in ['png','svg']:fig.savefig(DEST/(NAME+'_comparison.'+ext),dpi=220,facecolor='white')
    plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(13.1,5.3));fig.subplots_adjust(left=.07,right=.985,bottom=.24,top=.83,wspace=.19)
    for ax,case in zip(axes,CASES):
        for f in FS:
            a=datasets[case,f];ax.plot(a[:,0],a[:,1],color=COLORS[f],lw=1.9 if f==0 else 1.6)
        ax.set_xlim(-35,20);ax.set_ylim(-25,10);ax.grid(alpha=.2);ax.set_xlabel('Global x (mm)');ax.set_ylabel(r'$\beta$');ax.set_title(case)
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.04),ncol=5,frameon=False);fig.suptitle('Nose detail: the same representative at five frequencies',y=.97)
    for ext in ['png','svg']:fig.savefig(DEST/(NAME+'_nose.'+ext),dpi=220,facecolor='white')
    plt.close(fig)

def dat(name,items,datasets):
    lines=['TITLE="Single M3 representative at five fixed signed physical frequencies"','VARIABLES='+','.join('"'+k+'"' for k in VARS)]
    for case,f in items:
        a=datasets[case,f];lines += [f'ZONE T="{case}_f{f:+d}kHz", I={len(a)}, DATAPACKING=POINT','DT=('+' '.join(['DOUBLE']*len(VARS))+')']
        lines+=[' '.join('%.17g'%v for v in r) for r in a]
    (DEST/(name+'.dat')).write_text('\n'.join(lines)+'\n',encoding='ascii')

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
            row=[float(t) for t in s.split()];assert len(row)==len(variables) and all(math.isfinite(t) for t in row);zones[-1]['rows'].append(row)
    assert all(z['double'] for z in zones);return variables,zones

def macros():
    full=['#!MC 1410']
    for name in NAMES:
        _,zones=parse(DEST/(name+'.dat'));lines=['#!MC 1410',f'$!READDATASET \'"{name}.plt"\'','  READDATAOPTION = NEW','  RESETSTYLE = YES','  INITIALPLOTTYPE = XYLINE','$!PLOTTYPE = XYLINE','$!DELETELINEMAPS']
        for j,z in enumerate(zones,1):
            f=int(z['rows'][0][2]);pattern='DASHED' if len(zones)==10 and z['rows'][0][10]==2 else 'SOLID'
            lines+=['$!CREATELINEMAP',f'$!LINEMAP [{j}]',f"  NAME = '{z['name']}'",f'  ASSIGN{{ZONE={j} XAXISVAR=1 YAXISVAR=2}}',f'  LINES{{SHOW=YES COLOR={TEC_COLORS[f]} LINEPATTERN={pattern}}}','  SYMBOLS{SHOW=NO}']
        lines+=['$!ACTIVELINEMAPS = ['+','.join(str(j) for j in range(1,len(zones)+1))+']','$!VIEW FIT',f'$!SAVELAYOUT "{name}.lay"']
        (DEST/(name+'.mcr')).write_text('\n'.join(lines)+'\n',encoding='ascii')
        full += [f'$!INCLUDEMACRO "{name}.mcr"',f'$!WRITEDATASET "{name}_roundtrip.dat"','  BINARY = NO','  USEPOINTFORMAT = YES','  PRECISION = 17','  INCLUDETEXT = NO','  INCLUDEGEOM = NO','  INCLUDECUSTOMLABELS = NO']
    full+=['$!QUIT'];(DEST/'tecplot_roundtrip.mcr').write_text('\n'.join(full)+'\n',encoding='ascii')

def verify_native():
    assert 'ERROR' not in (DEST/'batch.log').read_text().upper();receipts=[]
    for name in NAMES:
        source=DEST/(name+'.dat');back=DEST/(name+'_roundtrip.dat');binary=DEST/(name+'.plt')
        va,za=parse(source);vb,zb=parse(back);assert va==vb and za==zb
        assert binary.read_bytes()[:8]==b'#!TDV112' and (DEST/(name+'.lay')).exists()
        zeros=0
        for z in za:
            a=np.array(z['rows']);f=a[0,2];assert f in FS and np.all(a[:,2]==f) and np.all(a[:,3]==6.0727867021443425*f/50)
            if f:assert np.array_equal(a[0,:4],a[-1,:4]) and a[-1,9]==2
            else:zeros+=1;assert not np.array_equal(a[0,:2],a[-1,:2])
        receipts.append(dict(file=binary.name,sha256=sha(binary),source_sha256=sha(source),roundtrip_sha256=sha(back),zones=len(za),points=sum(len(z['rows']) for z in za),double_values=sum(len(z['rows'])*len(va) for z in za),all_values_exact=True,zero_frequency_open_zones=zeros,nonzero_frequency_zones_closed=True))
    (DEST/'tecplot_validation.json').write_text(json.dumps(dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),tool='Tecplot 360 EX 2024 R1 / Preplot',binary_magic='#!TDV112',batch_log_sha256=sha(DEST/'batch.log'),files=receipts),indent=2)+'\n');print('NATIVE_VERIFIED',json.dumps(receipts),flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');p.add_argument('--render',action='store_true');args=p.parse_args()
    if args.verify:verify_native();return
    datasets,selection=load();DEST.mkdir(exist_ok=True)
    if not args.render:
        assert not (DEST/'selection_manifest.json').exists(),'Preserve completed selection'
        for (case,f),a in datasets.items():
            with (DEST/f'{case}_f{f:+d}kHz.csv').open('w',newline='') as stream:
                w=csv.writer(stream,lineterminator='\n');w.writerow(VARS);w.writerows(a)
        dat(NAME,list(datasets),datasets)
        for case in CASES:dat(case+'_M3_xbeta_five_frequencies',[(case,f) for f in FS],datasets)
        (DEST/'selection_manifest.json').write_text(json.dumps(selection,indent=2)+'\n');macros()
        shutil.copy2(NONZERO/'frequency_conversion.json',DEST/'frequency_conversion.json')
    plot(datasets);shutil.copy2(__file__,DEST/'overlay.py')
    print('OVERLAY',len(datasets),'curves',sum(len(a) for a in datasets.values()),'points',flush=True)
if __name__=='__main__':main()
