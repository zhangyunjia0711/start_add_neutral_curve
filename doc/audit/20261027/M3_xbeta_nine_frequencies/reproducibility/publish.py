"""Publish nine frequencies of one tracked M3 family; no added companion."""
import argparse,csv,datetime,hashlib,json,re,shutil,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
BASE=Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate/75°')
CALC=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/neutral_xbeta_extra_frequencies_20261027')
PRIOR=BASE/'outputs/M3_xbeta_frequency_overlay_20261027'
DEST=BASE/'outputs/M3_xbeta_nine_frequencies_20261027'
FS=[-100,-75,-50,-25,0,25,50,75,100]
NEW=[-75,-25,25,75]
CASES=['bluntplate','rough']
COLORS={-100:'#762a83',-75:'#5954a6',-50:'#2166ac',-25:'#149fb3',0:'#161616',25:'#658e28',50:'#1b9e77',75:'#e48120',100:'#d73027'}
TEC={-100:'PURPLE',-75:'CUSTOM1',-50:'BLUE',-25:'CYAN',0:'BLACK',25:'CUSTOM2',50:'GREEN',75:'ORANGE',100:'RED'}
VARS=['x_mm','beta','f_kHz','omega_nd','alpha_r','alpha_i','growth','x_flat_mm','station_coordinate','point_kind','case_id','right_residual','raw_MAC']
NAME='M3_xbeta_nine_frequencies_same_branch'
NAMES=[NAME]+[case+'_M3_xbeta_nine_frequencies' for case in CASES]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def js(p,r):p.write_text(json.dumps(r,indent=2,ensure_ascii=False)+'\n')
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

def load():
    data={};sources=[];previous=json.loads((PRIOR/'delivery_manifest.json').read_text());seeds=json.loads((CALC/'outputs/neutral_seeds.json').read_text())['records']
    for case in CASES:
        for f in FS:
            if f not in NEW:
                path=PRIOR/f'{case}_f{f:+d}kHz.csv'
                record=next(r for r in previous['files'] if r['path']==path.name);assert sha(path)==record['sha256']
                with path.open() as stream:rows=[[float(r[k]) for k in VARS] for r in csv.DictReader(stream)]
                status='closed_verified' if f else 'open_at_domain_boundaries'
            else:
                path=CALC/'outputs'/f'{case}_f{f}kHz_S_final.json';d=json.loads(path.read_text());status=d['status']
                assert status in ['closed_verified','open_at_domain_boundaries']
                assert d['original_checks'] and d['physical_tips']
                rows=[[p['x_mm'],p['beta_r'],f,p['omega_r'],p['alpha_r'],p['alpha_i'],p['growth'],p['x_flat_mm'],p['xi'],p['kind'],CASES.index(case)+1,p['right_residual'],p['raw_MAC']] for p in d['points']]
                i=119 if case=='bluntplate' else 24
                native=sorted([r for r in rows if r[9]==0 and r[8]==i],key=lambda r:r[1]);assert len(native)==2
                refs=sorted([r for r in seeds if r['case']==case and r['f_kHz']==f],key=lambda r:r['beta_r']);assert len(refs)==2
                assert all(abs(a[1]-b['beta_r'])<1e-7 and abs(complex(a[4],a[5])-complex(b['alpha_r'],b['alpha_i']))<1e-7 for a,b in zip(native,refs))
                shutil.copy2(path,DEST/(path.stem+'_calculation_evidence.json'))
            a=np.array(rows);assert np.isfinite(a).all() and max(abs(a[:,5]))<1e-8 and np.all(a[:,3]==6.0727867021443425*f/50)
            if status=='closed_verified':assert np.array_equal(a[0,:4],a[-1,:4]) and a[-1,9]==2
            else:assert abs(a[0,0]-1100)<1e-8 and abs(a[-1,0]-1100)<1e-8 and not np.array_equal(a[0,:2],a[-1,:2])
            data[case,f]=a
            sources.append(dict(case=case,f_kHz=f,path=str(path),sha256=sha(path),status=status,points=len(a),x_range_mm=[float(a[:,0].min()),float(a[:,0].max())],beta_range=[float(a[:,1].min()),float(a[:,1].max())]))
    return data,dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),selection='One representative continued from the stationary beta=-7 seed; direct signed-frequency calculations for the four added frequencies',source_files=sources,new_frequencies_kHz=NEW,source_seeds=dict(path=str(CALC/'outputs/neutral_seeds.json'),sha256=sha(CALC/'outputs/neutral_seeds.json')))

def plot(data):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
    handles=[Line2D([0],[0],color=COLORS[f],lw=2,label=f'f = {f:+d} kHz' if f else 'f = 0') for f in FS]
    for view in ['comparison','upstream','nose']:
        if view=='comparison':
            fig,axes=plt.subplots(2,2,figsize=(13.5,10));fig.subplots_adjust(left=.07,right=.985,bottom=.18,top=.89,wspace=.19,hspace=.33)
            for col,case in enumerate(CASES):
                for row in [0,1]:
                    ax=axes[row,col]
                    for f in FS:
                        a=data[case,f];ax.plot(a[:,0],a[:,1],color=COLORS[f],lw=1.8 if f==0 else 1.4,zorder=10 if f==0 else 2)
                    ax.set(xlim=(-35,1120 if row==0 else 160),ylim=(-26,10),xlabel='Global x (mm)',ylabel=r'$\beta$ (nondimensional)',title=case+(': full computed extent' if row==0 else ': upstream detail'));ax.grid(alpha=.2)
        else:
            fig,axes=plt.subplots(1,2,figsize=(13.5,6.2));fig.subplots_adjust(left=.07,right=.985,bottom=.25,top=.83,wspace=.19)
            for ax,case in zip(axes,CASES):
                for f in FS:
                    a=data[case,f];ax.plot(a[:,0],a[:,1],color=COLORS[f],lw=1.8 if f==0 else 1.5,zorder=10 if f==0 else 2)
                ax.set(xlim=(-35,160 if view=='upstream' else 5),ylim=(-26,10),xlabel='Global x (mm)',ylabel=r'$\beta$ (nondimensional)',title=case+': '+view+' detail');ax.grid(alpha=.2)
        fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.055),ncol=5,frameon=False)
        fig.suptitle('M3 x-beta neutral curves at nine fixed physical frequencies',y=.97,fontsize=14)
        fig.text(.5,.018,'One continuously tracked representative. Open boundaries at x = 1100 mm are retained; no artificial closure.',ha='center',fontsize=9)
        for ext in ['png','svg']:fig.savefig(DEST/(NAME+'_'+view+'.'+ext),dpi=220,facecolor='white')
        plt.close(fig)

def dat(name,items,data):
    with (DEST/(name+'.dat')).open('w') as stream:
        stream.write('TITLE="Nine fixed physical frequencies; one M3 representative"\nVARIABLES='+','.join('"'+s+'"' for s in VARS)+'\n')
        for case,f in items:
            a=data[case,f];stream.write(f'ZONE T="{case}_f{f:+d}kHz", I={len(a)}, DATAPACKING=POINT\nDT=('+','.join(['DOUBLE']*len(VARS))+')\n');np.savetxt(stream,a,fmt='%.17g')

def macros():
    full=['#!MC 1410']
    for name in NAMES:
        _,zones=parse(DEST/(name+'.dat'))
        full += [f'$!READDATASET \'"{name}.plt"\'',' READDATAOPTION=NEW',' RESETSTYLE=YES',' INITIALPLOTTYPE=XYLINE','$!PLOTTYPE=XYLINE','$!DELETELINEMAPS']
        # Use built-in colors only, avoiding installation-specific custom color tables.
        for j,z in enumerate(zones,1):
            f=int(z['rows'][0][2]);color={-75:'PURPLE',25:'GREEN',75:'RED'}.get(f,TEC[f]);pattern='DOTTED' if f in [-75,25] else ('DASHED' if f==75 else 'SOLID')
            full += ['$!CREATELINEMAP',f'$!LINEMAP [{j}]',f" NAME='{z['name']}'",f' ASSIGN{{ZONE={j} XAXISVAR=1 YAXISVAR=2}}',f' LINES{{SHOW=YES COLOR={color} LINEPATTERN={pattern}}}',' SYMBOLS{SHOW=NO}']
        full += ['$!ACTIVELINEMAPS=['+','.join(map(str,range(1,len(zones)+1)))+']','$!VIEW FIT',f'$!SAVELAYOUT "{name}.lay"',f'$!WRITEDATASET "{name}_roundtrip.dat"',' BINARY=NO',' USEPOINTFORMAT=YES',' PRECISION=17',' INCLUDETEXT=NO',' INCLUDEGEOM=NO',' INCLUDECUSTOMLABELS=NO']
    full+=['$!QUIT'];(DEST/'tecplot_roundtrip.mcr').write_text('\n'.join(full)+'\n',encoding='ascii')

def verify():
    assert 'ERROR' not in (DEST/'batch.log').read_text().upper();receipts=[]
    selection=json.loads((DEST/'selection_manifest.json').read_text())
    status={(r['case'],r['f_kHz']):r['status'] for r in selection['source_files']}
    for name in NAMES:
        source=DEST/(name+'.dat');back=DEST/(name+'_roundtrip.dat');binary=DEST/(name+'.plt')
        va,za=parse(source);vb,zb=parse(back);assert va==vb==VARS and za==zb
        assert binary.read_bytes()[:8]==b'#!TDV112' and (DEST/(name+'.lay')).exists()
        opened=0;closed=0
        for z in za:
            a=np.array(z['rows']);f=int(a[0,2]);case=CASES[int(a[0,10])-1]
            assert np.all(a[:,2]==f) and np.all(a[:,3]==6.0727867021443425*f/50)
            if status[case,f]=='closed_verified':closed+=1;assert np.array_equal(a[0,:4],a[-1,:4]) and a[-1,9]==2
            else:opened+=1;assert abs(a[0,0]-1100)<1e-8 and abs(a[-1,0]-1100)<1e-8 and not np.array_equal(a[0,:2],a[-1,:2])
        receipts.append(dict(file=binary.name,sha256=sha(binary),source_sha256=sha(source),roundtrip_sha256=sha(back),zones=len(za),points=sum(len(z['rows']) for z in za),double_values=sum(len(z['rows'])*len(va) for z in za),all_values_exact=True,open_zones=opened,closed_zones=closed))
    js(DEST/'tecplot_validation.json',dict(tool='Tecplot 360 EX 2024 R1 / Preplot',binary_magic='#!TDV112',batch_log_sha256=sha(DEST/'batch.log'),files=receipts))
    print('NATIVE_VERIFIED',[(r['file'],r['zones'],r['points']) for r in receipts],flush=True)

def main():
    if '--verify' in sys.argv:verify();return
    DEST.mkdir(exist_ok=True);data,manifest=load()
    for (case,f),a in data.items():np.savetxt(DEST/f'{case}_f{f:+d}kHz.csv',a,delimiter=',',header=','.join(VARS),comments='',fmt='%.17g')
    js(DEST/'selection_manifest.json',manifest)
    dat(NAME,list(data),data)
    for case in CASES:dat(case+'_M3_xbeta_nine_frequencies',[(case,f) for f in FS],data)
    plot(data);macros();shutil.copy2(__file__,DEST/'publish.py')
    conversion=json.loads((CALC/'outputs/frequency_conversion.json').read_text());conversion['frequencies']=[dict(f_kHz=f,omega_nd=6.0727867021443425*f/50) for f in FS];js(DEST/'frequency_conversion.json',conversion)
    for name in ['bootstrap_receipt.json','neutral_seeds.json','run_summary.json','input_manifest.json','new_frequency_audit.json']:
        source=CALC/'outputs'/name
        if source.exists():shutil.copy2(source,DEST/name)
    print('PUBLISHED',len(data),'curves',sum(len(a) for a in data.values()),'points',flush=True)
if __name__=='__main__':main()
