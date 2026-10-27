"""Publish original and extended signed-frequency neutral curves, without mirrors."""
import csv,datetime,hashlib,json,re,shutil,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
WIN=Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate/75°')
CALC=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/neutral_xbeta_extended_frequencies_20261027')
PRIOR=WIN/'outputs/M3_xbeta_nine_frequencies_20261027'
DEST=WIN/'outputs/M3_xbeta_extended_frequencies_20261027'
NEW=[-112,-110,125,150,175,200]
FS=[-112,-110,-100,-75,-50,-25,0,25,50,75,100,125,150,175,200]
CASES=['bluntplate','rough'];LABELS=['Baseline blunt plate','Optimized nose: continuous curvature']
COLORS={-112:'#00a6a6',-110:'#4c63d2',-100:'#762a83',-75:'#5954a6',-50:'#2166ac',-25:'#149fb3',0:'#161616',25:'#658e28',50:'#1b9e77',75:'#e48120',100:'#d73027',125:'#8c564b',150:'#e377c2',175:'#bc8f00',200:'#7f3c8d'}
STYLES={-112:'--',-110:'-',125:'--',150:'-',175:'--',200:'-'}
VARS=['x_mm','beta','f_kHz','omega_nd','alpha_r','alpha_i','growth','x_flat_mm','station_coordinate','point_kind','case_id','right_residual','raw_MAC']
NAME='M3_xbeta_extended_frequencies'
NAMES=[NAME]+[case+'_M3_xbeta_extended_frequencies' for case in CASES]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def js(p,d):p.write_text(json.dumps(d,indent=2,ensure_ascii=False)+'\n')
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
            row=[float(x) for x in s.split()];assert len(row)==len(variables) and np.isfinite(row).all();zones[-1]['rows'].append(row)
    assert all(z['double'] for z in zones);return variables,zones
def load():
    previous=json.loads((PRIOR/'delivery_manifest.json').read_text());audit=json.loads((CALC/'outputs/new_frequency_audit.json').read_text());data={};sources=[]
    for case in CASES:
        for f in FS:
            if f not in NEW:
                path=PRIOR/f'{case}_f{f:+d}kHz.csv';record=next(r for r in previous['files'] if r['path']==path.name);assert sha(path)==record['sha256']
                with path.open() as stream:a=np.array([[float(r[k]) for k in VARS] for r in csv.DictReader(stream)])
                status='closed_verified' if f else 'open_at_domain_boundaries';shutil.copy2(path,DEST/path.name)
            else:
                path=CALC/'outputs'/f'{case}_f{f}kHz_S_final.json';r=next(r for r in audit['results'] if r['case']==case and r['f_kHz']==f);assert sha(path)==r['final_sha256'];d=json.loads(path.read_text());status=d['status']
                a=np.array([[p['x_mm'],p['beta_r'],f,p['omega_r'],p['alpha_r'],p['alpha_i'],p['growth'],p['x_flat_mm'],p['xi'],p['kind'],CASES.index(case)+1,p['right_residual'],p['raw_MAC']] for p in d['points']])
                np.savetxt(DEST/f'{case}_f{f:+d}kHz.csv',a,delimiter=',',header=','.join(VARS),comments='',fmt='%.17g')
                shutil.copy2(path,DEST/(path.stem+'_calculation_evidence.json'))
            assert np.isfinite(a).all() and max(abs(a[:,5]))<1e-8 and np.all(a[:,3]==6.0727867021443425*f/50)
            if status=='closed_verified':assert np.array_equal(a[0,:4],a[-1,:4]) and a[-1,9]==2
            else:assert abs(a[0,0]-1100)<1e-8 and abs(a[-1,0]-1100)<1e-8
            data[case,f]=a;sources.append(dict(case=case,f_kHz=f,path=str(path),sha256=sha(path),status=status,points=len(a),x_range_mm=[float(a[:,0].min()),float(a[:,0].max())],beta_range=[float(a[:,1].min()),float(a[:,1].max())]))
    return data,dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),source_files=sources,frequencies_kHz=FS,new_frequencies_kHz=NEW,model='Same frozen local tangent-plane LST, zero explicit curvature',selection='One representative connected to the stationary beta=-7 seed; no added conjugate representative',negative_frequency_limit='Continuation at selected native sections stopped near -112 to -114 kHz. No exhaustive no-root claim at -125/-150/-200 kHz.')
def draw(ax,data,case,frequencies,xlim,ylim,title):
    for f in frequencies:
        a=data[case,f];ax.plot(a[:,0],a[:,1],color=COLORS[f],ls=STYLES.get(f,'-'),lw=1.7 if f==0 else 1.5,label=f'{f:+d}' if f else '0',zorder=5 if f==0 else 2)
    ax.set(xlim=xlim,ylim=ylim,xlabel='Global x (mm)',ylabel=r'$\beta$ (nondimensional)',title=title);ax.grid(alpha=.2)
def save(fig,suffix):
    for ext in ['png','svg']:fig.savefig(DEST/(NAME+'_'+suffix+'.'+ext),dpi=220,facecolor='white')
    plt.close(fig)
def plot(data):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
    fig,axes=plt.subplots(2,2,figsize=(13.5,9.3));fig.subplots_adjust(left=.065,right=.985,bottom=.09,top=.91,wspace=.19,hspace=.47)
    for col,case in enumerate(CASES):
        draw(axes[0,col],data,case,[100,125,150,175,200],(-35,85),(-6,21),LABELS[col]+' / positive f')
        draw(axes[1,col],data,case,[-100,-110,-112],(-35,-15),(-26,-16),LABELS[col]+' / negative f')
        axes[0,col].legend(title='f (kHz)',ncol=3,fontsize=9,loc='upper right',framealpha=.9)
        axes[1,col].legend(title='f (kHz)',ncol=3,fontsize=9,loc='lower right',framealpha=.9)
    fig.suptitle('M3 neutral curves: extended physical-frequency range',fontsize=15)
    fig.text(.5,.018,'Same continuously tracked representative; curves include corrected physical-profile turning points. Zero explicit curvature in LST.',ha='center',fontsize=9)
    save(fig,'new_range_comparison')
    fig,axes=plt.subplots(2,2,figsize=(14,10.2));fig.subplots_adjust(left=.065,right=.985,bottom=.21,top=.9,wspace=.19,hspace=.32)
    for col,case in enumerate(CASES):
        draw(axes[0,col],data,case,FS,(-35,1120),(-26,21),LABELS[col]+' / full extent')
        draw(axes[1,col],data,case,FS,(-35,250),(-26,21),LABELS[col]+' / upstream')
    handles=[Line2D([0],[0],color=COLORS[f],ls=STYLES.get(f,'-'),lw=1.7,label=f'{f:+d} kHz' if f else '0 kHz') for f in FS]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.055),ncol=5,frameon=False)
    fig.suptitle('M3 x-beta neutral curves at 15 fixed physical frequencies',fontsize=15)
    fig.text(.5,.018,'One representative; f = 0 remains open at x = 1100 mm. Same zero-explicit-curvature model as the prior nine-frequency results.',ha='center',fontsize=9)
    save(fig,'all_frequencies')
    fig,axes=plt.subplots(1,2,figsize=(13.5,6.5));fig.subplots_adjust(left=.065,right=.985,bottom=.24,top=.85,wspace=.19)
    pos=[f for f in FS if f>0]
    for col,case in enumerate(CASES):draw(axes[col],data,case,pos,(-35,250),(-16,21),LABELS[col])
    fig.legend(handles=[Line2D([0],[0],color=COLORS[f],ls=STYLES.get(f,'-'),lw=1.7,label=f'{f:+d} kHz') for f in pos],loc='lower center',bbox_to_anchor=(.5,.035),ncol=4,frameon=False)
    fig.suptitle('M3 positive physical frequencies: 25 to 200 kHz',fontsize=15)
    save(fig,'positive_frequency_detail')
def dat(name,items,data):
    with (DEST/(name+'.dat')).open('w') as stream:
        stream.write('TITLE="Extended signed physical frequencies; one M3 representative"\nVARIABLES='+','.join('"'+s+'"' for s in VARS)+'\n')
        for case,f in items:
            a=data[case,f];stream.write(f'ZONE T="{case}_f{f:+d}kHz", I={len(a)}, DATAPACKING=POINT\nDT=('+','.join(['DOUBLE']*len(VARS))+')\n');np.savetxt(stream,a,fmt='%.17g')
def macros():
    out=['#!MC 1410'];colors=['CYAN','BLUE','PURPLE','PURPLE','BLUE','CYAN','BLACK','GREEN','GREEN','RED','RED','BLACK','PURPLE','GREEN','BLUE'];patterns=['DASHED','DOTTED','SOLID','DASHED','SOLID','DASHED','SOLID','DOTTED','SOLID','DASHED','SOLID','DASHED','DOTTED','DASHED','DASHED']
    for name in NAMES:
        _,zones=parse(DEST/(name+'.dat'));out += [f'$!READDATASET \'"{name}.plt"\'',' READDATAOPTION=NEW',' RESETSTYLE=YES',' INITIALPLOTTYPE=XYLINE','$!PLOTTYPE=XYLINE','$!DELETELINEMAPS']
        for j,z in enumerate(zones,1):
            k=FS.index(int(z['rows'][0][2]));out += ['$!CREATELINEMAP',f'$!LINEMAP [{j}]',f" NAME='{z['name']}'",f' ASSIGN{{ZONE={j} XAXISVAR=1 YAXISVAR=2}}',f' LINES{{SHOW=YES COLOR={colors[k]} LINEPATTERN={patterns[k]}}}',' SYMBOLS{SHOW=NO}']
        out += ['$!ACTIVELINEMAPS=['+','.join(map(str,range(1,len(zones)+1)))+']','$!VIEW FIT',f'$!SAVELAYOUT "{name}.lay"',f'$!WRITEDATASET "{name}_roundtrip.dat"',' BINARY=NO',' USEPOINTFORMAT=YES',' PRECISION=17',' INCLUDETEXT=NO',' INCLUDEGEOM=NO',' INCLUDECUSTOMLABELS=NO']
    out+=['$!QUIT'];(DEST/'tecplot_roundtrip.mcr').write_text('\n'.join(out)+'\n',encoding='ascii')
def verify():
    assert 'ERROR' not in (DEST/'batch.log').read_text().upper();files=[]
    for name in NAMES:
        a=DEST/(name+'.dat');b=DEST/(name+'_roundtrip.dat');binary=DEST/(name+'.plt');va,za=parse(a);vb,zb=parse(b)
        assert va==vb==VARS and za==zb and binary.read_bytes()[:8]==b'#!TDV112' and (DEST/(name+'.lay')).exists()
        files.append(dict(file=binary.name,sha256=sha(binary),source_sha256=sha(a),roundtrip_sha256=sha(b),zones=len(za),points=sum(len(z['rows']) for z in za),all_values_exact=True))
    js(DEST/'tecplot_validation.json',dict(tool='Tecplot 360 EX 2024 R1 / Preplot',binary_magic='#!TDV112',files=files));print('NATIVE_VERIFIED',files,flush=True)
def main():
    if '--verify' in sys.argv:verify();return
    DEST.mkdir(parents=True,exist_ok=True);data,manifest=load();js(DEST/'selection_manifest.json',manifest)
    dat(NAME,list(data),data)
    for case in CASES:dat(case+'_M3_xbeta_extended_frequencies',[(case,f) for f in FS],data)
    plot(data);macros();shutil.copy2(__file__,DEST/'publish.py')
    for name in ['new_frequency_audit.json','bootstrap_receipt.json','neutral_seeds.json','seed_extension_summary.json','frequency_conversion.json','run_summary.json','input_manifest.json']:
        shutil.copy2(CALC/'outputs'/name,DEST/name)
    print('PUBLISHED',len(data),sum(len(a) for a in data.values()),flush=True)
if __name__=='__main__':main()
