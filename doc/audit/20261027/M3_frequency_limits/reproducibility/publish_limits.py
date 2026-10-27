"""Deliver signed M3 frequency limits, connected pockets, and validation figures."""
import csv,datetime,hashlib,json,re,shutil,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
WIN=Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate/75°')
CALC=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/neutral_m3_frequency_limits_20261027')
DEST=WIN/'outputs/M3_frequency_limits_20261027'
PRIOR=WIN/'outputs/M3_xbeta_extended_frequencies_20261027'
CASES=['bluntplate','rough'];LABELS=['Baseline blunt plate','Optimized continuous-curvature nose']
VARS=['x_mm','beta','f_kHz','omega_nd','alpha_r','alpha_i','growth','x_flat_mm','station_coordinate','point_kind','case_id','right_residual','raw_MAC']
W=6.0727867021443425/50
NAME='M3_frequency_limits_neutral_curves';NAMES=[NAME]+[c+'_M3_frequency_limits' for c in CASES]
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
        elif s.startswith('DT=') or s.startswith('DT ='):assert len(re.findall('DOUBLE',s))==len(variables);zones[-1]['double']=True
        elif zones and re.match(r'^[+-]?(?:[0-9]|\.[0-9])',s):
            row=[float(v) for v in s.split()];assert len(row)==len(variables) and np.isfinite(row).all();zones[-1]['rows'].append(row)
    assert all(z['double'] for z in zones);return variables,zones
def rows_from_json(case,f,d):
    return np.array([[p['x_mm'],p['beta_r'],f,p['omega_r'],p['alpha_r'],p['alpha_i'],p['growth'],p['x_flat_mm'],p['xi'],p['kind'],CASES.index(case)+1,p['right_residual'],p['raw_MAC']] for p in d['points']])
def load():
    curves=[];sources=[];previous=json.loads((PRIOR/'delivery_manifest.json').read_text())
    for case in CASES:
        for f,region in [(200.,'positive'),(-112.,'nose')]:
            file=PRIOR/f'{case}_f{f:+g}kHz.csv';r=next(r for r in previous['files'] if r['path']==file.name);assert sha(file)==r['sha256']
            with file.open() as stream:a=np.array([[float(row[k]) for k in VARS] for row in csv.DictReader(stream)])
            curves.append(dict(case=case,f_kHz=f,region=region,rows=a,source=str(file)));sources.append(dict(case=case,f_kHz=f,region=region,path=str(file),sha256=sha(file),prior=True))
        positive=[250.,275.,290.,291.6] if case=='bluntplate' else [230.,250.,264.53]
        for f in positive:
            file=CALC/'outputs'/f'{case}_f{f:g}kHz_neutral.json';d=json.loads(file.read_text());assert d['status']=='closed_verified';a=rows_from_json(case,f,d);curves.append(dict(case=case,f_kHz=f,region='positive',rows=a,source=str(file)));sources.append(dict(case=case,f_kHz=f,region='positive',path=str(file),sha256=sha(file),prior=False))
        for f in [-112.,-150.,-200.]:
            file=CALC/'outputs'/f'{case}_f{f}kHz_S_final.json';d=json.loads(file.read_text());assert d['status']=='closed_verified';a=rows_from_json(case,f,d);curves.append(dict(case=case,f_kHz=f,region='downstream',rows=a,source=str(file)));sources.append(dict(case=case,f_kHz=f,region='downstream',path=str(file),sha256=sha(file),prior=False))
        f=-240.2 if case=='bluntplate' else -237.14;file=CALC/'outputs'/f'{case}_secondary_f{f:g}kHz_neutral.json';d=json.loads(file.read_text());assert d['status']=='closed_verified';curves.append(dict(case=case,f_kHz=f,region='downstream',rows=rows_from_json(case,f,d),source=str(file)));sources.append(dict(case=case,f_kHz=f,region='downstream',path=str(file),sha256=sha(file),prior=False))
    limits=[]
    for case in CASES:
        for sign in [1,-1]:
            label=f'{case}_{"positive" if sign>0 else "negative_secondary"}';file=CALC/'outputs'/f'{label}_refined_limit.json';d=json.loads(file.read_text());r=d['critical'];limits.append(dict(case=case,sign=sign,f_kHz=r['f_kHz'],omega_nd=r['omega_r'],x_mm=r['x_mm'],xi=r['xi'],beta=r['beta_r'],alpha_r=r['alpha_r'],right_residual=r['right_residual'],source=str(file),sha256=sha(file),record=d))
    for c in curves:
        a=c['rows'];assert np.isfinite(a).all() and np.max(abs(a[:,5]))<1e-8 and np.allclose(a[:,3],W*c['f_kHz'],rtol=1e-13,atol=1e-13);assert np.array_equal(a[0,:4],a[-1,:4]) and a[-1,9]==2
    return curves,limits,sources
def color(c):
    if c['region']=='nose':return '#762a83'
    f=c['f_kHz']
    if f>0:return {200.:'#762a83',230.:'#1b9e77',250.:'#2166ac',275.:'#e6ab02',290.:'#d95f02',291.6:'#c51b7d',264.53:'#c51b7d'}[f]
    return {-112.:'#762a83',-150.:'#2166ac',-200.:'#1b9e77',-240.2:'#c51b7d',-237.14:'#c51b7d'}[f]
def save(fig,name):
    for ext in ['png','svg']:fig.savefig(DEST/(name+'.'+ext),dpi=220,facecolor='white')
    plt.close(fig)
def plot(curves,limits):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','axes.formatter.useoffset':False})
    neg=[c for c in curves if c['f_kHz']<0];xmax=max(55.,np.ceil(max(c['rows'][:,0].max() for c in neg)/10)*10+5)
    fig,axes=plt.subplots(2,2,figsize=(14,10.2));fig.subplots_adjust(left=.07,right=.985,top=.91,bottom=.09,hspace=.36,wspace=.19)
    for col,case in enumerate(CASES):
        for row,sign in [(0,1),(1,-1)]:
            ax=axes[row,col];seen=set()
            for c in [c for c in curves if c['case']==case and np.sign(c['f_kHz'])==sign]:
                f=c['f_kHz'];a=c['rows'];label=f'{f:+g} kHz' if f not in seen else None;seen.add(f);ax.plot(a[:,0],a[:,1],color=color(c),lw=1.5,label=label)
            r=next(r for r in limits if r['case']==case and r['sign']==sign);ax.plot(r['x_mm'],r['beta'],'k*',ms=10,label=f'limit {r["f_kHz"]:+.3f} kHz');ax.set(xlim=(-30,35) if sign>0 else (-35,xmax),ylim=(6,24) if sign>0 else (-38,-10),xlabel='Global x (mm)',ylabel=r'$\beta$ (nondimensional)',title=LABELS[col]+(' / positive f' if sign>0 else ' / negative f'));ax.grid(alpha=.2);ax.legend(fontsize=8,loc='upper right' if sign>0 else 'upper left',framealpha=.92)
    fig.suptitle('M3 neutral curves continued to frequency limits',fontsize=15)
    fig.text(.5,.018,'One tracked representative; no added conjugate copy. Negative f includes the newly resolved downstream pocket. Zero explicit curvature in LST.',ha='center',fontsize=9);save(fig,'M3_frequency_limits_comparison')
    fig,axes=plt.subplots(2,2,figsize=(13.5,9.5));fig.subplots_adjust(left=.08,right=.985,top=.91,bottom=.09,hspace=.36,wspace=.22)
    for col,case in enumerate(CASES):
        for row,sign in [(0,1),(1,-1)]:
            r=next(r for r in limits if r['case']==case and r['sign']==sign);cs=[c for c in curves if c['case']==case and np.sign(c['f_kHz'])==sign];c=(max if sign>0 else min)(cs,key=lambda c:c['f_kHz']);a=c['rows'];ax=axes[row,col];ax.plot(a[:,0]-r['x_mm'],a[:,1]-r['beta'],color='#2166ac',lw=1.6,label=f'last loop: {c["f_kHz"]:+g} kHz');ax.plot(0,0,'r*',ms=11,label=f'limit: {r["f_kHz"]:+.6f} kHz');ax.set(xlabel=r'$x-x_c$ (mm)',ylabel=r'$\beta-\beta_c$',title=LABELS[col]+(' / positive end' if sign>0 else ' / negative end'));ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Last resolved closed loops and extinction points',fontsize=15);fig.text(.5,.018,'Additional points solve reassembled interpolated local profiles; these subgrid loops are not CFD mesh refinement or a mesh-convergence study.',ha='center',fontsize=9);save(fig,'M3_frequency_limits_last_loops')
    fig,axes=plt.subplots(2,2,figsize=(12.5,8.5));fig.subplots_adjust(left=.09,right=.985,top=.91,bottom=.1,hspace=.34,wspace=.25)
    for col,case in enumerate(CASES):
        for row,sign in [(0,1),(1,-1)]:
            r=next(r for r in limits if r['case']==case and r['sign']==sign);tests=r['record']['frequency_side_checks'];xy=sorted([(sign*(t['f_kHz']-r['f_kHz']),t['maximum_growth']) for t in tests]+[(0.,0.)]);a=np.array(xy);ax=axes[row,col];ax.plot(a[:,0],a[:,1],'o-',color='#2166ac',ms=4);ax.axhline(0,color='k',lw=.8);ax.axvline(0,color='k',lw=.8);ax.axvspan(-1,0,color='#1b9e77',alpha=.08);ax.axvspan(0,1,color='#2166ac',alpha=.06);ax.set(xlabel='Outward frequency offset (kHz)',ylabel=r'Local $\max_\beta[-\mathrm{Im}(\alpha)]$',title=LABELS[col]+f' / {r["f_kHz"]:+.3f} kHz');ax.grid(alpha=.2)
    fig.suptitle('Growth-sign verification at the critical positions',fontsize=15);fig.text(.5,.02,'Local beta maxima at each critical position: negative offsets are inside the window; positive offsets are outside. Not an N-factor.',ha='center',fontsize=9);save(fig,'M3_frequency_limits_growth_checks')
    fig,axes=plt.subplots(1,2,figsize=(14,5.8));fig.subplots_adjust(left=.07,right=.985,top=.86,bottom=.22,wspace=.17)
    kinds=[('positive_native','#d73027','positive merger'),('negative_native','#762a83','nose merger'),('negative_secondary_native','#2166ac','downstream merger'),('negative_tail','#1b9e77','tail merger')]
    for ax,case,title in zip(axes,CASES,LABELS):
        for kind,col,label in kinds:
            d=json.loads((CALC/'outputs'/f'{case}_{kind}_fold_curve.json').read_text());rs=d['records'];ax.plot([r['x_mm'] for r in rs],[r['f_kHz'] for r in rs],color=col,lw=1.3,label=label)
        for kind in ['secondary','tail']:
            d=json.loads((CALC/'outputs'/f'{case}_negative_{kind}_arc.json').read_text());rs=d['records'];ax.plot([r['point']['x_mm'] for r in rs],[r['point']['f_kHz'] for r in rs],color='#777777',ls='--',lw=1,label='operator guide' if kind=='secondary' else None)
        for r in limits:
            if r['case']==case:ax.plot(r['x_mm'],r['f_kHz'],'k*',ms=9)
        ax.set(xlim=(-35,1100),ylim=(-260,315),xlabel='Global x (mm)',ylabel='Merger frequency (kHz)',title=title);ax.grid(alpha=.2)
    handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=5,bbox_to_anchor=(.5,.025),frameon=False);fig.suptitle('Neutral-merger branches used to locate the frequency extrema',fontsize=14);save(fig,'M3_frequency_limits_merger_map')
def dat(name,items):
    with (DEST/(name+'.dat')).open('w') as stream:
        stream.write('TITLE="Tracked M3 frequency limits and connected pockets"\nVARIABLES='+','.join('"'+v+'"' for v in VARS)+'\n')
        for title,a in items:stream.write(f'ZONE T="{title}", I={len(a)}, DATAPACKING=POINT\nDT=('+','.join(['DOUBLE']*len(VARS))+')\n');np.savetxt(stream,a,fmt='%.17g')
def macros():
    out=['#!MC 1410']
    colors=['PURPLE','BLUE','GREEN','RED','CYAN','BLACK'];patterns=['SOLID','DASHED','DOTTED']
    for name in NAMES:
        _,zones=parse(DEST/(name+'.dat'));out += [f'$!READDATASET \'"{name}.plt"\'',' READDATAOPTION=NEW',' RESETSTYLE=YES',' INITIALPLOTTYPE=XYLINE','$!PLOTTYPE=XYLINE','$!DELETELINEMAPS']
        for j,z in enumerate(zones,1):
            critical=z['rows'][0][9]==9;out+=['$!CREATELINEMAP',f'$!LINEMAP [{j}]',f" NAME='{z['name']}'",f' ASSIGN{{ZONE={j} XAXISVAR=1 YAXISVAR=2}}',f' LINES{{SHOW={"NO" if critical else "YES"} COLOR={colors[(j-1)%6]} LINEPATTERN={patterns[((j-1)//6)%3]}}}',f' SYMBOLS{{SHOW={"YES" if critical else "NO"}}}']
        out += ['$!ACTIVELINEMAPS=['+','.join(map(str,range(1,len(zones)+1)))+']','$!VIEW FIT',f'$!SAVELAYOUT "{name}.lay"',f'$!WRITEDATASET "{name}_roundtrip.dat"',' BINARY=NO',' USEPOINTFORMAT=YES',' PRECISION=17',' INCLUDETEXT=NO',' INCLUDEGEOM=NO',' INCLUDECUSTOMLABELS=NO']
    out+=['$!QUIT'];(DEST/'tecplot_roundtrip.mcr').write_text('\n'.join(out)+'\n',encoding='ascii')
def verify():
    assert 'ERROR' not in (DEST/'batch.log').read_text().upper();files=[]
    for name in NAMES:
        source=DEST/(name+'.dat');back=DEST/(name+'_roundtrip.dat');binary=DEST/(name+'.plt');va,za=parse(source);vb,zb=parse(back);assert va==vb==VARS and za==zb and binary.read_bytes()[:8]==b'#!TDV112';files.append(dict(file=binary.name,sha256=sha(binary),source_sha256=sha(source),roundtrip_sha256=sha(back),zones=len(za),points=sum(len(z['rows']) for z in za),all_values_exact=True))
    js(DEST/'tecplot_validation.json',dict(tool='Tecplot 360 EX 2024 R1 / Preplot',files=files));print('NATIVE_VERIFIED',files,flush=True)
def main():
    if '--verify' in sys.argv:verify();return
    DEST.mkdir(parents=True,exist_ok=True);curves,limits,sources=load();items=[]
    for c in curves:
        title=f'{c["case"]}_{c["region"]}_f{c["f_kHz"]:+g}kHz';np.savetxt(DEST/(title+'.csv'),c['rows'],delimiter=',',header=','.join(VARS),comments='',fmt='%.17g');items.append((title,c['rows']))
    rows=[]
    for r in limits:
        rows.append({k:v for k,v in r.items() if k!='record'});a=np.array([[r['x_mm'],r['beta'],r['f_kHz'],r['omega_nd'],r['alpha_r'],0.,0.,r['x_mm']-(0 if r['case']=='bluntplate' else 14.858163),r['xi'],9.,CASES.index(r['case'])+1,r['right_residual'],1.]])
        items.append((r['case']+('_positive_limit' if r['sign']>0 else '_negative_limit'),a))
    with (DEST/'frequency_limits.csv').open('w',newline='') as stream:w=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    js(DEST/'frequency_limits.json',dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),limits=rows,scope='Numerical extinction estimates for the tracked M3 representative and resolved connected pockets; no conjugate copy added',model='Frozen local wall coordinates, zero explicit curvature; original CFD grids',interpolation='Actual reassembled local-profile refinement; not CFD mesh refinement or a grid-convergence uncertainty estimate'))
    dat(NAME,items)
    for case in CASES:dat(case+'_M3_frequency_limits',[(name,a) for name,a in items if name.startswith(case+'_')])
    js(DEST/'selection_manifest.json',dict(curves=sources,curve_count=len(curves),critical_points=len(limits),points=sum(len(c['rows']) for c in curves),model='Unchanged zero-explicit-curvature local-wall LST',no_added_conjugate=True));plot(curves,limits);macros();shutil.copy2(__file__,DEST/'publish_limits.py');print('PUBLISHED',len(curves),'curves',sum(len(c['rows']) for c in curves),'points; four limits',flush=True)
if __name__=='__main__':main()
