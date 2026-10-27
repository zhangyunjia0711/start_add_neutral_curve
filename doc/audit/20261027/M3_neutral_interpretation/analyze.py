"""Quantify existing neutral-boundary geometry; no new eigenvalue solves."""
import csv,datetime,hashlib,json,shutil
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
BASE=Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate/75°')
SOURCE=BASE/'outputs/M3_xbeta_nine_frequencies_20261027'
DEST=BASE/'outputs/M3_neutral_interpretation_20261027'
FS=[-100,-75,-50,-25,0,25,50,75,100]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def crossing(a,x):
    vals=[]
    for p,q in zip(a[:-1],a[1:]):
        xa,xb=p['x_mm'],q['x_mm']
        if (xa<=x<xb) or (xb<=x<xa):
            t=(x-xa)/(xb-xa);vals.append(float(p['beta']+t*(q['beta']-p['beta'])))
    vals=sorted(vals);roots=[]
    for b in vals:
        if not roots or abs(roots[-1]-b)>1e-8:roots.append(b)
    assert len(roots)==2,roots
    return roots
def main():
    DEST.mkdir(exist_ok=True)
    manifest=json.loads((SOURCE/'delivery_manifest.json').read_text());records=[];provenance=[]
    casepath=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°/neutral_xbeta_four_omega_20261026/outputs/cases.json')
    cases=json.loads(casepath.read_text());flat={c:cases[c]['flat_start_mm'] for c in ['bluntplate','rough']}
    provenance.append(dict(path=str(casepath),sha256=sha(casepath)))
    for f in FS:
        entry={'f_kHz':f}
        for case in ['bluntplate','rough']:
            file=SOURCE/f'{case}_f{f:+d}kHz.csv';r=next(r for r in manifest['files'] if r['path']==file.name);assert sha(file)==r['sha256']
            a=np.genfromtxt(file,delimiter=',',names=True);beta=crossing(a,-25.)
            for k,v in dict(x_onset_mm=float(a['x_mm'].min()),x_extent_mm=float(a['x_mm'].max()),flat_start_mm=flat[case],x_extent_from_flat_start_mm=float(a['x_mm'].max())-flat[case],closed=bool(f),sample_x_mm=-25.,beta_lower=beta[0],beta_upper=beta[1],beta_width=beta[1]-beta[0]).items():entry[case+'_'+k]=v
            provenance.append(dict(path=str(file),sha256=sha(file)))
        entry['delta_x_onset_rough_minus_bluntplate_mm']=entry['rough_x_onset_mm']-entry['bluntplate_x_onset_mm']
        entry['delta_x_close_rough_minus_bluntplate_mm']=entry['rough_x_extent_mm']-entry['bluntplate_x_extent_mm'] if f else None
        entry['delta_x_close_from_flat_start_mm']=entry['rough_x_extent_from_flat_start_mm']-entry['bluntplate_x_extent_from_flat_start_mm'] if f else None
        entry['delta_beta_width']=entry['rough_beta_width']-entry['bluntplate_beta_width']
        entry['beta_width_change_percent']=100*entry['delta_beta_width']/entry['bluntplate_beta_width']
        records.append(entry)
    with (DEST/'neutral_boundary_metrics.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    (DEST/'metrics_provenance.json').write_text(json.dumps(dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),sources=provenance,method='Original solved x extrema; beta crossings at exactly global x=-25 mm by linear interpolation of existing curve segments; no new LST solve',zero_frequency_censored_at_domain_end=True,records=records),ensure_ascii=False,indent=2)+'\n')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
    fig,axes=plt.subplots(2,2,figsize=(12,8.4));fig.subplots_adjust(left=.085,right=.98,bottom=.14,top=.86,hspace=.34,wspace=.27)
    for case,color in [('bluntplate','#2166ac'),('rough','#c44e32')]:
        xend=[r[case+'_x_extent_mm'] if r['f_kHz'] else np.nan for r in records]
        axes[0,0].plot(FS,xend,'o-',color=color,label=case)
        axes[1,0].plot(FS,[r[case+'_beta_width'] for r in records],'o-',color=color,label=case)
    axes[0,0].plot([0],[1100],'^',mfc='none',mec='black',ms=8);axes[0,0].annotate('f = 0: both open at 1100 mm',xy=(0,1100),xytext=(-90,1010),fontsize=9)
    axes[0,0].set(ylabel='Downstream neutral closure x (mm)',xlabel='f (kHz)',title='Where does the neutral domain end?');axes[0,0].legend(frameon=False)
    xs=[r['f_kHz'] for r in records if r['f_kHz']]
    axes[0,1].bar(np.array(xs)-4,[r['delta_x_close_rough_minus_bluntplate_mm'] for r in records if r['f_kHz']],width=8,color='#9567a3',label='Global x difference')
    flatrows=[r for r in records if r['f_kHz'] and r['rough_x_extent_from_flat_start_mm']>0 and r['bluntplate_x_extent_from_flat_start_mm']>0]
    axes[0,1].bar(np.array([r['f_kHz'] for r in flatrows])+4,[r['delta_x_close_from_flat_start_mm'] for r in flatrows],width=8,color='#b8ad69',label='From each flat-section start')
    axes[0,1].axhline(0,color='.5',lw=.5);axes[0,1].legend(frameon=False,fontsize=8)
    axes[0,1].set(ylabel='Closure-position difference\nrough - bluntplate (mm)',xlabel='f (kHz)',title='Reference origin matters; f = 0 excluded')
    axes[1,0].set(ylabel='Unstable beta width at global x = -25 mm',xlabel='f (kHz)',title='Local interval width; values read from curve segments');axes[1,0].legend(frameon=False)
    axes[1,1].bar(FS,[r['beta_width_change_percent'] for r in records],width=14,color='#608a76')
    axes[1,1].set(ylabel='Beta width increase: rough vs bluntplate (%)',xlabel='f (kHz)',title='At global x = -25 mm; interpolation estimate')
    for ax in axes.flat:ax.set_xticks(FS);ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
    fig.suptitle('Quantifying the neutral boundaries of the tracked M3 family',y=.96,fontsize=14)
    fig.text(.5,.90,'Boundary geometry describes where growth is possible. Growth magnitude and N factors are not contained here.',ha='center',fontsize=10)
    fig.text(.5,.035,'Different wall geometries and x grids. Flat-section starts: 0 and 14.858163 mm. A coordinate shift is not a causal decomposition.',ha='center',fontsize=9)
    for ext in ['png','svg']:fig.savefig(DEST/('neutral_boundary_quantitative_comparison.'+ext),dpi=200,facecolor='white')
    plt.close(fig);shutil.copy2(__file__,DEST/'analyze.py')
    print('METRICS',[(r['f_kHz'],round(r['delta_x_close_rough_minus_bluntplate_mm'],3) if r['f_kHz'] else None,round(r['delta_x_close_from_flat_start_mm'],3) if r['f_kHz'] else None,round(r['beta_width_change_percent'],3)) for r in records],flush=True)
if __name__=='__main__':main()
