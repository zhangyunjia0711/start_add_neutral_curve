"""Scientific figures and Tecplot ASCII; 75-degree data are read-only reuse."""
import csv,datetime,hashlib,json,os,shutil
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
WIN=Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate')
ROOT=Path(os.environ.get('M3_SWEEP_ROOT','/home/zaza/attachment_line/Han2019_SweptBluntPlate/neutral_m3_sweep_consistent_20261028'))
DEST=Path(os.environ.get('M3_SWEEP_DEST',str(WIN/'outputs/M3_sweep_angles_20261028')))
OLD75=Path(os.environ.get('M3_SWEEP_OLD75',str(WIN/'75°/outputs/M3_all_neutral_curves_20261027')))
NEW_ANGLES=[30,45,50,55,60,65,70,73]
ANGLES=NEW_ANGLES+[75]
BASEVARS=['x_mm','beta','f_kHz','omega_nd','alpha_r','alpha_i','growth','x_flat_mm','station_coordinate','point_kind','case_id','right_residual','raw_MAC']
VARS=BASEVARS+['sweep_deg','normal_method']
COLORS={-240.2:'#0a3154',-237.14:'#385a7c',-200.:'#007b8a',-150.:'#498d6d',-112.:'#00a6a6',-110.:'#4c63d2',-100.:'#762a83',-75.:'#5954a6',-50.:'#2166ac',-25.:'#149fb3',0.:'#161616',25.:'#658e28',50.:'#1b9e77',75.:'#e48120',100.:'#d73027',125.:'#8c564b',150.:'#e377c2',175.:'#bc8f00',200.:'#7f3c8d',230.:'#006d77',250.:'#a44a3f',264.53:'#a52a6a',275.:'#d29d00',290.:'#e65b00',291.6:'#c51b7d'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def js(p,obj):p.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')
def load():
    (DEST/'curves').mkdir(parents=True,exist_ok=True);curves=[]
    for angle in NEW_ANGLES:
        for p in sorted((ROOT/'outputs').glob(f's{angle}_f*kHz_*_final.json')):
            d=json.loads(p.read_text());assert d['status'] in ['closed_verified','open_at_domain_boundaries']
            rows=[]
            for r in d['points']:
                assert abs(r['growth'])<1e-8 and r['right_residual']<1e-7
                rows.append([r['x_mm'],r['beta_r'],d['f_kHz'],d['omega_r'],r['alpha_r'],r['alpha_i'],r['growth'],r['x_flat_mm'],r['xi'],r['kind'],1,r['right_residual'],r['raw_MAC'],angle,1])
            a=np.array(rows);assert np.isfinite(a).all()
            path=DEST/'curves'/(p.stem+'.csv')
            with path.open('w',newline='') as f:w=csv.writer(f);w.writerow(VARS);w.writerows(a)
            curves.append(dict(angle=angle,f_kHz=d['f_kHz'],component=d['component'],zone=p.stem,rows=a,output_csv=str(path.relative_to(DEST)),source=dict(path=str(p),sha256=sha(p)),closed=d['status']=='closed_verified',method='wall_local_consistent',original_check_count=len(d['original_checks']),physical_tip_count=len(d['physical_tips'])))
    # No solver import or operation occurs for 75 degrees. Preserve every old curve byte-for-byte.
    selection=json.loads((OLD75/'selection_manifest.json').read_text());manifest=json.loads((OLD75/'delivery_manifest.json').read_text())
    for r in selection['curves']:
        if r['case']!='bluntplate':continue
        p=OLD75/r['output_csv'];expected=next(q['sha256'] for q in manifest['files'] if q['path']==r['output_csv']);assert sha(p)==expected
        with p.open() as f:
            reader=csv.DictReader(f);assert reader.fieldnames==BASEVARS;a=np.array([[float(row[k]) for k in BASEVARS] for row in reader])
        a=np.column_stack([a,np.full(len(a),75.),np.zeros(len(a))]);assert len(a)==r['points']
        path=DEST/'curves'/('s75_existing_'+p.name);shutil.copy2(p,path);assert sha(path)==sha(p)
        curves.append(dict(angle=75,f_kHz=r['f_kHz'],component=r['region'],zone='s75_existing_'+r['zone'],rows=a,output_csv=str(path.relative_to(DEST)),source=dict(path=str(p),sha256=sha(p)),closed=r['closed'],method='existing_75_legacy_wall_local',original_check_count=None,physical_tip_count=None))
    return sorted(curves,key=lambda r:(r['angle'],r['f_kHz'],r['component']))
def dat(name,curves):
    p=DEST/(name+'.dat')
    with p.open('w',encoding='ascii') as f:
        f.write('TITLE="M3 x-beta neutral curves versus sweep angle"\nVARIABLES='+','.join('"'+v+'"' for v in VARS)+'\n')
        for r in curves:
            f.write(f'ZONE T="{r["zone"]}", I={len(r["rows"])}, J=1, K=1, DATAPACKING=POINT\nDT=('+','.join(['DOUBLE']*len(VARS))+')\n')
            for row in r['rows']:f.write(' '.join(format(v,'.17g') for v in row)+'\n')
    return p
def plot(curves,shared_beta=False):
    allbeta=np.concatenate([r['rows'][:,1] for r in curves]);global_ylim=(float(allbeta.min()-3),float(allbeta.max()+3))
    freqs=sorted(set(r['f_kHz'] for r in curves));handles=[Line2D([0],[0],color=COLORS.get(f,'#888888'),lw=2 if f==0 else 1.4,ls='-' if f>=0 else '--',label=f'{f:+g}' if f else '0') for f in freqs]
    def axis(ax,angle,window,common=False):
        selected=[r for r in curves if r['angle']==angle]
        if selected:
            values=np.concatenate([r['rows'][:,1] for r in selected]);pad=max(1.,.04*(values.max()-values.min()));ylim=(values.min()-pad,values.max()+pad)
        else:ylim=(-1,1)
        for r in sorted(selected,key=lambda r:r['f_kHz']==0):
            a=r['rows'];ax.plot(a[:,0],a[:,1],color=COLORS.get(r['f_kHz'],'#888888'),lw=1.3 if r['f_kHz'] else 2.,ls='-' if r['f_kHz']>=0 else '--')
        ax.set_xlim(window);ax.set_ylim(global_ylim if shared_beta else ylim);ax.grid(alpha=.2);ax.set_xlabel('Global x (mm)');ax.set_ylabel(r'$\beta$ (nondimensional)')
        suffix=' / retained data' if angle==75 else ''
        ax.set_title(f'{angle}°'+suffix+(' / common 1.1 m window' if common else ' / upstream detail'),fontsize=11)
        if not selected:ax.text(.5,.5,'Calculation pending',ha='center',transform=ax.transAxes)
    footer='New calculations: current 1.1m FPP fields with consistent normal derivatives. 75°: original delivered curves, reused without recomputation.\n70° input extends to 3.5 m; comparison panels show the common 1.1 m window. Input processing/discretization differ from retained 75° data.'
    fig,axes=plt.subplots(6,3,figsize=(20,24))
    for index,angle in enumerate(ANGLES):
        row=(index//3)*2;col=index%3;axis(axes[row,col],angle,(-36,1105),True);axis(axes[row+1,col],angle,(-36,200))
    fig.legend(handles=handles,title='Signed cyclic frequency f (kHz); same colors in every angle',loc='lower center',bbox_to_anchor=(.5,.027),ncol=9,fontsize=9,title_fontsize=10)
    fig.suptitle('M3 x–β neutral curves at different sweep angles',fontsize=19,y=.993);fig.text(.5,.009,footer,ha='center',fontsize=9)
    fig.subplots_adjust(left=.06,right=.985,top=.97,bottom=.096,hspace=.43,wspace=.25)
    stem='M3_sweep_angles_shared_beta_scale' if shared_beta else 'M3_sweep_angles_comparison'
    for ext in ['png','svg']:fig.savefig(DEST/(stem+'.'+ext),dpi=160)
    plt.close(fig)
    if shared_beta:return
    for detail in [False,True]:
        fig,axes=plt.subplots(3,3,figsize=(20,15))
        for index,angle in enumerate(ANGLES):axis(axes.flat[index],angle,(-36,200) if detail else (-36,1105),not detail)
        fig.suptitle('M3 fixed-frequency x–β neutral curves / '+('upstream comparison' if detail else 'common 1.1 m window'),fontsize=19,y=.986)
        fig.legend(handles=handles,title='Signed cyclic frequency f (kHz)',loc='lower center',bbox_to_anchor=(.5,.035),ncol=9,fontsize=9,title_fontsize=10);fig.text(.5,.009,footer,ha='center',fontsize=9)
        fig.subplots_adjust(left=.06,right=.985,top=.94,bottom=.15,hspace=.36,wspace=.25)
        stem='M3_sweep_angles_upstream' if detail else 'M3_sweep_angles_common_1p1m'
        for ext in ['png','svg']:fig.savefig(DEST/(stem+'.'+ext),dpi=160)
        plt.close(fig)
    for angle in ANGLES:
        selected=[r for r in curves if r['angle']==angle]
        if not selected:continue
        xmax=max(1100,max(r['rows'][:,0].max() for r in selected));fig,axes=plt.subplots(1,2,figsize=(15,7))
        axis(axes[0],angle,(-36,xmax+5),False);axes[0].set_title(f'{angle}° / complete computed extent'+(' / retained data' if angle==75 else ''))
        local_freqs=sorted(set(r['f_kHz'] for r in selected))
        local_handles=[Line2D([0],[0],color=COLORS.get(f,'#888888'),lw=2 if f==0 else 1.4,ls='-' if f>=0 else '--',label=f'{f:+g}' if f else '0') for f in local_freqs]
        axis(axes[1],angle,(-36,200));fig.legend(handles=local_handles,loc='lower center',bbox_to_anchor=(.5,.015),ncol=9,fontsize=8,title='f (kHz)');fig.subplots_adjust(left=.07,right=.985,top=.9,bottom=.24,wspace=.25)
        fig.suptitle(f'M3 x–β neutral curves at {angle}°',fontsize=16);fig.savefig(DEST/f'M3_{angle}deg_xbeta.png',dpi=150);plt.close(fig)

def main():
    curves=load();assert curves
    for name,subset in [('M3_sweep_angles_all_neutral_curves',curves)]+[(f'M3_{a}deg_xbeta',[r for r in curves if r['angle']==a]) for a in ANGLES]:
        if subset:dat(name,subset)
    plot(curves);plot(curves,True)
    records=[]
    for r in curves:
        a=r['rows'];records.append(dict(**{k:v for k,v in r.items() if k!='rows'},points=len(a),x_range_mm=[float(a[:,0].min()),float(a[:,0].max())],beta_range=[float(a[:,1].min()),float(a[:,1].max())],max_growth_error=float(np.max(abs(a[:,6]))),max_right_residual=float(a[:,11].max())))
    js(DEST/'selection_manifest.json',dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),curves=records,method_codes={'0':'Retained original 75-degree wall_local result','1':'Current FPP fields with consistent projected-normal derivatives'},new_angles=NEW_ANGLES,reused_angle=75,pending_angles=[63,68,80]))
    for name in ['frequency_conversion.json','consistent_validation.json','build_receipt.json','consistent_normal_derivatives.patch','consistent_mode_complete.patch','source_application.json','final_numerical_audit.json','provenance_corrections.json']:shutil.copy2(ROOT/'outputs'/name,DEST/name)
    inputs=json.loads((ROOT/'outputs/input_selection.json').read_text());inputs['selected']=[r for r in inputs['selected'] if r['angle'] in NEW_ANGLES]
    inputs['retained_75']=dict(source_directory=str(OLD75),selection_manifest_sha256=sha(OLD75/'selection_manifest.json'),delivery_manifest_sha256=sha(OLD75/'delivery_manifest.json'),recomputed=False,csv_copies_verified=True)
    js(DEST/'input_selection.json',inputs)
    frequency_status=[]
    for angle in ANGLES:
        for f in sorted(COLORS):
            chosen=[r for r in records if r['angle']==angle and r['f_kHz']==f]
            frequency_status.append(dict(sweep_deg=angle,f_kHz=f,status=('reused_verified' if angle==75 else 'new_verified') if chosen else ('not_in_retained_75_delivery' if angle==75 else 'unresolved_no_delivered_curve'),components=len(chosen),points=sum(r['points'] for r in chosen)))
    with (DEST/'frequency_status.csv').open('w',newline='') as file:
        writer=csv.DictWriter(file,fieldnames=list(frequency_status[0]));writer.writeheader();writer.writerows(frequency_status)
    js(DEST/'frequency_scope.json',dict(requested_common_signed_frequencies_kHz=sorted(COLORS),frequency_status=frequency_status,scope='Every verified M3 component delivered by this calculation; unsuccessful seed search is unresolved, not proof of stability. Existing 75-degree data alone do not supply absent frequencies.',pending_fpp_angles=[63,68,80]))
    js(DEST/'plot_validation.json',dict(curves=len(curves),points=sum(len(r['rows']) for r in curves),angles={str(a):sum(r['angle']==a for r in curves) for a in ANGLES},original75_csv_exact_copies=all(sha(DEST/r['output_csv'])==r['source']['sha256'] for r in curves if r['angle']==75),no_conjugate_copies_added=True,neutral_curves_only=True,plotted_zone_ids=[r['zone'] for r in records],comparison_x_window_mm=[-36,1105],upstream_x_window_mm=[-36,200],per_angle_figures_include_complete_computed_extent=True))
    print('PUBLISHED',len(curves),sum(len(r['rows']) for r in curves),DEST,flush=True)
if __name__=='__main__':main()
