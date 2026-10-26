"""Scientific plots and Tecplot ASCII inputs from audited physical points."""
import csv,json,shutil
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import ROOT,CASES,snapshot

DEST=Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate/75°/outputs/M3_omega0_xbeta_20261026')
COLORS={'bluntplate':'#2166ac','rough':'#c44e32'}
VARS=['x_mm','beta','omega','alpha_r','alpha_i','growth','x_flat_mm','station_coordinate','point_kind','case_id','conjugate','right_residual','raw_MAC']

def array(case,points,conjugate=False):
    out=[]
    for p in points:
        out.append([p['x_mm'],-p['beta_r'] if conjugate else p['beta_r'],0.,-p['alpha_r'] if conjugate else p['alpha_r'],p['alpha_i'],p['growth'],p['x_flat_mm'],p['xi'],p['kind'],1 if case=='bluntplate' else 2,int(conjugate),p['right_residual'],p['raw_MAC']])
    return np.array(out)

def bands(data,coordinate='x_mm'):
    points=[q for q in data['points'] if q['kind']==0];rows=[]
    for i in sorted({int(q['xi']) for q in points}):
        a=sorted([q for q in points if int(q['xi'])==i],key=lambda q:q['beta_r']);assert len(a)==2
        rows.append([a[0][coordinate],a[0]['beta_r'],a[1]['beta_r']])
    return np.array(rows)

def curves(ax,items,coordinate='x_mm',signed=False,fill=True):
    for case,data in items:
        p=data['points'];x=np.array([q[coordinate] for q in p]);b=np.array([q['beta_r'] for q in p]);color=COLORS[case]
        ax.plot(x,b,color=color,lw=1.8,label=case)
        if fill:
            bd=bands(data,coordinate);ax.fill_between(bd[:,0],bd[:,1],bd[:,2],color=color,alpha=.09)
        for tip in data['physical_tips']:ax.plot(tip[coordinate],tip['beta_r'],'o',color=color,ms=4)
        if signed:
            ax.plot(x,-b,color=color,lw=1.4,ls='--')
            if fill:ax.fill_between(bd[:,0],-bd[:,1],-bd[:,2],color=color,alpha=.06)
    ax.set_xlabel(r'Global $x$ (mm)' if coordinate=='x_mm' else r'$x-x_{\mathrm{flat}}$ (mm)')
    ax.set_ylabel(r'$\beta$ (nondimensional)');ax.grid(alpha=.2);ax.set_axisbelow(True)

def save(fig,name):
    fig.savefig(DEST/(name+'.png'),dpi=220,facecolor='white');fig.savefig(DEST/(name+'.svg'),facecolor='white');plt.close(fig)

def dat(name,items):
    lines=['TITLE="M3 stationary neutral boundaries: omega=0"','VARIABLES='+','.join('"'+x+'"' for x in VARS)]
    for case,data,conjugate in items:
        values=array(case,data['points'],conjugate);zone=case+('_beta_positive_conjugate' if conjugate else '_beta_negative')
        lines.append(f'ZONE T="{zone}", I={len(values)}, DATAPACKING=POINT')
        lines.append('DT=('+(' '.join(['DOUBLE']*len(VARS)))+')')
        lines.extend(' '.join('%.17g'%v for v in row) for row in values)
    (DEST/(name+'.dat')).write_text('\n'.join(lines)+'\n',encoding='ascii')

def main():
    DEST.mkdir(exist_ok=True);snapshot();datasets={case:json.loads((ROOT/'outputs'/f'{case}_final.json').read_text()) for case in CASES};items=list(datasets.items())
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
    for name,selection in [('M3_omega0_xbeta_comparison',items)]+[(case+'_M3_omega0_xbeta',[(case,data)]) for case,data in items]:
        fig,ax=plt.subplots(1,2,figsize=(12.5,5.2),gridspec_kw={'width_ratios':[1.45,1]})
        fig.subplots_adjust(left=.07,right=.985,bottom=.18,top=.82,wspace=.24)
        for a in ax:curves(a,selection)
        ax[0].set_xlim(-35,1100);ax[0].set_ylim(-16,0);ax[0].set_title('(a) Full computed domain');ax[0].legend(loc='lower right',frameon=False)
        ax[1].set_xlim(-35,45);ax[1].set_ylim(-16,0);ax[1].set_title('(b) Nose and upstream region')
        fig.suptitle(r'M3 neutral boundaries at $\omega=0$',y=.97)
        fig.text(.5,.035,'Shading: M3 unstable band. The boundary is open at x = 1100 mm; no artificial closure.',ha='center',fontsize=9,color='#444444')
        save(fig,name)
    fig,ax=plt.subplots(figsize=(10,5.4),layout='constrained');curves(ax,items,signed=True);ax.set_xlim(-35,1100);ax.set_ylim(-16,16);ax.legend(frameon=False);ax.set_title(r'$\omega=0$: negative-$\beta$ branch and verified positive-$\beta$ conjugate');save(fig,'M3_omega0_xbeta_signed')
    fig,ax=plt.subplots(1,2,figsize=(12.5,4.7),layout='constrained')
    for a in ax:curves(a,items,'x_flat_mm')
    ax[0].set_xlim(-50,1100);ax[1].set_xlim(-50,45)
    for a in ax:a.set_ylim(-16,0)
    ax[0].legend(frameon=False);fig.suptitle(r'$\omega=0$: coordinates relative to each flat-plate origin');save(fig,'M3_omega0_xbeta_flat_origin')
    combined=[]
    for case,data in items:
        for conj in [False,True]:
            values=array(case,data['points'],conj)
            with (DEST/(case+('_positive_conjugate' if conj else '')+'_M3_omega0_xbeta.csv')).open('w',newline='') as f:
                w=csv.writer(f);w.writerow(VARS);w.writerows(values)
            combined.append((case,data,conj))
        dat(case+'_M3_omega0_xbeta',[(case,data,False),(case,data,True)])
        shutil.copy2(ROOT/'outputs'/f'{case}_final.json',DEST/f'{case}_calculation_evidence.json')
    # Negative branches first so that Tecplot maps 1 and 2 directly compare cases.
    dat('M3_omega0_xbeta_comparison',[(c,d,False) for c,d in items]+[(c,d,True) for c,d in items])
    summary=[]
    for case,data in items:
        native=[q for q in data['points'] if q['kind']==0];tail=sorted([q for q in native if int(q['xi'])==CASES[case]['stations']],key=lambda q:abs(q['beta_r']))
        tip=data['physical_tips'][0];maxb=min(native,key=lambda q:q['beta_r'])
        summary.append(dict(case=case,station_count=data['native_station_count'],native_points=data['native_point_count'],profile_tips=len(data['physical_tips']),x_min_mm=tip['x_mm'],beta_at_xmin=tip['beta_r'],x_max_mm=max(q['x_mm'] for q in native),
            tail_abs_beta_band=[abs(q['beta_r']) for q in tail],maximum_abs_beta=abs(maxb['beta_r']),x_at_maximum_abs_beta_mm=maxb['x_mm'],maximum_growth_error=data['maximum_growth_error'],maximum_right_residual=data['maximum_right_residual'],status=data['status']))
    (DEST/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');(ROOT/'outputs/summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    for name in ['input_manifest.json','cases.json','pilot.json','preparation.json','band_audit.json']:
        if (ROOT/'outputs'/name).exists():shutil.copy2(ROOT/'outputs'/name,DEST/name)
    scripts=DEST/'reproduction';scripts.mkdir(exist_ok=True)
    for name in ['bootstrap.py','common.py','run_adapter.py','solver.py','prepare.py','contour.py','physical.py','finalize.py','publish.py','audit.py']:
        shutil.copy2(ROOT/'work'/name,scripts/name)
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':main()
