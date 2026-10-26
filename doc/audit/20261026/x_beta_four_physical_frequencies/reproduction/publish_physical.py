"""Publish only corrected native points and actual-profile fold points."""
import csv,json,shutil
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from common import ROOT,snapshot

DEST=Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate/75°/outputs/M3_xbeta_four_omega_20261026')
COLORS={'bluntplate':'#2166ac','rough':'#c44e32'}
VARS=['x_mm','beta','omega_nd','f_kHz','alpha_r','alpha_i','growth','x_flat_mm','station_coordinate','point_kind','case_id','component_id','conjugate','right_residual','raw_MAC']
FREQUENCIES=[-100,-50,50,100]

def values(case,f,c,data):
    sign=1 if f>0 else -1
    return np.array([[p['x_mm'],sign*p['beta_r'],sign*p['omega_r'],f,sign*p['alpha_r'],p['alpha_i'],p['growth'],p['x_flat_mm'],p['xi'],p['kind'],1 if case=='bluntplate' else 2,1 if c=='A' else 2,int(f<0),p['right_residual'],p['raw_MAC']] for p in data['points']])

def curve(ax,items,f,zoom=False):
    for case,c,data in items:
        a=values(case,f,c,data);ax.plot(a[:,0],a[:,1],color=COLORS[case],lw=1.65,ls='-' if c=='A' else '--')
    ax.set_xlabel(r'Global $x$ (mm)');ax.set_ylabel(r'$\beta$ (nondimensional)');ax.grid(alpha=.2);ax.set_axisbelow(True)
    if zoom:ax.set_xlim(-35,5)
    ax.set_title(f'$f={f:+d}$ kHz; '+r'$\omega_{nd}=$'+f'{np.sign(f)*abs(items[0][2]["omega_r"]):.6f}',fontsize=11)

def save(fig,name):
    fig.savefig(DEST/(name+'.png'),dpi=220,facecolor='white');fig.savefig(DEST/(name+'.svg'),facecolor='white');plt.close(fig)

def dat(name,items):
    lines=['TITLE="M3 neutral curves at signed physical frequencies"','VARIABLES='+','.join('"'+s+'"' for s in VARS)]
    for case,f,c,data in items:
        a=values(case,f,c,data);zone=f'{case}_f{f:+d}kHz_{c}'+('_conjugate' if f<0 else '')
        lines += [f'ZONE T="{zone}", I={len(a)}, DATAPACKING=POINT','DT=('+' '.join(['DOUBLE']*len(VARS))+')']
        lines += [' '.join('%.17g'%v for v in row) for row in a]
    (DEST/(name+'.dat')).write_text('\n'.join(lines)+'\n',encoding='ascii')

def main():
    DEST.mkdir(exist_ok=True);snapshot()
    data={(case,f,c):json.loads((ROOT/'outputs'/f'{case}_f{f}kHz_{c}_final.json').read_text()) for case in COLORS for f in [50,100] for c in ['A','B']}
    assert all(d['negative_frequency_checks'] for d in data.values())
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
    legend=[Line2D([0],[0],color=v,lw=2,label=k) for k,v in COLORS.items()]+[Line2D([0],[0],color='#444',ls=s,label=f'M3 component {c}') for c,s in [('A','-'),('B','--')]]
    allitems=[]
    for f in FREQUENCIES:
        items=[(case,c,data[case,abs(f),c]) for case in COLORS for c in ['A','B']]
        fig,ax=plt.subplots(1,2,figsize=(12.8,5.1));fig.subplots_adjust(left=.07,right=.985,bottom=.19,top=.83,wspace=.23)
        curve(ax[0],items,f);curve(ax[1],items,f,True)
        ax[0].text(.98,.97,'Full traced closed components',transform=ax[0].transAxes,fontsize=9,ha='right',va='top')
        ax[1].text(.98,.97,'Nose detail',transform=ax[1].transAxes,fontsize=9,ha='right',va='top')
        fig.legend(handles=legend,loc='lower center',bbox_to_anchor=(.5,.025),ncol=4,frameon=False)
        fig.suptitle('Two-case M3 neutral curves: fixed physical frequency',y=.97)
        name=f'M3_xbeta_f{f:+d}kHz_comparison';save(fig,name)
        dat(name,[(case,f,c,d) for case,c,d in items]);allitems += [(case,f,c,d) for case,c,d in items]
        for case,c,d in items:
            with (DEST/f'{case}_f{f:+d}kHz_{c}.csv').open('w',newline='') as stream:
                w=csv.writer(stream,lineterminator='\n');w.writerow(VARS);w.writerows(values(case,f,c,d))
    fig,axes=plt.subplots(2,2,figsize=(12.8,9.3));fig.subplots_adjust(left=.075,right=.985,bottom=.13,top=.9,wspace=.2,hspace=.3)
    for ax,f in zip(axes.flat,FREQUENCIES):curve(ax,[(case,c,data[case,abs(f),c]) for case in COLORS for c in ['A','B']],f)
    fig.legend(handles=legend,loc='lower center',bbox_to_anchor=(.5,.025),ncol=4,frameon=False)
    fig.suptitle('M3 neutral curves at f = -100, -50, +50, +100 kHz',y=.975)
    save(fig,'M3_xbeta_four_frequencies_comparison');dat('M3_xbeta_four_frequencies',allitems)
    for case in COLORS:dat(case+'_M3_xbeta_four_frequencies',[i for i in allitems if i[0]==case])
    for (case,f,c),d in data.items():shutil.copy2(ROOT/'outputs'/f'{case}_f{f}kHz_{c}_final.json',DEST/f'{case}_f{f}kHz_{c}_calculation_evidence.json')
    for name in ['frequency_conversion.json','physical_audit.json','input_manifest.json','cases.json','physical_pilot.json','physical_survey_summary.json','contour_summary.json','finalize_summary.json','preparation.json']:
        shutil.copy2(ROOT/'outputs'/name,DEST/name)
    repro=DEST/'reproduction';repro.mkdir(exist_ok=True)
    for name in ['bootstrap.py','common.py','run_adapter.py','core.py','solver.py','prepare.py','physical_pilot.py','physical_survey.py','survey.py','contour.py','contour3d.py','physical.py','tip_correction.py','finalize.py','run_contours.py','run_finalize.py','audit_physical.py','publish_physical.py']:
        shutil.copy2(ROOT/'work'/name,repro/name)
    summary=[dict(case=case,f_kHz=f,component=c,status=d['status'],x_range_mm=d['x_range_mm'],beta_range=d['negative_beta_range'],native_points=d['native_point_count'],physical_tips=len(d['physical_tips']),growth_error=d['maximum_growth_error'],right_residual=d['maximum_right_residual']) for (case,f,c),d in data.items()]
    (DEST/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':main()
