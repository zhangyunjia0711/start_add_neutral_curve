"""Display the frequency coordinate that the original x-beta projection hides."""
import csv,json,shutil
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
BASE=Path('/home/zaza/attachment_line/Han2019_SweptBluntPlate/75°')
SOURCE=BASE/'neutral_xbeta_four_omega_20261026'
AUDIT=BASE/'neutral_frequency_connection_20261026'
DEST=Path('/mnt/e/work/calc/Instability of attachment line/Han2019_SweptBluntPlate/75°/outputs/neutral_frequency_connection_20261026')

def main():
    report=json.loads((DEST/'connection_report.json').read_text());paths=report['paths'];summary=[]
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
    fig,axes=plt.subplots(2,2,figsize=(12.6,9.2));fig.subplots_adjust(left=.075,right=.985,bottom=.16,top=.89,hspace=.35,wspace=.25)
    blue='#2166ac';purple='#7b3294'
    crossings=[]
    for row,case in enumerate(['bluntplate','rough']):
        ax=axes[row,0];ar=axes[row,1]
        da=json.loads((SOURCE/'outputs'/f'{case}_f50kHz_A_final.json').read_text());db=json.loads((SOURCE/'outputs'/f'{case}_f50kHz_B_final.json').read_text())
        ax.plot([p['x_mm'] for p in da['points']],[p['beta_r'] for p in da['points']],color=blue,lw=1.8,label='A at +50 kHz')
        ax.plot([p['x_mm'] for p in db['points']],[-p['beta_r'] for p in db['points']],color=purple,lw=1.8,ls='--',label='Conjugate of B: now at -50 kHz')
        ax.set_xlim(-35,10);ax.set_ylim(-21,5);ax.set_xlabel('Global x (mm)');ax.set_ylabel(r'$\beta$');ax.grid(alpha=.2)
        ax.set_title(f'{case}: x-beta projection, nose detail')
        # Explicitly identify approximate projected intersections; they are not
        # intersections in full (x,beta,f) space.
        natives={}
        for key,d in [('A',da),('Bc',db)]:
            group={}
            for p in d['points']:
                if p['kind']!=0:continue
                group.setdefault(int(p['xi']),[]).append(p)
            natives[key]={i:sorted([(p['beta_r'] if key=='A' else -p['beta_r'],p['x_mm']) for p in ps]) for i,ps in group.items()}
        stations=sorted(set(natives['A'])&set(natives['Bc']))
        for ka in [0,1]:
            for kb in [0,1]:
                for ia,ib in zip(stations[:-1],stations[1:]):
                    if ib!=ia+1:continue
                    aa=natives['A'][ia][ka];ab=natives['A'][ib][ka];ba=natives['Bc'][ia][kb];bb=natives['Bc'][ib][kb]
                    d0=aa[0]-ba[0];d1=ab[0]-bb[0]
                    if d0*d1<0:
                        t=-d0/(d1-d0);x=aa[1]+t*(ab[1]-aa[1]);b=aa[0]+t*(ab[0]-aa[0])
                        crossings.append(dict(case=case,estimate_only=True,x_mm=x,beta=b,native_bracket=[ia,ib],A_boundary=ka,conjugate_B_boundary=kb,frequency_A=50,frequency_conjugate_B=-50))
                        ax.plot(x,b,'x',color='#222',ms=5,mew=1.)
        selected=sorted([p for p in paths if p['case']==case],key=lambda p:p['side'])
        for p in selected:
            records=p['records'];f=np.array([r['f_kHz'] for r in records]);b=np.array([r['beta'] for r in records])
            ar.plot(f,b,color='#444',lw=1.6,label=('Lower neutral beta' if p['side']=='lower' else 'Upper neutral beta'))
            ar.plot(f[0],b[0],'o',color=blue,ms=6);ar.plot(f[-1],b[-1],'s',color=purple,ms=6)
            checkpoints=[records[c['path_index']] for c in p['original_checks'] if abs(c['f_kHz'])<49]
            ar.plot([v['f_kHz'] for v in checkpoints],[v['beta'] for v in checkpoints],'o',mfc='white',mec='#444',ms=4)
            with (DEST/f'{case}_{p["side"]}_continuous_neutral_path.csv').open('w',newline='') as stream:
                names=['f_kHz','omega_nd','beta','alpha_r','alpha_i','growth','x_mm','right_residual','previous_vector_MAC'];w=csv.DictWriter(stream,fieldnames=names,extrasaction='ignore',lineterminator='\n');w.writeheader();w.writerows(records)
            summary.append(dict(case=case,side=p['side'],x_mm=p['x_mm'],records=len(records),maximum_growth_error=p['maximum_growth_error'],minimum_neighbor_MAC=p['minimum_neighbor_MAC'],endpoint=p['negative_endpoint_matches_conjugate_B'],checkpoints=[dict(f_kHz=c['f_kHz'],beta=c['beta']) for c in p['original_checks']]))
        ar.axvline(0.,color='#aaa',lw=.8);ar.set_xlim(-53,53);ar.set_ylim(-21,5);ar.grid(alpha=.2);ar.set_xlabel('Signed physical frequency f (kHz)');ar.set_ylabel(r'Neutral $\beta$')
        ar.set_title(f'Continuous neutral paths at x = {selected[0]["x_mm"]:.3f} mm')
        ar.text(.04,.96,'Both boundaries stay neutral while f varies',transform=ar.transAxes,va='top',fontsize=9)
    fig.suptitle('A(+50 kHz) connects to the conjugate of B(+50 kHz) through f = 0',y=.97,fontsize=13)
    handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,.043),ncol=2,fontsize=9,frameon=False)
    fig.text(.5,.018,'Left: different frequency planes projected together; crosses mark estimated projection intersections. Right: verified neutral continuation.',ha='center',fontsize=9)
    for ext in ['png','svg']:fig.savefig(DEST/f'neutral_frequency_connection.{ext}',dpi=220,facecolor='white')
    plt.close(fig)
    (DEST/'projection_crossings.json').write_text(json.dumps(crossings,indent=2)+'\n')
    (DEST/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    shutil.copy2(__file__,DEST/'publish.py');print('PUBLISHED',len(crossings),'estimated projected crossings',flush=True)
if __name__=='__main__':main()
