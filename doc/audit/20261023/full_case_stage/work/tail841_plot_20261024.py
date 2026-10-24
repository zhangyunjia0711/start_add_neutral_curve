#!/usr/bin/env python3
"""Final actual-point-only plot, preserving unknown beta-stationary identities."""
import json,shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from full_driver import ROOT,sha
from local_peak_x_continue import stamp

def main():
    docs=[json.loads((ROOT/('outputs/tail_lower_'+s+'_841_fine_20261024.json')).read_text()) for s in ['beta4','beta6']]
    boundaries=json.loads((ROOT/'outputs/tail841_all_boundaries_20261024.json').read_text())
    zero=json.loads((ROOT/'outputs/tail841_zero_peak_audit_20261024.json').read_text())
    fig,axes=plt.subplots(3,1,figsize=(9,9),sharex=True)
    for doc,marker in zip(docs,['o','x']):
        points=[doc['upper_seed']]+[r['point'] for r in doc['rows']]
        for ax,key in zip(axes,['growth','beta_r','alpha_r']):
            ax.scatter([p['omega_r'] for p in points],[p[key] for p in points],marker=marker,s=11,label=doc['source']+' origin')
    for row in zero['rows']:
        for peak in row.get('refined_peaks',[]):
            if peak['status']!='peak_verified':continue
            p=peak['point']
            for ax,key in zip(axes,['growth','beta_r','alpha_r']):
                ax.scatter(p['omega_r'],p[key],s=18,facecolors='none',edgecolors='black',marker='s')
    unique={x['direction']:x for x in boundaries['results']}
    for item in unique.values():
        p=item['result']['point']
        for ax in axes:ax.axvline(p['omega_r'],c='gray',ls=':',lw=.8)
    axes[0].axhline(0,c='black',lw=.7);axes[0].set_ylabel('Spatial growth -Im(alpha)')
    axes[1].set_ylabel('Local peak beta');axes[2].set_ylabel('Re(alpha)');axes[2].set_xlabel('Real omega (nondimensional)')
    axes[0].legend();axes[0].set_title('Native i841: finite local-beta-maximum samples')
    for ax in axes:ax.grid(alpha=.2)
    fig.text(.08,.015,'Actual solved points only; peak identities and uncomputed intervals are not connected. Squares: extra low-frequency peaks.',fontsize=8)
    fig.tight_layout(rect=[0,.03,1,1]);figure=ROOT/'outputs/tail841_finite_peak_paths_20261024.png';fig.savefig(figure,dpi=180);plt.close(fig)
    source=ROOT/'work/tail841_plot_20261024.py';frozen=ROOT/'outputs/tail841_finite_ridge_sources_20261024/tail841_plot_20261024.py'
    assert not frozen.exists();shutil.copy2(source,frozen)
    p=ROOT/'outputs/tail841_finite_ridge_audit_20261024.json';audit=json.loads(p.read_text())
    audit['sources'].append(dict(path=str(frozen),sha256=sha(frozen)))
    audit['figure_plot_policy']='Actual solved scatter points only. Do not connect unknown stationary-point identity or uncomputed frequency intervals.'
    audit['figure_visual_review']=stamp();p.write_text(json.dumps(audit,indent=2)+'\n')

if __name__=='__main__':main()
