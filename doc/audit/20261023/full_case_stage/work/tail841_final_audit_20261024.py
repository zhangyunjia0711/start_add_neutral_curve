#!/usr/bin/env python3
"""Freeze finite tail ridge evidence, native zero witness, tables and scientific plot."""
import csv,json,re,shutil,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from full_driver import ROOT,CASE,sha,snapshot
from local_peak_x_continue import stamp

def main():
    before=snapshot();out=ROOT/'outputs/tail841_finite_ridge_audit_20261024.json'
    if out.exists():raise RuntimeError('Existing final audit preserved')
    paths=[ROOT/('outputs/tail_lower_'+s+'_841_fine_20261024.json') for s in ['beta4','beta6']]
    docs=[json.loads(p.read_text()) for p in paths];assert all(d['finished'] for d in docs)
    bpath=ROOT/'outputs/tail841_all_boundaries_20261024.json';boundaries=json.loads(bpath.read_text())
    # A fresh solve at each distinct neutral root is shared; compute the actual
    # difference to each independent refinement rather than reuse its first value.
    for item in boundaries['results']:
        n=item['actual_native_fresh'];p=item['result']['point']
        item['fresh_alpha_difference']=abs(complex(n['alpha_r']-p['alpha_r'],n['alpha_i']-p['alpha_i']))
    boundaries['derived_difference_audit_time']=stamp();bpath.write_text(json.dumps(boundaries,indent=2)+'\n')
    zpath=ROOT/'outputs/tail841_zero_peak_audit_20261024.json';zero=json.loads(zpath.read_text());assert zero['finished']
    sys.path.insert(0,str(ROOT/'work/flat_neutral_normal_validation_20261023'))
    from validate_flat import solve
    native=solve('tail841_zero_positive_native_20261024',CASE/'M1/in',400,841,docs[1]['last_actual_peak'],
        here=ROOT/'work/tail841_boundary_native_20261024',binary=Path(docs[1]['solver']),sha=docs[1]['solver_sha256'])
    assert native['passed'] and native['growth']>0
    rows=[]
    for d in docs:
        for r in d['rows']:
            p=r['point'];v=r['verification']
            rows.append(dict(source=d['source'],omega=p['omega_r'],beta=p['beta_r'],alpha_r=p['alpha_r'],
                alpha_i=p['alpha_i'],growth=p['growth'],g_beta=r['g_beta'],curvature=r['curvature'],
                g_omega=r['g_omega'],beta_h2_correction=r['beta_h2_correction'],
                beta_prediction_error=r['beta_prediction_error'],verification_run=v['name']))
    table=ROOT/'outputs/tail841_finite_peak_paths_20261024.csv'
    with table.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    fig,axes=plt.subplots(3,1,figsize=(9,9),sharex=True)
    for d,marker in zip(docs,['o','x']):
        ps=[d['upper_seed']]+[r['point'] for r in d['rows']]
        x=[p['omega_r'] for p in ps]
        for ax,key in zip(axes,['growth','beta_r','alpha_r']):
            ax.plot(x,[p[key] for p in ps],marker=marker,ms=3,lw=.7,label=d['source']+' origin')
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
    axes[0].legend();axes[0].set_title('Native i841: finite local-beta-maximum paths')
    for ax in axes:ax.grid(alpha=.2)
    fig.text(.08,.015,'Markers are actual solves; lines guide the eye. No all-beta/all-mode stability claim. Squares: extra low-frequency peaks.',fontsize=8)
    fig.tight_layout(rect=[0,.03,1,1]);figure=ROOT/'outputs/tail841_finite_peak_paths_20261024.png';fig.savefig(figure,dpi=180);plt.close(fig)
    prefixes=['tail_lower_beta4_841_20261024','tail_lower_beta6_841_20261024','tail_lower_beta4_841_fine_20261024',
        'tail_lower_beta6_841_fine_20261024','tail841_beta4_low_lobe_upper_20261024','tail841_beta6_low_lobe_upper_20261024','tail841_zero_audit_20261024']
    runs={p.name:p for prefix in prefixes for p in (ROOT/'runs').glob(prefix+'*') if (p/'summary.json').exists()}
    sums=[json.loads((p/'summary.json').read_text()) for p in runs.values()]
    overlaps=[]
    for p in runs.values():
        log=p/'out/neutral_run.log'
        if log.exists():overlaps.extend(float(x) for x in re.findall(r'vector_overlap=\s*([\d.E+\-]+)',log.read_text()))
    frozen=ROOT/'outputs/tail841_finite_ridge_sources_20261024';frozen.mkdir(exist_ok=False)
    sources=[]
    for name in ['tail_lower_ridge_841.py','tail_lower_ridge_841_fine.py','tail_lower_postprocess_20261024.py',
                 'tail_zero_peak_audit_20261024.py','tail841_final_audit_20261024.py']:
        p=ROOT/'work'/name;q=frozen/name;shutil.copy2(p,q);sources.append(dict(path=str(q),sha256=sha(q)))
    # Only generated Python bytecode from these scripts is unnecessary; preserve
    # all actual solver evidence and the five reproduction sources.
    removed=[]
    for name in ['tail_lower_ridge_841','tail_lower_ridge_841_fine','tail_lower_postprocess_20261024','tail_zero_peak_audit_20261024']:
        p=ROOT/'work/__pycache__'/(name+'.cpython-38.pyc')
        assert p.resolve().is_relative_to(ROOT.resolve()) if hasattr(p,'is_relative_to') else str(p.resolve()).startswith(str(ROOT.resolve())+'/work/')
        if p.exists():p.unlink();removed.append(str(p))
    report=dict(created=stamp(),source_paths=[dict(path=str(p),sha256=sha(p)) for p in paths+[bpath,zpath]],
        finite_paths=[dict(source=d['source'],accepted_frequency_rows=len(d['rows']),reached_zero=d['reached_zero'],
            stop_reason=d['stop_reason'],last_actual_peak=d['last_actual_peak']) for d in docs],
        unique_neutral_boundaries=[dict(direction=k,point=x['result']['point'],h2_J=x['result']['verification']['J'],
            h2_correction=x['result']['verification_correction'],actual_native_fresh=x['actual_native_fresh']) for k,x in unique.items()],
        actual_zero_positive_witness=native,extra_low_frequency_peaks=[dict(omega=r['omega'],peaks=[p.get('point') for p in r.get('refined_peaks',[]) if p['status']=='peak_verified']) for r in zero['rows']],
        qualifications=['Negative growth between the two newly verified boundaries was sampled on these finite paths; no all-beta/all-mode stable-hole proof.',
          'Two source origins recover the same lower/low-lobe-upper roots; provenance is not two physical modes or a proved peak-merger position.',
          'Near zero, two local beta maxima occur at some inserted frequencies. Large last beta displacement does not certify persistence of an old beta stationary-point identity.',
          'Beta4 old-identity continuation remains unknown below its last verified omega. Beta6 reaches an actual positive-growth zero-frequency peak; finite insertions supply other positive-growth witnesses.'],
        cost=dict(unique_cached_run_directories=len(sums),solver_wall_seconds_serial_sum=sum(s['seconds'] for s in sums),
            solves=sum(s['solves'] or 0 for s in sums),rejected=sum(s['rejected'] or 0 for s in sums),
            failed_requested_field_records=sum(s['failed_field_records'] for s in sums),petsc_errors=sum(s['petsc_errors'] for s in sums),
            minimum_recorded_vector_overlap=min(overlaps)),
        original_inputs_unchanged=snapshot()==before,table=str(table),figure=str(figure),sources=sources,
        cleanup_removed=removed,finished=True)
    out.write_text(json.dumps(report,indent=2)+'\n');print('FINAL AUDIT '+json.dumps(dict(output=str(out),cost=report['cost'],original_inputs_unchanged=report['original_inputs_unchanged'])),flush=True)

if __name__=='__main__':main()
