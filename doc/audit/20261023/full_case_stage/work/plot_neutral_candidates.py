#!/usr/bin/env python3
"""Export actual continued candidates and explicit coverage; never fill gaps."""
import argparse,csv,datetime,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from full_driver import ROOT,sha
from centered_survey import coarse_stations


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--phase',action='append',required=True)
    p.add_argument('--output',default='neutral_candidate_progress')
    a=p.parse_args();all_points=[];chains=[];inputs=[];failures=[];sources=[]
    for phase in a.phase:
        for f in sorted((ROOT/'outputs').glob(phase+'_*_progress.json')):
            d=json.loads(f.read_text());spec=d.get('spec',{})
            if 'seed_provenance' not in spec:continue
            source=spec['seed_provenance']['source'];sources.append(source)
            points=[dict(spec['seed'],station=spec['seed_station'],source=source,chain=spec['name'],
                kind='actual_seed',solver_sha256=spec['solver_sha256'],result_file=str(f))]
            for r in d['records']:
                if r.get('status')!='candidate_converged':continue
                v=dict(r['point'],station=r['station'],source=source,chain=spec['name'],kind='refined_candidate',
                    solver_sha256=spec['solver_sha256'],result_file=str(f),
                    g_beta=r['result']['derivative'],g_beta_beta=r['result']['curvature'],
                    fd_h=r['result']['verification']['h'])
                points.append(v)
            for v in points:v['x_mm']=35*v['x']
            points.sort(key=lambda v:v['station']);all_points.extend(points);chains.append((source,points))
            inputs.append(dict(path=str(f),sha256=sha(f),finished=d.get('finished'),accepted=len(points)-1))
            failures.extend(dict(v,source=source,chain=spec['name'],result_file=str(f)) for v in d['failures'])
    assert all(abs(v['growth'])<=1e-7 for v in all_points),'non-neutral point in candidate export'
    base=ROOT/'outputs'/a.output
    fields=['source','station','x_mm','omega_r','beta_r','alpha_r','alpha_i','growth','g_beta','g_beta_beta','fd_h',
        'kind','solver_sha256','chain','result_file']
    with base.with_suffix('.csv').open('w',newline='') as out:
        w=csv.DictWriter(out,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(all_points)
    grid=json.loads((ROOT/'outputs/seed_catalog.json').read_text())['grids']['main_841']['x_nondim']
    coverage=[]
    for s in sorted(set(sources)):
        for i in coarse_stations(grid):
            points=[v for v in all_points if v['station']==i and v['source']==s]
            coverage.append(dict(source=s,station=i,x_mm=35*grid[i-1],
                status='neutral_candidate' if points else 'unknown_no_accepted_neutral_candidate',
                candidate_count=len(points),omega=points[-1]['omega_r'] if points else None,
                beta=points[-1]['beta_r'] if points else None))
    with base.with_name(base.name+'_coverage.csv').open('w',newline='') as out:
        w=csv.DictWriter(out,fieldnames=list(coverage[0]));w.writeheader();w.writerows(coverage)
    colors={'M1p':'#2878b5','M2p':'#e07b24','M1c':'#2b9162','M3near':'#9266b3'}
    fig,axs=plt.subplots(3,1,figsize=(11,9),constrained_layout=True)
    labeled=set()
    for source,points in chains:
        x=[v['x_mm'] for v in points];omega=[v['omega_r'] for v in points]
        kw=dict(color=colors.get(source),linewidth=1.35,marker='.',markersize=3.1)
        label=source if source not in labeled else None;labeled.add(source)
        axs[0].plot(x,omega,label=label,**kw);axs[1].plot(x,omega,**kw)
        axs[2].plot(x,[v['beta_r'] for v in points],**kw)
    axs[0].set(xlim=(-35,1100),ylabel=r'$\omega$ (dimensionless)',title='75-degree swept blunt plate: continued neutral candidates')
    axs[0].legend(ncol=4,loc='upper right',frameon=False)
    axs[1].set(xlim=(-35,160),ylabel=r'$\omega$ (dimensionless)',title='Nose and early flat-plate detail')
    axs[2].set(xlim=(-35,1100),ylabel=r'$\beta_*$ (dimensionless)',xlabel=r'Global wall $x$ (mm)')
    for ax in axs:
        ax.grid(True,alpha=.22);ax.axvline(0,color='.5',linewidth=.7)
    fig.suptitle('Individual branches only; unresolved portions remain open. Competition audit incomplete.',fontsize=10,color='#9d382b')
    fig.savefig(base.with_suffix('.png'),dpi=180);fig.savefig(base.with_suffix('.svg'));plt.close(fig)
    metadata=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),
        phases=a.phase,inputs=inputs,actual_records=len(all_points),
        unique_source_stations=len({(v['source'],v['station']) for v in all_points}),
        max_abs_growth=max(abs(v['growth']) for v in all_points),coverage=coverage,failures=failures,
        qualification='Actual single-branch candidates, not a certified finite-domain multi-branch envelope.',
        coordinate_model='Frozen wall-local tangent plane, zero explicit curvature; x remains global wall x.',
        frequency_units='nondimensional; no Hz conversion applied',cleanup='No source or calculation files deleted.')
    base.with_suffix('.json').write_text(json.dumps(metadata,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:metadata[k] for k in ['created','actual_records','unique_source_stations','max_abs_growth']}),flush=True)


if __name__=='__main__':main()
