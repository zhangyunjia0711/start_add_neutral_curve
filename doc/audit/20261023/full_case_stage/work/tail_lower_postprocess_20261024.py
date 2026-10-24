#!/usr/bin/env python3
"""Refine both sign-change directions and export actual native tail boundary seeds."""
import json,sys,time
from pathlib import Path
from full_driver import ROOT,CASE,sha,snapshot
from local_peak_x_continue import stamp
from fd_envelope import refine

def main():
    snapshot();solver=ROOT/'outputs/bin/START_LST_wall_local_v10'
    target=ROOT/'outputs/tail841_all_boundaries_20261024.json'
    if target.exists():raise RuntimeError('Existing boundary export preserved')
    report=dict(created=stamp(),station=841,results=[],path_audits=[],
        qualification='Neutral roots of finite observed local-beta maxima, not a complete all-beta/all-mode stable-region proof')
    for source in ['beta4','beta6']:
        path=ROOT/('outputs/tail_lower_'+source+'_841_fine_20261024.json')
        doc=json.loads(path.read_text());assert doc['finished']
        points=[doc['upper_seed']]+[r['point'] for r in doc['rows']]
        crossings=[]
        for index,(a,b) in enumerate(zip(points[:-1],points[1:])):
            if a['growth']*b['growth']>=0 or min(abs(a['growth']),abs(b['growth']))<1e-7:continue
            direction='lower' if a['growth']>0 else 'low_lobe_upper'
            if direction=='lower':
                q=next((q for q in doc['lower_roots'] if q.get('validated_lower_boundary')),None)
                if q is None:raise RuntimeError('Expected automatic lower-root evidence missing')
            else:
                fraction=a['growth']/(a['growth']-b['growth'])
                w=a['omega_r']+fraction*(b['omega_r']-a['omega_r'])
                beta=a['beta_r']+fraction*(b['beta_r']-a['beta_r'])
                row=doc['rows'][index]
                q=refine(841,a,w,beta,solver,'wall_local',label='tail841_'+source+'_'+direction+'_20261024',
                    max_iterations=20,omega_bounds=(b['omega_r'],a['omega_r']),beta_bounds=row['beta_window'],h=.00125)
            q=dict(q);q['descending_omega_sign_bracket']=[a,b]
            expected_sign=1 if direction=='lower' else -1
            q['validated_directional_boundary']=bool(q['status']=='candidate_converged' and
                q['verification']['J'][1][1]<0 and q['verification']['J'][0][0]*expected_sign>0)
            name='tail841_lobe_lower' if direction=='lower' else 'tail841_low_lobe_upper'
            item=dict(job=dict(station=841,source=name,name=name+'_'+source+'_i841'),result=q,
                      origin_peak_path=source,origin_path=str(path),origin_sha256=sha(path),direction=direction)
            report['results'].append(item);crossings.append(dict(direction=direction,
                actual_bracket=[a,b],status=q['status'],validated=q['validated_directional_boundary']))
            target.write_text(json.dumps(report,indent=2)+'\n')
        report['path_audits'].append(dict(source=source,path=str(path),sha256=sha(path),
            accepted_frequency_rows=len(doc['rows']),crossings=crossings,reached_zero=doc['reached_zero'],
            last_actual_peak=doc['last_actual_peak'],stop_reason=doc['stop_reason']))
    # Fresh original-H5 solves preserve native primitive vectors for each unique boundary.
    sys.path.insert(0,str(ROOT/'work/flat_neutral_normal_validation_20261023'))
    from validate_flat import solve
    unique={}
    for item in report['results']:
        q=item['result']
        if not q['validated_directional_boundary']:continue
        point=q['point'];key=item['direction']
        if key not in unique:
            native=solve('tail841_'+key+'_native_20261024',CASE/'M1/in',400,841,point,
                here=ROOT/'work/tail841_boundary_native_20261024',binary=solver,sha=sha(solver))
            difference=abs(complex(native['alpha_r']-point['alpha_r'],native['alpha_i']-point['alpha_i']))
            assert native['passed'] and abs(native['growth'])<1e-7 and difference<1e-7
            unique[key]=dict(native=native,alpha_difference=difference)
        item['actual_native_fresh']=unique[key]['native'];item['fresh_alpha_difference']=unique[key]['alpha_difference']
    report['finished']=True;report['original_inputs_unchanged']=bool(snapshot())
    target.write_text(json.dumps(report,indent=2)+'\n')
    print('BOUNDARIES '+json.dumps([dict(source=x['job']['source'],origin=x['origin_peak_path'],
        status=x['result']['status'],point=x['result'].get('point'),validated=x['result']['validated_directional_boundary'])
        for x in report['results']]),flush=True)

if __name__=='__main__':main()
