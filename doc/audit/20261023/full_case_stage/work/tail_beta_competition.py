#!/usr/bin/env python3
"""Compare both actually continued tail beta peaks at common x and omega.

The two original neutral candidates are independently transported with MAC,
then scanned over the same finite beta interval. Every local sampled maximum
and unresolved sample is retained; no smooth-x line is treated as mode identity.
"""
import concurrent.futures,csv,datetime,json,math
from pathlib import Path
from types import SimpleNamespace
from full_driver import ROOT,sha,snapshot
from trace_envelope_x import transfer
from ridge_path_survey import point_run

SOLVER=ROOT/'outputs/bin/START_LST_wall_local_v7'
SHA='7aa6016570f406587156f229107c0e5e4f76c8a85ba89dda267039ff1d688d79'
TARGET_OMEGA=4.84
SOURCE=ROOT/'outputs/neutral_x_cont_v7a_M2p_downstream_progress.json'
TARGET=ROOT/'outputs/tail_beta_competition_20261023.json'

def stamp():return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()

def curve(station,label,record):
    result=dict(station=station,origin_label=label,origin_station=record['station'],
                origin_point=record['point'],omega=TARGET_OMEGA,beta_window=[3,8],status='unknown')
    args=SimpleNamespace(station=station,solver=SOLVER,solver_sha=SHA,timeout=90,
                         label=f'tail_compete_20261023_i{station}_{label}')
    try:
        actual,transport_run=transfer(record['station'],station,record['point'],SOLVER,args.label)
        result.update(transport_run=transport_run,transported_actual=actual)
        beta=actual['beta_r']
        center_run=point_run(args,actual,TARGET_OMEGA,beta-1e-9,beta+1e-9,3)
        result['frequency_transport']=center_run
        if not center_run['complete']:raise RuntimeError('Common-frequency center unresolved')
        origin=min(center_run['accepted_samples'],key=lambda v:abs(v['beta_r']-beta))
        result['actual_common_frequency_origin']=origin
        arms=[]
        for lo,hi,reverse in [(3,beta,True),(beta,8,False)]:
            arms.append(point_run(args,origin,TARGET_OMEGA,lo,hi,
                max(2,math.ceil((hi-lo)/.1-1e-12)+1),reverse))
        result['arms']=arms
        unique={round(v['beta_r'],10):v for arm in arms for v in arm['accepted_samples']}
        values=sorted(unique.values(),key=lambda v:v['beta_r'])
        result['samples']=values;result['local_sampled_peaks']=[]
        for j in range(1,len(values)-1):
            a,b,c=values[j-1:j+2]
            if b['growth']>=a['growth'] and b['growth']>=c['growth']:
                result['local_sampled_peaks'].append(b)
        result['observed_winner']=max(values,key=lambda v:v['growth']) if values else None
        result['status']='finite_window_sampled' if all(arm['complete'] for arm in arms) else 'partial_window_unknown'
    except Exception as error:result['reason']=repr(error)
    return result

def main():
    snapshot();assert sha(SOLVER)==SHA
    source=json.loads(SOURCE.read_text())
    records={v['station']:v for v in source['records'] if v.get('status')=='candidate_converged'}
    seeds={'old_beta6':records[783],'new_beta4':records[793]}
    output=dict(created=stamp(),solver=str(SOLVER),solver_sha256=SHA,
        source_path=str(SOURCE),source_sha256=sha(SOURCE),seeds=seeds,
        omega=TARGET_OMEGA,beta_window=[3,8],stations=[783,788,793],curves=[],same_beta_comparisons=[],
        qualification='Finite two-origin common-parameter competition audit; no full-spectrum or global envelope guarantee')
    def save():
        output['updated']=stamp();tmp=TARGET.with_suffix('.tmp')
        tmp.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n');tmp.replace(TARGET)
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        pending={pool.submit(curve,i,label,record):(i,label) for i in [783,788,793] for label,record in seeds.items()}
        for future in concurrent.futures.as_completed(pending):
            result=future.result();output['curves'].append(result);save()
            print('TAIL CURVE '+json.dumps(dict(station=result['station'],origin=result['origin_label'],
                status=result['status'],peaks=[(v['beta_r'],v['growth']) for v in result.get('local_sampled_peaks',[])])),flush=True)
    # Exact common-beta repeats distinguish two actual eigenvalues from mere
    # differing beta grids; lambda differences are measured, not interpolated.
    for station in [783,788,793]:
        group=[v for v in output['curves'] if v['station']==station and v.get('actual_common_frequency_origin')]
        for beta in [seeds['old_beta6']['point']['beta_r'],seeds['new_beta4']['point']['beta_r']]:
            probes=[]
            for item in group:
                args=SimpleNamespace(station=station,solver=SOLVER,solver_sha=SHA,timeout=90,
                    label=f'tail_compete_20261023_i{station}_{item["origin_label"]}_samebeta')
                run=point_run(args,item['actual_common_frequency_origin'],TARGET_OMEGA,beta-1e-9,beta+1e-9,3)
                point=min(run['accepted_samples'],key=lambda v:abs(v['beta_r']-beta)) if run['complete'] else None
                probes.append(dict(origin_label=item['origin_label'],run=run,point=point))
            comparison=dict(station=station,omega=TARGET_OMEGA,beta=beta,probes=probes)
            good=[p['point'] for p in probes if p['point']]
            if len(good)==2:
                comparison['delta_alpha']=abs(complex(good[0]['alpha_r']-good[1]['alpha_r'],good[0]['alpha_i']-good[1]['alpha_i']))
                comparison['delta_growth']=abs(good[0]['growth']-good[1]['growth'])
            output['same_beta_comparisons'].append(comparison);save()
    columns=['station','origin_label','beta_r','omega_r','alpha_r','alpha_i','growth','source_run']
    with TARGET.with_suffix('.csv').open('w',newline='') as out:
        writer=csv.DictWriter(out,fieldnames=columns);writer.writeheader()
        for item in output['curves']:
            for point in item.get('samples',[]):
                writer.writerow(dict([(k,point.get(k)) for k in columns],station=item['station'],origin_label=item['origin_label']))
    output['finished']=True;save();snapshot()

if __name__=='__main__':main()
