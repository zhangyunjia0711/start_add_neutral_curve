#!/usr/bin/env python3
"""Bounded, actual local-beta-maximum paths below two native tail upper roots."""
import argparse, datetime, json, math, time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from full_driver import ROOT, sha, snapshot
from fd_envelope import stencil, refine
from local_peak_x_continue import predicted_point, execute, seed

def stamp():
    return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()

def bounded(x,lo,hi): return float(max(lo,min(hi,x)))

def actual(args,origin,w,b,label):
    old=args.label; args.label=old+'_'+label
    try: return predicted_point(args,841,origin,float(w),float(b))
    finally: args.label=old

def fd(args,state,w,b,label,h=.00125,tight=False):
    return stencil(841,state,float(w),float(b),args.solver,'wall_local',h=h,
                   tight=tight,label=args.label+'_'+label)

def correct_peak(args,origin,w,prediction,window,label):
    result=dict(status='unknown',omega=float(w),prediction=prediction,
                beta_window=window,trace=[],all_observed_peaks=[])
    state,transport=actual(args,origin,w,prediction,label+'_predict')
    result['parameter_transport_run']=transport
    b=state['beta_r']; h=.00125
    try:
        for iteration in range(12):
            if not window[0]+2*h<b<window[1]-2*h:
                raise RuntimeError('Peak iterate reaches actual adaptive identity window')
            cur=fd(args,state,w,b,label+'_peak'); result['trace'].append(cur)
            state=cur['center']; slope=cur['F'][1]; curvature=cur['J'][1][1]
            if curvature>=0: raise RuntimeError('Predicted state is not locally concave in beta')
            displacement=-slope/curvature
            cur['beta_newton_displacement']=displacement
            if abs(displacement)<2e-6:
                check=fd(args,state,w,b,label+'_h2',h=h/2,tight=True)
                result['verification']=check
                kappa=check['J'][1][1]
                db=-check['F'][1]/kappa if kappa else float('inf')
                result['beta_h2_correction']=db
                # Actual neighbours certify sampled strict local maximum, not a
                # full beta-domain search. Source labels/MAC alone are not used.
                same=[v for v in check['samples'] if abs(v['omega_r']-w)<1e-10]
                ordered=sorted(same,key=lambda x:x['beta_r'])
                shape=(len(ordered)==3 and ordered[1]['growth']>ordered[0]['growth']
                       and ordered[1]['growth']>ordered[2]['growth'])
                result['local_shape_samples']=ordered
                result['strict_sampled_local_maximum']=shape
                if kappa<0 and abs(db)<1e-5 and shape:
                    result.update(status='peak_verified',point=check['center'],
                        g_beta=check['F'][1],curvature=kappa,g_omega=check['J'][0][0],
                        ridge_beta_slope=-check['J'][1][0]/kappa,
                        beta_prediction_error=check['center']['beta_r']-prediction)
                    result['all_observed_peaks']=[check['center']]
                    return result
                if abs(db)<1e-5 and not shape:
                    # Refine to the peak before testing strict neighbours.
                    b+=db; state=check['center']; h/=2
                    if h<.00015625: raise RuntimeError('No strict peak shape at minimum finite stencil')
                    continue
            trust=bounded(math.sqrt(.02/abs(curvature)),.01,.25)
            delta=bounded(displacement,-trust,trust)
            cur['curvature_based_beta_trust']=trust
            trials=[]; accepted=False
            for line in range(5):
                nb=b+delta/(2**line)
                if not window[0]+2*h<nb<window[1]-2*h:
                    trials.append(dict(beta=nb,status='outside adaptive peak identity window'));continue
                trial=fd(args,state,w,nb,label+'_line',h=h)
                trials.append(trial)
                if trial['J'][1][1]<0 and abs(trial['F'][1])<abs(slope):
                    b=nb;state=trial['center'];accepted=True;break
            cur['line_search']=trials
            if not accepted: raise RuntimeError('No decreasing beta-peak step in bounded line search')
        raise RuntimeError('Fixed-frequency beta-peak corrector iteration limit')
    except Exception as error:
        result['reason']=repr(error);result['last_actual']=state
        return result

def recovery_scan(args,origin,w,prediction,window,label):
    center,name=actual(args,origin,w,prediction,label+'_center')
    rows={}; runs=[name]
    for side,lo,hi in [(-1,window[0],prediction),(1,prediction,window[1])]:
        n=min(65,max(17,int(math.ceil((hi-lo)/.02))+1))
        values,run_name=execute(args,dict(station=841,seed=seed(center),
            omega_min=w,omega_max=w,omega_points=1,beta_min=lo,beta_max=hi,
            beta_points=n,reverse_inner=side<0),label+'_scan'+str(side))
        runs.append(run_name)
        for value in values:rows[round(value['beta_r'],12)]=value
    ordered=sorted(rows.values(),key=lambda x:x['beta_r'])
    peaks=[v for j,v in enumerate(ordered) if 0<j<len(ordered)-1
           and v['growth']>ordered[j-1]['growth'] and v['growth']>ordered[j+1]['growth']]
    evidence=dict(runs=runs,samples=ordered,all_observed_peaks=peaks,
                  boundary_samples=[ordered[0],ordered[-1]],window=window)
    if not peaks:return None,evidence
    ranked=sorted(peaks,key=lambda x:abs(x['beta_r']-prediction))
    if len(ranked)>1 and abs(ranked[1]['beta_r']-prediction)<2*abs(ranked[0]['beta_r']-prediction)+.02:
        evidence['reason']='Multiple comparable nearby beta peaks; identity ambiguous'
        return None,evidence
    evidence['selection']='Unique nearest predicted peak, actual local shape still required'
    return ranked[0],evidence

def main():
    p=argparse.ArgumentParser();p.add_argument('--source',choices=['beta4','beta6'],required=True)
    p.add_argument('--label',required=True);p.add_argument('--resume',action='store_true');p.add_argument('--max-rows',type=int,default=100)
    args=p.parse_args();args.solver=ROOT/'outputs/bin/START_LST_wall_local_v10';args.solver_sha=sha(args.solver);args.eigen_jump=4.
    assert args.solver_sha=='634b45c5d30e076f5b86ccb136e0734c1905beedde108477dceb8a75adc3853c'
    assert 1<=args.max_rows<=100
    source=ROOT/('outputs/local_peak_tail_beta4_down_v10_hsmall_progress.json' if args.source=='beta4'
                 else 'outputs/local_peak_tail_beta6_v10_hsmall_progress.json')
    record=next(r for r in json.loads(source.read_text())['records'] if r['station']==841)
    upper=record['point'];path=ROOT/('outputs/'+args.label+'.json');before=snapshot();started=time.time()
    doc=dict(created=stamp(),source=args.source,station=841,x_mm=1100,upper_seed=upper,
        upper_provenance=dict(path=str(source),sha256=sha(source)),solver=str(args.solver),solver_sha256=args.solver_sha,
        policy=dict(frequency_step_max=.1,accepted_frequency_rows_max=args.max_rows,bounded_recoveries_max=1,
          beta_identity='Actual local maximum, previous actual peak/J prediction, adaptive local window; MAC is additional eigenvector continuity only',
          beta_position_tolerance=1e-5,h=.00125,h2=.000625,original_polynomial_residual=1e-7,minimum_overlap=.95),
        rows=[],failures=[],lower_roots=[],recovery_used=False,finished=False,
        qualification='Two finite peak paths; no completeness across beta or eigenmodes, and no interpolation through unknown intervals')
    if path.exists():
        if not args.resume:raise RuntimeError('Existing phase preserved; select --resume or new label')
        old=json.loads(path.read_text());assert old['source']==args.source and old['solver_sha256']==args.solver_sha
        if old['finished']:print('FINISHED ARCHIVE '+str(path));return
        doc=old
    def save():
        doc['updated']=stamp();doc['seconds_this_invocation']=time.time()-started
        tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(doc,indent=2,allow_nan=False)+'\n');tmp.replace(path)
    previous=upper
    if not doc['rows']:
        initial=fd(args,upper,upper['omega_r'],upper['beta_r'],'initial',h=.000625,tight=True)
        assert initial['J'][1][1]<0 and initial['J'][0][0]<0
        doc['initial_actual_verification']=initial; save()
    else:previous=doc['rows'][-1]['point']
    J=doc['rows'][-1]['verification']['J'] if doc['rows'] else doc['initial_actual_verification']['J']
    step=.05 if doc['recovery_used'] else .1
    while previous['omega_r']>1e-10 and len(doc['rows'])<args.max_rows:
        w=max(0.,previous['omega_r']-step)
        k=-J[1][0]/J[1][1];pred=previous['beta_r']+k*(w-previous['omega_r'])
        width=bounded(math.sqrt(.02/abs(J[1][1])),.025,.25)
        halfwidth=bounded(max(4*width,2*abs(pred-previous['beta_r'])+.05),.15,1.2)
        window=[min(previous['beta_r'],pred)-halfwidth,max(previous['beta_r'],pred)+halfwidth]
        label='row'+str(len(doc['rows'])+1)
        output=correct_peak(args,previous,w,pred,window,label)
        output.update(previous_actual_peak=previous,ridge_predictor_k=k,frequency_step=w-previous['omega_r'],
                      curvature_width_scale=width,adaptive_halfwidth=halfwidth)
        if output['status']!='peak_verified':
            doc['failures'].append(output);save()
            if doc['recovery_used']:
                doc['stop_reason']='Unknown peak after the one bounded recovery';break
            doc['recovery_used']=True;step=.05
            rw=max(0.,previous['omega_r']-step);rp=previous['beta_r']+k*(rw-previous['omega_r'])
            rwindow=[min(previous['beta_r'],rp)-halfwidth,max(previous['beta_r'],rp)+halfwidth]
            try:
                peak,evidence=recovery_scan(args,previous,rw,rp,rwindow,label+'_recovery')
                output['bounded_recovery']=evidence;save()
                if peak is None:
                    doc['stop_reason']='Unknown local peak identity in bounded recovery scan';break
                output=correct_peak(args,peak,rw,peak['beta_r'],rwindow,label+'_recovered')
                output.update(previous_actual_peak=previous,ridge_predictor_k=k,frequency_step=rw-previous['omega_r'],
                    recovery_evidence=evidence,adaptive_halfwidth=halfwidth)
                if output['status']!='peak_verified':
                    doc['failures'].append(output);doc['stop_reason']='Unknown after bounded recovery corrector';break
            except Exception as error:
                doc['failures'].append(dict(status='unknown',reason=repr(error),target_omega=rw))
                doc['stop_reason']='Bounded recovery solve failed';break
        point=output['point'];doc['rows'].append(output)
        if previous['growth']>1e-7 and point['growth']<-1e-7:
            frac=previous['growth']/(previous['growth']-point['growth'])
            nw=previous['omega_r']+frac*(point['omega_r']-previous['omega_r'])
            nb=previous['beta_r']+frac*(point['beta_r']-previous['beta_r'])
            lower=refine(841,previous,nw,nb,args.solver,'wall_local',label=args.label+'_lower',max_iterations=20,
                omega_bounds=(point['omega_r'],previous['omega_r']),beta_bounds=output['beta_window'],h=.00125)
            lower['positive_to_negative_bracket_descending_omega']=[previous,point]
            lower['validated_lower_boundary']=bool(lower['status']=='candidate_converged' and
                lower['verification']['J'][0][0]>0 and lower['verification']['J'][1][1]<0)
            doc['lower_roots'].append(lower)
        print('TAIL RIDGE '+json.dumps(dict(source=args.source,n=len(doc['rows']),point=point,
            curvature=output['curvature'],g_omega=output['g_omega'])),flush=True)
        previous=point;J=output['verification']['J'];save()
    if previous['omega_r']<=1e-10:doc['stop_reason']='Reached actual omega zero on this locally verified peak path'
    elif len(doc['rows'])>=args.max_rows:doc['stop_reason']='Accepted frequency-row budget exhausted; remainder unknown'
    doc.update(finished=True,last_actual_peak=previous,reached_zero=previous['omega_r']<=1e-10,
               accepted_frequency_rows=len(doc['rows']),original_inputs_unchanged=(snapshot()==before))
    save();print('TAIL DONE '+json.dumps({k:doc[k] for k in ['source','stop_reason','accepted_frequency_rows','reached_zero']}),flush=True)

if __name__=='__main__':main()
