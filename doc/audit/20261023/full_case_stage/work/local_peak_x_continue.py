#!/usr/bin/env python3
"""Continue a local beta-stationary neutral root, preserving its neighbourhood.

Use physical-x secants for omega/beta prediction, true iloc/MAC transport, and
a corrector constrained to a fixed beta window around the preceding actual
root. No corrector can walk out to a remote competing peak. Failed intervals
are subdivided at original integer stations; an unresolved adjacent step stops.
"""
import argparse,datetime,hashlib,json,math
from pathlib import Path
import numpy as np
from full_driver import ROOT,run,table,sha,snapshot
from fd_envelope import stencil

OPTIONS=['-pep_tol','1e-12','-pep_max_it','1000','-pep_nev','1','-pep_ncv','8',
         '-st_pc_factor_mat_solver_type','superlu']

def stamp():return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
def seed(p):return {k:float(p[k]) for k in ['omega_r','beta_r','alpha_r','alpha_i']}
def eigen(p):return complex(p['alpha_r'],p['alpha_i'])

def execute(args,spec,suffix):
    spec.update(mode='M1',solver=str(args.solver),solver_sha256=args.solver_sha,
        coordinates='wall_local',seed_shift=1e-4,task=33,anchor_row_seed=True,
        sample_only=True,adaptive=True,vector_tracking=True,minimum_overlap=.95,
        retry_eigenpairs=8,minimum_step=1e-6,step=.08,jump=args.eigen_jump,options=OPTIONS,timeout=180)
    digest=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()[:12]
    spec['name']=args.label+'_'+suffix+'_'+digest
    summary=run(spec)
    if summary['exit'] not in [0,1] or summary['petsc_errors'] or summary['unused_options'] or summary['coordinate_metadata_missing']:
        raise RuntimeError('Runtime/configuration error: '+spec['name'])
    if summary['failed_field_records'] or summary['rejected']:
        raise RuntimeError('Unresolved requested transport sample: '+spec['name'])
    values=table(ROOT/'runs'/spec['name']/'out/neutral_field.plt')
    if not values or any(v['status']!=0 for v in values):raise RuntimeError('Incomplete transport: '+spec['name'])
    return values,spec['name']

def transport(args,origin,target,point):
    if origin==target:return dict(point),None
    beta=point['beta_r']
    values,name=execute(args,dict(station=target,seed_station=origin,seed=seed(point),
        outer_parameter='iloc',inner_parameter='beta',
        outer_min=min(origin,target),outer_max=max(origin,target),outer_points=2,
        inner_min=beta,inner_max=beta+1e-8,inner_points=2,
        reverse_outer=target<origin,station_predictor=True),f'x{origin}_{target}')
    actual=[v for v in values if abs(v['iloc']-target)<1e-8 and abs(v['beta_r']-beta)<1e-10]
    if not actual:raise RuntimeError('No actual target station: '+name)
    return actual[0],name

def predicted_point(args,station,actual,omega,beta):
    values,name=execute(args,dict(station=station,seed=seed(actual),
        omega_min=omega,omega_max=omega,omega_points=1,
        beta_min=beta-1e-9,beta_max=beta+1e-9,beta_points=3),f'predict_i{station}')
    return min(values,key=lambda v:abs(v['beta_r']-beta)),name

def nearby_peak(args,station,actual,window,suffix):
    """Retain only a unique observed maximum within an explicitly local window."""
    w,b=actual['omega_r'],actual['beta_r']
    left=max(window[0]+args.h,b-args.peak_reseed_radius)
    right=min(window[1]-args.h,b+args.peak_reseed_radius)
    evidence=dict(bounds=[left,right],levels=[])
    for level in range(getattr(args,'peak_refinement_levels',2)+1):
        rows={};runs=[];b=actual['beta_r']
        for side,lo,hi in [(-1,left,b),(1,b,right)]:
            values,name=execute(args,dict(station=station,seed=seed(actual),
                omega_min=w,omega_max=w,omega_points=1,beta_min=lo,beta_max=hi,
                beta_points=33,reverse_inner=side<0),suffix+f'_L{level}_arm{side}')
            runs.append(name)
            for v in values:rows[round(v['beta_r'],13)]=v
        ordered=sorted(rows.values(),key=lambda v:v['beta_r'])
        indices=[j for j,v in enumerate(ordered) if 0<j<len(ordered)-1
                 and v['growth']>ordered[j-1]['growth'] and v['growth']>ordered[j+1]['growth']]
        peaks=[ordered[j] for j in indices]
        evidence['levels'].append(dict(bounds=[left,right],runs=runs,samples=ordered,peaks=peaks))
        if len(peaks)!=1:return None,evidence
        j=indices[0];actual=ordered[j];left=ordered[j-1]['beta_r'];right=ordered[j+1]['beta_r']
    evidence['peaks']=peaks
    return actual,evidence

def correct(args,station,actual,omega,beta,window,label):
    output=dict(status='unknown',station=station,beta_window=window,trace=[],
        qualification='Locally constrained beta-peak root; no full beta/spectrum completeness claim')
    state=dict(actual);w,b=float(omega),float(beta);h=args.h
    try:
        for iteration in range(args.max_iterations):
            if not window[0]+h<b<window[1]-h:raise RuntimeError('Stencil leaves anchored beta-peak window')
            current=stencil(station,state,w,b,args.solver,'wall_local',h=h,label=label)
            output['trace'].append(current);state=current['center']
            F=np.array(current['F']);J=np.array(current['J'])
            if not np.all(np.isfinite(J)) or np.linalg.cond(J)>1e10:
                raise RuntimeError('Ill-conditioned local stationary-root Jacobian')
            if J[1,1]>=0:
                if args.peak_reseed_radius>0:
                    peak,evidence=nearby_peak(args,station,state,window,label+f'_reseed{iteration}')
                    current['local_peak_reseed']=evidence
                    if peak is not None:
                        if abs(peak['beta_r']-b)<1e-11 and abs(peak['growth']-state['growth'])<1e-9:
                            raise RuntimeError('Nearby sampled peak does not resolve the positive-curvature stencil')
                        state=peak;w,b=peak['omega_r'],peak['beta_r'];continue
                raise RuntimeError('No unique nearby beta maximum recovered from positive-curvature predictor')
            slope_tol=1e-6 if args.strict_gradient else max(1e-6,abs(J[1,1])*args.beta_location_tol)
            current['beta_gradient_tolerance']=slope_tol
            if abs(F[0])<=1e-8 and abs(F[1])<=slope_tol:
                verify=stencil(station,state,w,b,args.solver,'wall_local',h=h/2,tight=True,label=label+'_verify')
                output['verification']=verify
                correction=np.linalg.solve(np.array(verify['J']),-np.array(verify['F']))
                verify_slope_tol=1e-6 if args.strict_gradient else max(1e-6,abs(verify['J'][1][1])*args.beta_location_tol)
                if verify['J'][1][1]<0 and abs(verify['F'][0])<=1e-7 and abs(verify['F'][1])<=verify_slope_tol and abs(correction[0])<1e-4 and abs(correction[1])<args.beta_location_tol:
                    output.update(status='candidate_converged',point=verify['center'],
                        derivative=verify['F'][1],curvature=verify['J'][1][1],
                        verification_correction=correction.tolist(),beta_gradient_tolerance=verify_slope_tol,
                        beta_location_tolerance=args.beta_location_tol)
                    return output
                h/=2;state=verify['center'];continue
            delta=np.linalg.solve(J,-F)
            scale=min(1.,args.omega_trust/max(abs(delta[0]),1e-100),args.beta_trust/max(abs(delta[1]),1e-100))
            current['unlimited_newton_delta']=delta.tolist();current['trust_scale']=scale
            merit=float(np.linalg.norm(F));trials=[];accepted=False
            for line in range(7):
                factor=scale*2**(-line);nw,nb=np.array([w,b])+factor*delta
                if not (args.omega_min<nw<args.omega_max and window[0]+h<nb<window[1]-h):
                    trials.append(dict(factor=factor,omega=float(nw),beta=float(nb),status='outside local window'))
                    continue
                trial=stencil(station,state,float(nw),float(nb),args.solver,'wall_local',h=h,label=label+'_line')
                new_merit=float(np.linalg.norm(trial['F']))
                trials.append(dict(factor=factor,stencil=trial,merit=new_merit))
                if trial['J'][1][1]<0 and new_merit<merit*(1-1e-4*factor):
                    w,b=float(nw),float(nb);state=trial['center'];accepted=True;break
            current['line_search']=trials
            if not accepted:
                if h/2>=args.minimum_stencil_h:
                    current['stencil_refinement_reason']='No decreasing admissible step with this finite-difference Jacobian'
                    h/=2;continue
                raise RuntimeError('No decreasing admissible local-peak Newton step at minimum stencil size')
        raise RuntimeError('Local corrector iteration limit')
    except Exception as error:
        output['reason']=repr(error);return output

def load_history(path,source):
    doc=json.loads(path.read_text());points={}
    if 'results' in doc:
        for record in doc['results']:
            if record['job']['source']==source and record['result']['status']=='candidate_converged':
                points[record['job']['station']]=record['result']['point']
    else:
        for record in doc.get('records',[]):
            if record.get('status')=='candidate_converged':points[record['station']]=record['point']
    return points

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--solver',type=Path,required=True);p.add_argument('--solver-sha',required=True)
    p.add_argument('--history-file',type=Path,default=ROOT/'outputs/tail_beta_roots_20261023.json')
    p.add_argument('--source',default='tail_beta6');p.add_argument('--history-stations',nargs='+',type=int,default=[783,788])
    p.add_argument('--stations',nargs='+',type=int,default=[793,803,813,841])
    p.add_argument('--validate-stations',nargs='*',type=int,default=[793])
    p.add_argument('--beta-halfwidth',type=float,default=.5);p.add_argument('--beta-trust',type=float,default=.15)
    p.add_argument('--omega-trust',type=float,default=.25);p.add_argument('--omega-min',type=float,default=0)
    p.add_argument('--omega-max',type=float,default=20);p.add_argument('--h',type=float,default=.0025)
    p.add_argument('--max-iterations',type=int,default=32)
    p.add_argument('--eigen-jump',type=float,default=4.,help='Absolute eigenvalue jump gate; MAC/residual requirements remain unchanged')
    p.add_argument('--peak-reseed-radius',type=float,default=0.,help='If positive-curvature predictor occurs, inspect this local beta radius and accept only a unique observed maximum')
    p.add_argument('--minimum-stencil-h',type=float,default=1e-5)
    p.add_argument('--beta-location-tol',type=float,default=1e-5,
        help='Engineering local beta-coordinate tolerance, checked by an actual h/2 stencil Newton correction')
    p.add_argument('--strict-gradient',action='store_true',help='Retain the stricter raw |g_beta|<=1e-6 tail-validation policy')
    p.add_argument('--resume',action='store_true',help='Resume this exact archived phase; otherwise an existing progress file is never overwritten')
    p.add_argument('--label',default='local_peak_tail_beta6_v10');args=p.parse_args()
    args.solver=args.solver.resolve();assert sha(args.solver)==args.solver_sha
    assert len(args.history_stations)>=2 and args.beta_halfwidth>2*args.h
    snapshot();known=load_history(args.history_file,args.source)
    history=[dict(station=i,point=known[i],provenance='existing independently verified root') for i in args.history_stations]
    grid=json.loads((ROOT/'outputs/seed_catalog.json').read_text())['grids']['main_841']['x_nondim']
    output=dict(created=stamp(),solver=str(args.solver),solver_sha256=args.solver_sha,
        source=args.source,history_path=str(args.history_file),history_sha256=sha(args.history_file),
        initial_history=history.copy(),requested=args.stations,records=[],failures=[],finished=False,
        policy=dict(beta_halfwidth=args.beta_halfwidth,beta_trust=args.beta_trust,omega_trust=args.omega_trust,
                    physical_x_secant_prediction=True,intermediate_native_stations=True,eigen_jump=args.eigen_jump,
                    peak_reseed_radius=args.peak_reseed_radius),
        numerical_tolerances=dict(growth=1e-7,beta_gradient=1e-6,beta_verification_correction=1e-5,
            omega_verification_correction=1e-4,initial_stencil_h=args.h,max_corrector_iterations=args.max_iterations,
            beta_location_tol=args.beta_location_tol,strict_gradient=args.strict_gradient,minimum_stencil_h=args.minimum_stencil_h,
            gradient_policy='raw 1e-6' if args.strict_gradient else 'max(1e-6, abs(beta_curvature)*beta_location_tol)'),
        qualification='Anchored local beta-stationary candidates; source label alone never certifies persistent peak or global envelope')
    target=ROOT/f'outputs/{args.label}_progress.json'
    if target.exists():
        if not args.resume:raise RuntimeError('Existing phase is preserved; use --resume or a new --label')
        archived=json.loads(target.read_text())
        for key in ['solver_sha256','source','history_sha256','requested','initial_history','policy']:
            assert archived[key]==output[key],'Resume configuration differs: '+key
        if archived.get('finished'):
            print('Existing phase finished; archived evidence preserved: '+str(target),flush=True);return
        assert archived.get('numerical_tolerances')==output['numerical_tolerances'],'Resume numerical policy differs; use a new label'
        output=archived
        history.extend(dict(station=r['station'],point=r['point'],provenance='resumed actual continuation')
                       for r in output['records'] if r['status']=='candidate_converged')
    def save():
        output['updated']=stamp();tmp=target.with_suffix('.tmp')
        tmp.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n');tmp.replace(target)
    def attempt(station,depth=0):
        previous=history[-1];earlier=history[-2];old=previous['point']
        x0,x1,xt=grid[earlier['station']-1],grid[previous['station']-1],grid[station-1]
        fraction=(xt-x1)/(x1-x0)
        wp=old['omega_r']+fraction*(old['omega_r']-earlier['point']['omega_r'])
        bp=old['beta_r']+fraction*(old['beta_r']-earlier['point']['beta_r'])
        window=[old['beta_r']-args.beta_halfwidth,old['beta_r']+args.beta_halfwidth]
        record=dict(station=station,from_station=previous['station'],prediction=dict(omega=wp,beta=bp,
            physical_x_fraction=fraction),beta_window=window,depth=depth,requested=station in args.stations)
        try:
            if not window[0]+args.h<bp<window[1]-args.h:
                raise RuntimeError('Physical-x parameter prediction leaves anchored peak window; split station interval')
            actual,transport_name=transport(args,previous['station'],station,old)
            record.update(transport_run=transport_name,transported_actual=actual)
            predicted,prediction_name=predicted_point(args,station,actual,wp,bp)
            record.update(prediction_transport_run=prediction_name,predicted_actual=predicted)
            result=correct(args,station,predicted,wp,bp,window,args.label+f'_i{station}')
            record['result']=result
            if result['status']!='candidate_converged':raise RuntimeError(result.get('reason','corrector unresolved'))
            point=result['point'];record.update(status='candidate_converged',point=point,
                actual_beta_change=point['beta_r']-old['beta_r'],prediction_error=dict(
                    omega=point['omega_r']-wp,beta=point['beta_r']-bp))
            if station in known:
                ref=known[station];record['independent_known_root_difference']=dict(
                    omega=abs(point['omega_r']-ref['omega_r']),beta=abs(point['beta_r']-ref['beta_r']),
                    alpha=abs(eigen(point)-eigen(ref)))
            output['records'].append(record);history.append(dict(station=station,point=point,provenance='new constrained continuation'))
            save();print('LOCAL PEAK '+json.dumps({k:record[k] for k in ['station','status','point','actual_beta_change']}),flush=True)
            return True
        except Exception as error:
            record.update(status='unknown',reason=repr(error));output['failures'].append(record);save()
            start=previous['station']
            if abs(station-start)<=1:return False
            middle=(station+start)//2
            if not attempt(middle,depth+1):return False
            return attempt(station,depth+1)
    for j,station in enumerate(args.stations):
        if any(r['station']==station and r['status']=='candidate_converged' for r in output['records']):continue
        if station in known and station not in args.validate_stations:
            history.append(dict(station=station,point=known[station],provenance='reused independent root'))
            output['records'].append(dict(station=station,status='candidate_converged',point=known[station],reused_existing=True));save();continue
        if not attempt(station):output['unresolved_requested']=args.stations[j:];break
    output['finished']=True;save();snapshot()

if __name__=='__main__':main()
