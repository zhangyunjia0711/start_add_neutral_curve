#!/usr/bin/env python3
"""Bounded actual normal-peak correction alternating with ridge frequency steps."""
import argparse,json,time
from pathlib import Path
import numpy as np
import full_driver as driver
from local_peak_x_continue import stamp
from ridge_shear_corrector import Sampler,accepted_root

ROOT=driver.ROOT

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--joint-report',type=Path,required=True);p.add_argument('--label',required=True)
    p.add_argument('--max-cycles',type=int,default=15);p.add_argument('--normal-iterations',type=int,default=8)
    p.add_argument('--h-omega',type=float,default=1e-6);p.add_argument('--h-eta',type=float,default=6.25e-8)
    p.add_argument('--path-step',type=float,default=5e-6);p.add_argument('--omega-trust',type=float,default=3e-4)
    a=p.parse_args();joint=json.loads(a.joint_report.read_text())
    assert joint['validation']['passed']
    source_path=Path(joint['source_report']);source=json.loads(source_path.read_text())
    next0=joint['next_attempt'];column=next0['column'];last=next0['trace'][-1]
    state=last['line_search'][-1]['actual']['center'] if last.get('line_search') else last['center']
    k=float(last['k']-last['J_shear'][1][0]/last['J_shear'][1][1]);window=next0['beta_window']
    output=dict(created=stamp(),joint_report=str(a.joint_report),joint_sha256=driver.sha(a.joint_report),
        source_report=str(source_path),source_sha256=driver.sha(source_path),script_sha256=driver.sha(Path(__file__)),
        source_script_sha256=driver.sha(ROOT/'work/ridge_shear_corrector.py'),column=column,auxiliary_only=True,
        physical_candidates=[],initial_actual=state,initial_k=k,beta_window=window,status='running',cycles=[],
        policy=dict(max_cycles=a.max_cycles,normal_iterations_per_cycle=a.normal_iterations,max_backtracks=4,
            h_omega=a.h_omega,h_eta=a.h_eta,omega_trust=a.omega_trust,minimum_MAC=.95,
            original_PEP_residual_gate=1e-7,path_step=a.path_step,normal_peak_location_tolerance=4*a.h_eta),
        qualification='Auxiliary local beta maximum plus ridge-frequency iteration; initialization only until native endpoint and h/2 proof')
    target=ROOT/f'outputs/{a.label}.json'
    if target.exists():raise RuntimeError('Existing alternating phase preserved')
    sampler=Sampler(a,source);started=time.monotonic()
    def save():
        output.update(updated=stamp(),seconds=time.monotonic()-started,run_names=sampler.runs)
        temp=target.with_suffix('.tmp');temp.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n');temp.replace(target)
    def normal(origin,w,b,suffix):
        center,name=sampler.point(column,origin,w,b,suffix+'_c');paths=[name];ends=[]
        for sign in [-1,1]:
            q,name=sampler.point(column,center,w,b+sign*a.h_eta,suffix+f'_e{sign}')
            paths.append(name);ends.append(q)
        gradient=(ends[1]['growth']-ends[0]['growth'])/(2*a.h_eta)
        curvature=(ends[1]['growth']-2*center['growth']+ends[0]['growth'])/(a.h_eta*a.h_eta)
        return dict(center=center,g_beta=gradient,curvature=curvature,paths=paths,
            actual_samples=ends,normal_Newton_displacement=-gradient/curvature if curvature else None)
    driver.snapshot();save()
    try:
        for cycle in range(a.max_cycles):
            c=dict(index=cycle,normal=[],frequency_trials=[]);output['cycles'].append(c);save()
            peak=False;bracket={}
            for ni in range(a.normal_iterations):
                n=normal(state,state['omega_r'],state['beta_r'],f'c{cycle}_n{ni}')
                c['normal'].append(n);state=n['center'];save()
                slope,curvature=n['g_beta'],n['curvature'];b=state['beta_r']
                if slope>0:bracket['left']=b
                elif slope<0:bracket['right']=b
                if curvature<0 and abs(n['normal_Newton_displacement'])<4*a.h_eta:
                    peak=True;c['normal_status']='actual_beta_maximum_located';break
                if curvature>=0:raise RuntimeError('Normal correction reached positive curvature before an actual maximum was located')
                # This scale is computed from actual local growth/curvature; it
                # is a numerical quadratic proposal, not a physical peak width.
                scale=max(8*a.h_eta,np.sqrt(2*max(abs(state['growth']),1e-4)/abs(curvature)))
                trust=min(2e-4,scale)
                delta=float(np.clip(n['normal_Newton_displacement'],-trust,trust))
                if 'left' in bracket and 'right' in bracket and bracket['left']<bracket['right']:
                    proposed=b+delta
                    if not bracket['left']<proposed<bracket['right']:
                        delta=.5*(bracket['left']+bracket['right'])-b
                n.update(quadratic_growth_change_scale=scale,normal_trust=trust,proposed_delta=delta,bracket=bracket.copy(),trials=[])
                moved=False
                for line in range(4):
                    db=delta/2**line;nb=b+db
                    if not window[0]<nb<window[1]:continue
                    trial=normal(state,state['omega_r'],nb,f'c{cycle}_n{ni}_line{line}')
                    trial['actual_delta_beta']=db;n['trials'].append(trial);save()
                    improves=(trial['center']['growth']>state['growth']-1e-9 and
                        (abs(trial['g_beta'])<abs(slope) or trial['g_beta']*slope<0))
                    if improves and trial['curvature']<0:
                        state=trial['center'];moved=True;break
                if not moved:raise RuntimeError('No admissible increasing-growth/decreasing-slope normal-peak step in four backtracks')
            if not peak:raise RuntimeError('Normal peak iteration limit before scalar frequency correction')
            s=sampler.stencil(column,state,state['omega_r'],state['beta_r'],k,a.h_omega,a.h_eta,f'c{cycle}_ridge')
            c['ridge_stencil']=s;state=s['center'];save()
            if accepted_root(s):
                verify=sampler.stencil(column,state,state['omega_r'],state['beta_r'],k,a.h_omega/2,a.h_eta/2,f'c{cycle}_verify',True)
                c['verification']=verify;save()
                if accepted_root(verify):
                    output.update(status='auxiliary_candidate_converged',point=verify['center'],result=verify);save();return
            if s['J_shear'][1][1]>=0:raise RuntimeError('Located normal maximum did not retain negative curvature in actual ridge stencil')
            k=float(k-s['J_shear'][1][0]/s['J_shear'][1][1]);c['updated_k']=k
            if not np.isfinite(k) or abs(k)>=100:raise RuntimeError('Nonfinite or out-of-range actual ridge orientation')
            # Recompute the actual tangent after reorientation, rather than
            # changing the coordinate label on the previous sampled values.
            tangent=sampler.stencil(column,state,state['omega_r'],state['beta_r'],k,a.h_omega,a.h_eta,f'c{cycle}_aligned')
            c['aligned_stencil']=tangent;save();gu=tangent['J_shear'][0][0];g=tangent['F'][0]
            if not np.isfinite(gu) or abs(gu)<1e-8:raise RuntimeError('Actual ridge growth derivative too small for bounded scalar correction')
            dw=float(np.clip(-g/gu,-a.omega_trust,a.omega_trust));c['scalar_frequency_proposal']=dict(g=g,g_u=gu,unlimited=-g/gu,limited=dw)
            moved=False
            for line in range(4):
                step=dw/2**line;nw=state['omega_r']+step;nb=state['beta_r']+k*step
                if not(-80<nw<80 and window[0]<nb<window[1]):continue
                point,name=sampler.point(column,state,nw,nb,f'c{cycle}_frequency{line}')
                trial=dict(actual=point,run=name,actual_delta_omega=step,actual_delta_beta=k*step)
                c['frequency_trials'].append(trial);save()
                if abs(point['growth'])<abs(g)*(1-1e-4/2**line):
                    state=point;moved=True;break
            if not moved:raise RuntimeError('No decreasing actual growth scalar-frequency step in four backtracks')
            output['last_actual']=state;save()
        raise RuntimeError('Fifteen-cycle bounded alternating correction limit')
    except Exception as error:output.update(status='stopped_unresolved',reason=repr(error),last_actual=state);save()
    finally:driver.snapshot();output['finished']=stamp();save()

if __name__=='__main__':main()
