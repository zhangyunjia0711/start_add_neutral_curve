#!/usr/bin/env python3
"""Finite normal derivative bracketing from an actual stopped ridge iterate."""
import argparse,json,time
from pathlib import Path
import numpy as np
import full_driver as driver
from local_peak_x_continue import stamp
from ridge_shear_corrector import Sampler

ROOT=driver.ROOT

def main():
    p=argparse.ArgumentParser();p.add_argument('--normal-report',type=Path,required=True)
    p.add_argument('--label',required=True);p.add_argument('--max-steps',type=int,default=15)
    p.add_argument('--path-step',type=float,default=2e-5);p.add_argument('--h-eta',type=float,default=6.25e-8)
    a=p.parse_args();prior=json.loads(a.normal_report.read_text());source=json.loads(Path(prior['source_report']).read_text())
    column=prior['column'];state=prior['last_actual'];window=prior['beta_window'];w=state['omega_r']
    target=ROOT/f'outputs/{a.label}.json'
    if target.exists():raise RuntimeError('Prior probe preserved')
    sampler=Sampler(a,source);started=time.monotonic()
    report=dict(created=stamp(),source_report=str(a.normal_report),source_sha256=driver.sha(a.normal_report),
        original_auxiliary_report=prior['source_report'],column=column,auxiliary_only=True,physical_candidates=[],
        omega_fixed=w,beta_window=window,initial_actual=state,steps=[],status='running',
        policy=dict(max_steps=a.max_steps,max_backtracks=4,h_eta=a.h_eta,path_step=a.path_step,
            minimum_MAC=.95,PEP_residual_gate=1e-7,normal_location_tolerance=4*a.h_eta),
        qualification='Finite actual normal peak bracket; no spectral completeness or physical endpoint claim')
    def save():
        report.update(updated=stamp(),seconds=time.monotonic()-started,run_names=sampler.runs)
        t=target.with_suffix('.tmp');t.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');t.replace(target)
    def normal(origin,beta,suffix):
        q,name=sampler.point(column,origin,w,beta,suffix+'_c');paths=[name];ends=[]
        for sign in [-1,1]:
            v,name=sampler.point(column,q,w,beta+sign*a.h_eta,suffix+f'_e{sign}');paths.append(name);ends.append(v)
        gb=(ends[1]['growth']-ends[0]['growth'])/(2*a.h_eta)
        gbb=(ends[1]['growth']-2*q['growth']+ends[0]['growth'])/(a.h_eta*a.h_eta)
        return dict(center=q,g_beta=gb,curvature=gbb,normal_displacement=-gb/gbb if gbb else None,paths=paths,ends=ends)
    driver.snapshot();save();left=None;right=None
    try:
        for j in range(a.max_steps):
            n=normal(state,state['beta_r'],f's{j}');report['steps'].append(n);state=n['center'];save()
            gb,gbb=n['g_beta'],n['curvature'];b=state['beta_r']
            if gb>0:left=n
            if gb<0:right=n
            if gbb<0 and abs(n['normal_displacement'])<4*a.h_eta:
                report.update(status='actual_normal_maximum_located',point=state,normal_result=n);save();return
            if left and right and left['center']['beta_r']<right['center']['beta_r']:
                blo,bhi=left['center']['beta_r'],right['center']['beta_r']
                proposal=(blo*right['g_beta']-bhi*left['g_beta'])/(right['g_beta']-left['g_beta'])
                if not blo+.1*(bhi-blo)<proposal<bhi-.1*(bhi-blo):proposal=.5*(blo+bhi)
                n['derivative_bracket']=[blo,bhi]
            else:
                growth_scale=np.sqrt(2*max(abs(state['growth']),1e-4)/max(abs(gbb),1e-30))
                increment=min(.005,max(1e-4,growth_scale,abs(n['normal_displacement'] or 0)))
                proposal=b+np.sign(gb)*increment
                proposal=float(np.clip(proposal,window[0]+4*a.h_eta,window[1]-4*a.h_eta))
                n.update(quadratic_growth_change_scale=float(growth_scale),unbracketed_increment=increment)
            if abs(proposal-b)<a.h_eta:
                report.update(status='unknown_reached_anchored_window_edge',last_actual=state);save();return
            n['proposal_beta']=float(proposal);n['trials']=[];moved=False
            for line in range(4):
                nb=b+(proposal-b)/2**line
                try:
                    trial=normal(state,float(nb),f's{j}_line{line}');trial['actual_delta_beta']=float(nb-b)
                    n['trials'].append(trial);save()
                    # A bracket sign change is retained even if growth does not
                    # improve: the final root must still have negative curvature.
                    sign_change=trial['g_beta']*gb<0
                    if sign_change or trial['center']['growth']>=state['growth']-1e-8:
                        if trial['g_beta']>0:left=trial
                        if trial['g_beta']<0:right=trial
                        state=trial['center'];moved=True;break
                except Exception as error:n['trials'].append(dict(beta=float(nb),reason=repr(error)));save()
            if not moved:raise RuntimeError('All four actual normal search steps rejected or decreased growth without derivative bracketing')
            report['last_actual']=state;save()
        report.update(status='unknown_fifteen_normal_steps_exhausted',last_actual=state);save()
    except Exception as error:report.update(status='stopped_unresolved',reason=repr(error),last_actual=state);save()
    finally:driver.snapshot();report['finished']=stamp();save()

if __name__=='__main__':main()
