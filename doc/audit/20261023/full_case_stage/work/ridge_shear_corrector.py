#!/usr/bin/env python3
"""Validate actual oblique parameter differences, then try one bounded next root.

All sample locations are solved in omega/beta. Auxiliary x is never physical.
The local change is beta=beta0+k*(omega-omega0)+eta; k is held during each FD.
"""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
import full_driver as driver
from local_peak_x_continue import stamp,seed,eigen

ROOT=driver.ROOT

class Sampler:
    def __init__(self,args,source):
        self.args=args;self.source=source;self.runs=[]
    def point(self,column,origin,w,b,suffix,tight=False):
        spec=dict(mode=self.source['template'],auxiliary_input_hashes=self.source['input_hashes'],
            station=column,solver=self.source['solver'],solver_sha256=self.source['solver_sha256'],
            coordinates='wall_local',seed=seed(origin),seed_shift=1e-4,task=33,
            omega_min=float(w),omega_max=float(w),omega_points=1,
            beta_min=float(b-2e-12),beta_max=float(b+2e-12),beta_points=3,
            anchor_row_seed=True,sample_only=True,adaptive=True,vector_tracking=True,
            minimum_overlap=.95,retry_eigenpairs=8,minimum_step=1e-10,step=self.args.path_step,
            jump=8,eigen_tol=1e-7,parameter_tol=1e-10,growth_tol=1e-9,timeout=90,
            options=['-pep_tol','1e-14' if tight else '1e-13','-pep_max_it','1000',
                '-pep_nev','1','-pep_ncv','16' if tight else '8','-st_pc_factor_mat_solver_type','superlu'])
        digest=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()[:12]
        spec['name']=self.args.label+'_'+suffix+'_'+digest
        result=driver.run(spec);self.runs.append(spec['name'])
        if result['exit'] not in (0,1) or result['petsc_errors'] or result['unused_options'] or result['coordinate_metadata_missing']:
            raise RuntimeError('Solver/configuration failure: '+spec['name'])
        rows=driver.table(ROOT/'runs'/spec['name']/'out/neutral_field.plt')
        if result['rejected'] or result['failed_field_records'] or len(rows)!=3 or any(x['status']!=0 for x in rows):
            raise RuntimeError('Actual parameter path/sample unresolved: '+spec['name'])
        point=min(rows,key=lambda x:abs(x['beta_r']-b))
        if abs(point['omega_r']-w)>2e-13 or abs(point['beta_r']-b)>2e-13:
            raise RuntimeError('Target coordinate not actually reached')
        return point,spec['name']
    def stencil(self,column,origin,w,b,k,hu,he,suffix,tight=False):
        values={};paths=[]
        center,name=self.point(column,origin,w,b,suffix+'_center',tight);paths.append(name);values[0,0]=center
        for i in [0,-1,1]:
            row=center
            if i:
                row,name=self.point(column,center,w+i*hu,b+k*i*hu,suffix+f'_u{i}',tight)
                paths.append(name);values[i,0]=row
            for j in [-1,1]:
                q,name=self.point(column,row,w+i*hu,b+k*i*hu+j*he,suffix+f'_u{i}_e{j}',tight)
                paths.append(name);values[i,j]=q
        g=lambda i,j:values[i,j]['growth']
        du=(g(1,0)-g(-1,0))/(2*hu);de=(g(0,1)-g(0,-1))/(2*he)
        cross=(g(1,1)-g(1,-1)-g(-1,1)+g(-1,-1))/(4*hu*he)
        curvature=(g(0,1)-2*g(0,0)+g(0,-1))/(he*he)
        J=np.array([[du,de],[cross,curvature]],dtype=float);F=np.array([g(0,0),de])
        delta=np.linalg.solve(J,-F)
        correction=[float(delta[0]),float(k*delta[0]+delta[1])]
        return dict(column=column,center=center,k=k,h_omega=hu,h_eta=he,
            transformation='omega=omega0+u; beta=beta0+k*u+eta',F=F.tolist(),J_shear=J.tolist(),
            J_original_inferred=[[float(du-k*de),float(de)],[float(cross-k*curvature),float(curvature)]],
            correction_shear=delta.tolist(),correction_omega_beta=correction,
            condition_shear=float(np.linalg.cond(J)),paths=paths,
            actual_samples=[dict(u=i*hu,eta=j*he,point=values[i,j]) for i in [-1,0,1] for j in [-1,0,1]])
    def gradient(self,column,origin,w,b,he,suffix):
        center,name=self.point(column,origin,w,b,suffix+'_c');names=[name];ends=[]
        for j in [-1,1]:
            q,name=self.point(column,center,w,b+j*he,suffix+f'_e{j}');ends.append(q);names.append(name)
        return dict(center=center,F=[center['growth'],(ends[1]['growth']-ends[0]['growth'])/(2*he)],paths=names)

def accepted_root(s):
    return bool(s['J_shear'][1][1]<0 and abs(s['F'][0])<1e-7
        and abs(s['correction_omega_beta'][0])<1e-4 and abs(s['correction_omega_beta'][1])<1e-5
        and abs(s['F'][1])<max(1e-6,abs(s['J_shear'][1][1])*1e-5))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-report',type=Path,required=True);p.add_argument('--column',type=int,required=True)
    p.add_argument('--label',required=True);p.add_argument('--h-omega',type=float,default=2e-6)
    p.add_argument('--h-eta',type=float,default=1.25e-7);p.add_argument('--path-step',type=float,default=5e-6)
    p.add_argument('--try-next',action='store_true');p.add_argument('--max-iterations',type=int,default=8)
    a=p.parse_args();source=json.loads(a.source_report.read_text())
    record=next(x for x in source['records'] if x['auxiliary_column']==a.column)
    point=record['point'];oldJ=np.array(record['result']['verification']['J']);k=float(-oldJ[1,0]/oldJ[1,1])
    target=ROOT/f'outputs/{a.label}.json'
    if target.exists():raise RuntimeError('Existing shear phase is preserved')
    assert driver.sha(source['solver'])==source['solver_sha256'];driver.snapshot();started=time.monotonic()
    sampler=Sampler(a,source)
    output=dict(created=stamp(),source_report=str(a.source_report),source_sha256=driver.sha(a.source_report),
        script_sha256=driver.sha(Path(__file__)),solver=source['solver'],solver_sha256=source['solver_sha256'],
        auxiliary_only=True,physical_candidates=[],actual_path_maximum_step=a.path_step,minimum_step=1e-10,
        original_MAC_threshold=.95,original_PEP_residual_gate=1e-7,
        status='running',validation={},next_attempt=None,qualification='Numerical local shear-coordinate audit and bounded auxiliary-root attempt; no physical endpoint or global spectrum claim')
    def save():
        output.update(updated=stamp(),seconds=time.monotonic()-started,run_names=sampler.runs)
        temp=target.with_suffix('.tmp');temp.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n');temp.replace(target)
    save()
    try:
        w,b=point['omega_r'],point['beta_r']
        coarse=sampler.stencil(a.column,point,w,b,k,a.h_omega,a.h_eta,'validation_h')
        output['validation']['h']=coarse;save()
        fine=sampler.stencil(a.column,coarse['center'],w,b,k,a.h_omega/2,a.h_eta/2,'validation_h2',True)
        output['validation']['h2']=fine
        rel=lambda x,y:float(abs(x-y)/max(abs(y),1e-30))
        changes=dict(ridge_growth_derivative_relative=rel(coarse['J_shear'][0][0],fine['J_shear'][0][0]),
            beta_curvature_relative=rel(coarse['J_shear'][1][1],fine['J_shear'][1][1]),
            center_alpha_difference=abs(eigen(coarse['center'])-eigen(fine['center'])),
            center_to_existing_alpha_difference=abs(eigen(fine['center'])-eigen(point)),
            remaining_mixed_derivative=fine['J_shear'][1][0],original_mixed_derivative=float(oldJ[1,0]))
        passed=accepted_root(coarse) and accepted_root(fine) and changes['ridge_growth_derivative_relative']<.05 and changes['beta_curvature_relative']<.05 and changes['center_to_existing_alpha_difference']<1e-7
        output['validation'].update(changes=changes,passed=bool(passed));save()
        if not passed:output['status']='stopped_known_root_validation_failed';save();return
        if not a.try_next:output['status']='known_root_validation_passed_only';save();return
        failure=next(x for x in source['failures'] if x['auxiliary_column']==a.column+1)
        state=failure['predicted_actual'];w,b=state['omega_r'],state['beta_r']
        beta_window=[point['beta_r']-.05,point['beta_r']+.05]
        next_report=dict(column=a.column+1,source_actual=failure['predicted_actual'],source_actual_report=str(a.source_report),
            beta_window=beta_window,trace=[],status='unknown')
        output['next_attempt']=next_report
        # Use the last actual root's ridge orientation; re-align only from newly solved derivatives.
        k=float(k-fine['J_shear'][1][0]/fine['J_shear'][1][1])
        for it in range(a.max_iterations):
            s=sampler.stencil(a.column+1,state,w,b,k,a.h_omega/2,a.h_eta/2,f'next_i{it}')
            next_report['trace'].append(s);save();state=s['center']
            if accepted_root(s):
                verify=sampler.stencil(a.column+1,state,w,b,k,a.h_omega/4,a.h_eta/4,f'next_i{it}_verify',True)
                next_report['verification']=verify;save()
                if accepted_root(verify):
                    next_report.update(status='auxiliary_candidate_converged',point=verify['center'])
                    output['status']='next_auxiliary_root_corrected_requires_native_endpoint';save();return
            J=np.array(s['J_shear']);F=np.array(s['F']);delta=np.linalg.solve(J,-F)
            scale=min(1.,2e-4/max(abs(delta[0]),1e-100),2e-4/max(abs(delta[1]),1e-100))
            weights=np.array([max(abs(J[0,0]),1.)*1e-4,max(abs(J[1,1]),1.)*1e-5])
            merit=float(np.linalg.norm(F/weights));s.update(unlimited_delta_shear=delta.tolist(),trust_scale=scale,merit=merit,line_search=[])
            moved=False
            for line in range(4):
                factor=scale/2**line;du,de=factor*delta;nw=w+du;nb=b+k*du+de
                if not(-80<nw<80 and beta_window[0]<nb<beta_window[1]):continue
                trial=sampler.gradient(a.column+1,state,float(nw),float(nb),a.h_eta/2,f'next_i{it}_line{line}')
                trial_merit=float(np.linalg.norm(np.array(trial['F'])/weights))
                s['line_search'].append(dict(factor=factor,merit=trial_merit,actual=trial))
                save()
                if trial_merit<merit*(1-1e-4*factor):
                    state=trial['center'];w,b=float(nw),float(nb);moved=True;break
            if not moved:raise RuntimeError('No decreasing admissible shear-coordinate step in bounded line search')
            if abs(J[1,1])>1e-20:
                revised=float(k-J[1,0]/J[1,1])
                if np.isfinite(revised) and abs(revised)<100:k=revised
        raise RuntimeError('Bounded shear corrector iteration limit')
    except Exception as error:
        output.update(status='stopped_unresolved',reason=repr(error));save()
    finally:
        driver.snapshot();output['finished']=stamp();save()

if __name__=='__main__':main()
