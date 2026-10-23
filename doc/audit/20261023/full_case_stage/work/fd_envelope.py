"""Refine a smooth interior beta-envelope: g=0 and dg/dbeta=0.

This follows one locally verified eigenmode. Global beta/family competition and
domain-boundary checks are separate and mandatory before accepting an envelope.
"""
import argparse,hashlib,json,math,pathlib
import numpy as np
from full_driver import ROOT,run,table,sha,snapshot

def stencil(station,seed,omega,beta,solver,coordinates='wall_local',h=.02,tight=False,label='fd'):
    spec=dict(station=station,mode='M1',seed={k:float(seed[k]) for k in ['omega_r','beta_r','alpha_r','alpha_i']},
        solver=str(solver),solver_sha256=sha(solver),coordinates=coordinates,task=33,seed_shift=1e-4,
        omega_min=omega-h,omega_max=omega+h,omega_points=3,beta_min=beta-h,beta_max=beta+h,beta_points=3,
        step=.2,adaptive=True,vector_tracking=True,minimum_overlap=.95,minimum_step=1e-6,
        sample_only=True,anchor_row_seed=True,retry_eigenpairs=8,parameter_tol=1e-8,growth_tol=1e-9,eigen_tol=1e-7,
        options=['-pep_tol','1e-13' if tight else '1e-12','-pep_max_it','1000',
                 '-pep_nev','1','-pep_ncv','16' if tight else '8','-st_pc_factor_mat_solver_type','superlu'],timeout=300)
    digest=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()[:12]
    spec['name']=f'{label}_i{station:03d}_{digest}'
    result=run(spec)
    if result['petsc_errors'] or result['unused_options'] or result['coordinate_metadata_missing'] or result['exit']==124:
        raise RuntimeError('Solver/coordinate error: '+spec['name'])
    values={}
    for v in table(ROOT/'runs'/spec['name']/'out/neutral_field.plt'):
        a=(v['omega_r']-omega)/h;b=(v['beta_r']-beta)/h
        ia,ib=round(a),round(b)
        if abs(a-ia)<1e-7 and abs(b-ib)<1e-7 and abs(ia)<=1 and abs(ib)<=1:
            if v['status']!=0:raise RuntimeError('Failed stencil sample: '+spec['name'])
            values[ia,ib]=v
    if len(values)!=9:raise RuntimeError(f'Incomplete stencil {spec["name"]}: {len(values)}/9')
    g=lambda i,j:values[i,j]['growth']
    gradient=np.array([(g(1,0)-g(-1,0))/(2*h),(g(0,1)-g(0,-1))/(2*h)])
    jac=np.array([[gradient[0],gradient[1]],[(g(1,1)-g(1,-1)-g(-1,1)+g(-1,-1))/(4*h*h),(g(0,1)-2*g(0,0)+g(0,-1))/(h*h)]])
    return dict(name=spec['name'],center=values[0,0],h=h,F=[g(0,0),gradient[1]],J=jac.tolist(),samples=list(values.values()))

def refine(station,seed,omega,beta,solver,coordinates='wall_local',label='fd',max_iterations=15,
           omega_bounds=(0,80),beta_bounds=(-40,40),h=.02,allow_intermediate_curvature=False,
           beta_position_tol=1e-5,omega_position_tol=1e-4):
    trace=[];state=dict(seed);w,b=float(omega),float(beta)
    output=dict(station=station,coordinates=coordinates,solver=str(solver),solver_sha256=sha(solver),
        status='unresolved',trace=trace,finite_domain=dict(omega=omega_bounds,beta=beta_bounds),
        position_tolerances=dict(beta=beta_position_tol,omega=omega_position_tol),
        note='Single-mode interior candidate; not yet an audited multi-family envelope boundary.')
    try:
        for iteration in range(max_iterations):
            point=stencil(station,state,w,b,solver,coordinates,h=h,label=label)
            trace.append(point);F=np.array(point['F']);J=np.array(point['J'])
            state=point['center']
            if J[1,1]>=0 and not allow_intermediate_curvature:
                raise RuntimeError('Stationary point is not a strict beta maximum')
            # Curvature varies by orders of magnitude between the nose and
            # plate. Test the implied beta displacement as well as growth;
            # an absolute derivative alone would demand inconsistent positions.
            slope_tol=max(1e-6,abs(J[1,1])*beta_position_tol)
            if abs(F[0])<=1e-8 and abs(F[1])<=slope_tol:
                if J[1,1]>=0:raise RuntimeError('Converged stationary point is not a strict beta maximum')
                verify=stencil(station,state,w,b,solver,coordinates,h=h/2,tight=True,label=label+'_verify')
                output['verification']=verify
                correction=np.linalg.solve(np.array(verify['J']),-np.array(verify['F']))
                verify_slope_tol=max(1e-6,abs(verify['J'][1][1])*beta_position_tol)
                if verify['J'][1][1]<0 and abs(verify['F'][0])<=1e-7 and abs(verify['F'][1])<=verify_slope_tol and abs(correction[0])<omega_position_tol and abs(correction[1])<beta_position_tol:
                    output.update(status='candidate_converged',point=verify['center'],derivative=verify['F'][1],curvature=verify['J'][1][1],verification_correction=correction.tolist())
                    break
                # A finite-difference bias cannot be corrected by repeatedly
                # solving the same coarse stencil. Continue on the finer one.
                if h>1e-4:
                    h/=2;point=verify;state=verify['center'];F=np.array(verify['F']);J=np.array(verify['J'])
            delta=np.linalg.solve(J,-F)
            scale=min(1.,2/max(abs(delta[0]),1e-100),2/max(abs(delta[1]),1e-100))
            delta*=scale
            new_w,new_b=np.array([w,b])+delta
            if not (omega_bounds[0]<=new_w<=omega_bounds[1] and beta_bounds[0]<new_b<beta_bounds[1]):
                raise RuntimeError('Newton step leaves declared domain; needs separate endpoint/onset treatment')
            # Retain the actual previous parameter state as the next path origin;
            # starting at an unverified predicted state could silently change mode.
            samples={(round((v['omega_r']-w)/h),round((v['beta_r']-b)/h)):v for v in point['samples']}
            eigen=lambda a,c:complex(samples[a,c]['alpha_r'],samples[a,c]['alpha_i'])
            predicted=eigen(0,0)+(eigen(1,0)-eigen(-1,0))/(2*h)*delta[0]+(eigen(0,1)-eigen(0,-1))/(2*h)*delta[1]
            point['next_eigen_prediction']=[predicted.real,predicted.imag]
            state=point['center']
            w,b=float(new_w),float(new_b)
        else:output['reason']='maximum iterations reached'
    except Exception as e:output['reason']=repr(e)
    return output

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--validate',action='store_true');p.add_argument('--label',default='fd_validation2');a=p.parse_args()
    if a.validate:
        snapshot()
        binary=ROOT/'outputs/bin/START_LST_vectors_v4'
        seed=dict(omega_r=20,beta_r=8.29004379818586,alpha_r=42.26896603780138,alpha_i=-.1130588763967948)
        result=refine(321,seed,20,8.29004379818586,binary,'global',label=a.label,omega_bounds=(16,28),beta_bounds=(-20,20))
        (ROOT/f'outputs/{a.label}.json').write_text(json.dumps(result,indent=2)+'\n')
        print('FD VALIDATION '+json.dumps({k:v for k,v in result.items() if k not in ['trace','verification']}),flush=True)
        snapshot()
