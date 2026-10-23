#!/usr/bin/env python3
"""Continue stationary neutral roots through verified auxiliary local profiles.

No auxiliary column is a physical x station. The final point is only an
initialization until independently checked on the unmodified native H5.
"""
import argparse,configparser,hashlib,json
from pathlib import Path
from types import SimpleNamespace
import full_driver as driver
import fd_envelope as fd
import local_peak_x_continue as lp

ROOT=driver.ROOT

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-manifest',type=Path,required=True)
    p.add_argument('--source',required=True);p.add_argument('--label',required=True)
    p.add_argument('--first-delta',type=float,nargs=2,default=[0.,0.],metavar=('OMEGA','BETA'),
        help='Parameter increment for the first auxiliary step only, based on actual native history')
    p.add_argument('--h',type=float,default=1e-6);p.add_argument('--minimum-stencil-h',type=float,default=1e-7)
    p.add_argument('--beta-halfwidth',type=float,default=.1);p.add_argument('--peak-reseed-radius',type=float,default=.005)
    p.add_argument('--eigen-jump',type=float,default=8);p.add_argument('--max-iterations',type=int,default=32)
    p.add_argument('--start-only',action='store_true');p.add_argument('--resume',action='store_true')
    a=p.parse_args();manifest=json.loads(a.input_manifest.read_text())
    assert manifest['status']=='prepared_and_start_endpoint_verified'
    assert manifest['start_endpoint_check']['passed']
    solver=Path(manifest['solver']);assert driver.sha(solver)==manifest['solver_sha256']
    inputs=Path(manifest['input_dir']);count=manifest['auxiliary_columns']
    hashes={str(inputs/n):driver.sha(inputs/n) for n in ['grid_vertical2d_N.h5','flow_vertical2d_N.h5']}
    template=ROOT/'work/aux_neutral_templates'/a.label;template.mkdir(parents=True,exist_ok=True)
    link=template/'in'
    if not link.exists():link.symlink_to(inputs,target_is_directory=True)
    assert link.resolve()==inputs.resolve()
    config=template/'config.ini';cfg=configparser.ConfigParser();cfg.read(manifest['config_template'])
    cfg['domain'].update({'in':str(count),'jn':'400','istart':'1','iend':str(count),'iloc':'1'})
    if config.exists():
        old=configparser.ConfigParser();old.read(config);assert dict(old['domain'])==dict(cfg['domain'])
    else:
        with config.open('w') as f:cfg.write(f)
    original_run=driver.run
    def auxiliary_run(spec):
        spec.update(mode=str(template),auxiliary_input_hashes=hashes,station_predictor=False)
        digest=hashlib.sha256(json.dumps(hashes,sort_keys=True).encode()).hexdigest()[:8]
        spec['name']+='_aux'+digest
        return original_run(spec)
    # Per-process adapters only: core driver, original configs and binaries stay unchanged.
    fd.run=auxiliary_run;lp.run=auxiliary_run
    args=SimpleNamespace(solver=solver,solver_sha=manifest['solver_sha256'],label=a.label,
        eigen_jump=a.eigen_jump,h=a.h,minimum_stencil_h=a.minimum_stencil_h,max_iterations=a.max_iterations,
        strict_gradient=False,beta_location_tol=1e-5,omega_min=-80,omega_max=80,
        omega_trust=.05,beta_trust=.05,peak_reseed_radius=a.peak_reseed_radius)
    output=dict(created=lp.stamp(),source=a.source,manifest_path=str(a.input_manifest),
        manifest_sha256=driver.sha(a.input_manifest),input_hashes=hashes,template=str(template),
        solver=str(solver),solver_sha256=args.solver_sha,auxiliary_columns=count,
        corrector_source_sha256=driver.sha(Path(__file__)),
        local_peak_source_sha256=driver.sha(Path(lp.__file__)),fd_source_sha256=driver.sha(Path(fd.__file__)),
        native_from=manifest['native_from'],native_to=manifest['native_to'],
        first_step_parameter_increment=a.first_delta,records=[],failures=[],finished=False,
        auxiliary_only=True,physical_candidates=[],
        policy=dict(initial_h=a.h,minimum_stencil_h=a.minimum_stencil_h,beta_halfwidth=a.beta_halfwidth,
            peak_reseed_radius=a.peak_reseed_radius,eigen_jump=a.eigen_jump,max_iterations=a.max_iterations,
            minimum_MAC=.95,maximum_PEP_residual=1e-7,beta_location_tolerance=1e-5,
            omega_location_tolerance=1e-4,negative_curvature_required=True),
        qualification='Auxiliary numerical profiles only. Final omega/beta/alpha is initialization, not a physical candidate until fresh native endpoint proof and native h/2 root verification.')
    target=ROOT/f'outputs/{a.label}.json'
    if target.exists():
        if not a.resume:raise RuntimeError('Existing auxiliary phase preserved; use --resume or a new label')
        prior=json.loads(target.read_text())
        for k in ['source','manifest_sha256','input_hashes','solver_sha256','first_step_parameter_increment','policy']:
            assert prior[k]==output[k],'Changed auxiliary resume policy: '+k
        if prior['finished']:print('Existing phase finished: '+str(target),flush=True);return
        output=prior
    def save():
        output['updated']=lp.stamp();tmp=target.with_suffix('.tmp')
        tmp.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n');tmp.replace(target)
    driver.snapshot();save()
    history=[r for r in output['records'] if r['status']=='candidate_converged']
    for column in range(1,count+1):
        if any(r['auxiliary_column']==column for r in history):continue
        record=dict(auxiliary_column=column,theta=(column-1)/(count-1),source=a.source,physical_station=False)
        try:
            if column==1:
                old=manifest['source_anchor_point'];wp,bp=old['omega_r'],old['beta_r']
                actual,name=lp.predicted_point(args,1,old,wp,bp)
                record.update(start_run=name,start_alpha_difference=abs(lp.eigen(actual)-lp.eigen(old)))
                if record['start_alpha_difference']>1e-8:raise RuntimeError('Auxiliary start eigenvalue not equivalent to native anchor')
            else:
                old=history[-1]['point']
                actual,name=lp.transport(args,column-1,column,old);record['transport_run']=name
                if len(history)>=2:
                    previous=history[-2]['point']
                    wp=old['omega_r']+(old['omega_r']-previous['omega_r'])
                    bp=old['beta_r']+(old['beta_r']-previous['beta_r'])
                else:wp,bp=old['omega_r']+a.first_delta[0],old['beta_r']+a.first_delta[1]
            window=[old['beta_r']-a.beta_halfwidth,old['beta_r']+a.beta_halfwidth]
            record.update(prediction=dict(omega=wp,beta=bp),beta_window=window,transported_actual=actual)
            if not window[0]+a.h<bp<window[1]-a.h:raise RuntimeError('Auxiliary prediction outside anchored local beta window')
            predicted,name=lp.predicted_point(args,column,actual,wp,bp)
            record.update(predicted_actual=predicted,prediction_run=name)
            result=lp.correct(args,column,predicted,wp,bp,window,a.label+f'_c{column:03d}')
            record['result']=result
            if result['status']!='candidate_converged':raise RuntimeError(result.get('reason','Auxiliary root corrector unresolved'))
            record.update(status='candidate_converged',point=result['point'])
            output['records'].append(record);history.append(record);save()
            print('AUX NEUTRAL '+json.dumps(dict(column=column,theta=record['theta'],point=result['point'],h2=result['verification']['h'])),flush=True)
        except Exception as error:
            record.update(status='unknown',reason=repr(error));output['failures'].append(record)
            output.update(status='stopped_unresolved_auxiliary_root',unresolved_auxiliary_columns=list(range(column,count+1)))
            break
        if a.start_only:output['status']='start_only_checked';break
    else:
        output.update(status='auxiliary_endpoint_root_corrected',endpoint_initialization=dict(
            source=a.source,auxiliary_column=count,native_target=manifest['native_to'],point=history[-1]['point'],
            result=history[-1]['result'],physical_candidate=False))
    output['finished']=True;save();driver.snapshot()

if __name__=='__main__':main()
