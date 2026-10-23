#!/usr/bin/env python3
"""Finite negative-beta path, followed by one bounded auxiliary/native transfer."""
import argparse,csv,hashlib,json,math,sys,time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import full_driver as driver
from local_peak_x_continue import stamp,seed,eigen
from ridge_shear_corrector import Sampler

ROOT=driver.ROOT
OPTIONS=['-pep_tol','1e-13','-pep_max_it','1000','-pep_nev','1','-pep_ncv','8',
         '-st_pc_factor_mat_solver_type','superlu']

def main():
    p=argparse.ArgumentParser();p.add_argument('--label',default='free_beta_M1p_aux29_to7_20261023')
    p.add_argument('--db',type=float,default=.02);p.add_argument('--beta-end',type=float,default=-7)
    p.add_argument('--transport',action='store_true');a=p.parse_args()
    source_path=ROOT/'outputs/normal_ridge_bracket_M1p29_fixed_20261023.json'
    aux_path=ROOT/'outputs/aux_neutral_M1p_17_16_s64_20261023.json'
    source=json.loads(source_path.read_text());aux=json.loads(aux_path.read_text());origin=source['last_actual']
    column=source['column'];b0=origin['beta_r'];w=origin['omega_r'];end=a.beta_end
    assert end<b0 and .0<a.db<=.05;count=math.ceil((b0-end)/a.db)+1
    target=ROOT/f'outputs/{a.label}.json'
    if target.exists():raise RuntimeError('Existing finite exploration is preserved')
    report=dict(created=stamp(),source=str(source_path),source_sha256=driver.sha(source_path),
        auxiliary_source=str(aux_path),auxiliary_sha256=driver.sha(aux_path),script_sha256=driver.sha(Path(__file__)),
        fixed_omega=w,start_actual=origin,finite_beta_domain=[end,b0],maximum_db=a.db,
        auxiliary_only=True,physical_candidates=[],status='running',
        qualification='Finite real-beta single-path exploration; unknowns, growing edges and other unsearched spectral branches remain unexhausted')
    started=time.monotonic()
    def save():
        report.update(updated=stamp(),seconds=time.monotonic()-started)
        tmp=target.with_suffix('.tmp');tmp.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');tmp.replace(target)
    def run(spec,suffix):
        spec.update(mode=aux['template'],auxiliary_input_hashes=aux['input_hashes'],
            solver=aux['solver'],solver_sha256=aux['solver_sha256'],coordinates='wall_local',seed_shift=1e-4,
            task=33,anchor_row_seed=True,sample_only=True,adaptive=True,vector_tracking=True,
            minimum_overlap=.95,retry_eigenpairs=8,station_predictor=False,minimum_step=1e-6,
            step=.002,jump=8,eigen_tol=1e-7,options=OPTIONS,timeout=180)
        digest=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()[:12]
        spec['name']=a.label+'_'+suffix+'_'+digest
        summary=driver.run(spec);rows=driver.table(ROOT/'runs'/spec['name']/'out/neutral_field.plt')
        fatal=bool(summary['exit'] not in (0,1) or summary['petsc_errors'] or summary['unused_options'] or summary['coordinate_metadata_missing'])
        return summary,rows,fatal
    driver.snapshot();save()
    try:
        summary,rows,fatal=run(dict(station=column,seed=seed(origin),omega_min=w,omega_max=w,omega_points=1,
            beta_min=end,beta_max=b0,beta_points=count,reverse_inner=True),'scan')
        good=[v for v in rows if v['status']==0 and not fatal]
        ordered=sorted(good,key=lambda q:q['beta_r'])
        peaks=[dict(point=q,left=ordered[j-1],right=ordered[j+1]) for j,q in enumerate(ordered)
            if 0<j<len(ordered)-1 and q['growth']>ordered[j-1]['growth'] and q['growth']>ordered[j+1]['growth']]
        complete=bool(not fatal and summary['rejected']==0 and summary['failed_field_records']==0 and len(rows)==count)
        edges=[ordered[0],ordered[-1]] if ordered else []
        report['scan']=dict(summary=summary,complete=complete,requested_points=count,accepted_points=len(good),
            actual_samples=rows,observed_internal_peaks=peaks,edge_points=edges,
            negative_boundary_still_growing=bool(len(ordered)>1 and abs(ordered[0]['beta_r']-end)<1e-10 and ordered[0]['growth']>ordered[1]['growth']))
        with (ROOT/f'outputs/{a.label}_scan.csv').open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=driver.COLUMNS+['row']);writer.writeheader();writer.writerows(rows)
        save()
        if not complete:raise RuntimeError('Finite beta scan incomplete; cannot propagate across its failed samples')
        report['peak_refinements']=[]
        sampler=Sampler(SimpleNamespace(label=a.label+'_peak',path_step=.002),aux)
        for j,peak in enumerate(peaks):
            left,right=peak['left']['beta_r'],peak['right']['beta_r'];state=peak['point'];trace=[]
            # Actual golden-section beta search, with the observed peak's true
            # eigenpair as its path anchor. All local maxima are kept separately.
            for iteration in range(24):
                if right-left<2e-6:break
                c=right-(right-left)*.6180339887498949;d=left+(right-left)*.6180339887498949
                pc,nc=sampler.point(column,state,w,c,f'p{j}_{iteration}_c')
                pd,nd=sampler.point(column,state,w,d,f'p{j}_{iteration}_d')
                trace.append(dict(left=left,right=right,c=pc,d=pd,runs=[nc,nd]))
                if pc['growth']>pd['growth']:right=d;state=pc
                else:left=c;state=pd
            h=min(1e-4,max((right-left),2e-6))
            stencils=[]
            for hs,tight in [(h,False),(h/2,True)]:
                center,n=sampler.point(column,state,w,state['beta_r'],f'p{j}_verify_c{hs:.2g}',tight)
                lo,nl=sampler.point(column,center,w,state['beta_r']-hs,f'p{j}_verify_l{hs:.2g}',tight)
                hi,nh=sampler.point(column,center,w,state['beta_r']+hs,f'p{j}_verify_r{hs:.2g}',tight)
                gb=(hi['growth']-lo['growth'])/(2*hs);gbb=(hi['growth']-2*center['growth']+lo['growth'])/(hs*hs)
                stencils.append(dict(h=hs,point=center,g_beta=gb,curvature=gbb,
                    location_correction=-gb/gbb if gbb else None,runs=[n,nl,nh]))
            report['peak_refinements'].append(dict(coarse=peak,trace=trace,verification=stencils,
                beta_peak_qualified=all(x['curvature']<0 and abs(x['location_correction'])<1e-5 for x in stencils),
                qualification='Actual local beta peak at fixed omega, not a neutral root'))
            save()
        chosen=max(good,key=lambda q:q['growth']);report['strong_positive_actual']=chosen;save()
        if chosen['growth']<=0:raise RuntimeError('No sampled positive-growth seed')
        if not a.transport:report['status']='finite_beta_scan_complete';save();return
        b=chosen['beta_r']
        summary,rows,fatal=run(dict(station=65,seed_station=column,seed=seed(chosen),
            outer_parameter='iloc',inner_parameter='beta',outer_min=column,outer_max=65,outer_points=66-column,
            inner_min=b,inner_max=b+1e-10,inner_points=2),'profile_transfer')
        prefix=[]
        for col in range(column,66):
            records=[v for v in rows if int(round(v['iloc']))==col]
            actual=next((v for v in records if abs(v['beta_r']-b)<1e-11 and v['status']==0),None)
            if fatal or len(records)!=2 or any(v['status']!=0 for v in records) or actual is None:break
            prefix.append(dict(auxiliary_column=col,theta=(col-1)/64,point=actual))
        report['profile_transfer']=dict(summary=summary,all_samples=rows,accepted_contiguous_prefix=prefix,
            completed_original_target=bool(len(prefix)==66-column and prefix[-1]['auxiliary_column']==65),
            qualification='Only a fixed omega/beta eigenpair transported through auxiliary numerical profiles')
        save()
        if not report['profile_transfer']['completed_original_target']:raise RuntimeError('Strong-growth auxiliary transfer stopped before native target')
        sys.path.insert(0,str(ROOT/'work/flat_neutral_normal_validation_20261023'))
        from validate_flat import solve,vectors
        here=ROOT/'work'/a.label;last=prefix[-1]['point'];solver=Path(aux['solver']);sha=aux['solver_sha256']
        actual_aux=solve('aux65_fixed_positive',Path(aux['input_hashes'].keys().__iter__().__next__()).parent,
            400,65,last,here=here,binary=solver,sha=sha,domain_in=65)
        native=solve('native16_fixed_positive',driver.CASE/'M1/in',400,16,last,here=here,binary=solver,sha=sha)
        av=vectors(Path(actual_aux['native_vector_file']),400)[actual_aux['raw_index']]
        nv=vectors(Path(native['native_vector_file']),400)[native['raw_index']]
        mac=float(abs(np.vdot(av,nv))**2/(np.vdot(av,av).real*np.vdot(nv,nv).real))
        delta=abs(complex(actual_aux['alpha_r']-native['alpha_r'],actual_aux['alpha_i']-native['alpha_i']))
        passed=bool(actual_aux['passed'] and native['passed'] and delta<1e-8 and mac>=1-1e-10 and native['growth']>0)
        report['native_endpoint_check']=dict(auxiliary=actual_aux,native=native,alpha_difference=delta,raw_primitive_MAC=mac,passed=passed)
        if not passed:raise RuntimeError('Fresh native positive-growth endpoint check failed')
        native_point=dict(omega_r=w,beta_r=b,alpha_r=native['alpha_r'],alpha_i=native['alpha_i'],growth=native['growth'],iloc=16)
        report.update(status='native_positive_growth_seed_verified',native_positive_growth_seed=native_point,
            native_seed_qualification='Fresh original-H5 fixed-parameter positive growth only; not a neutral root or global beta maximum')
        save()
    except Exception as error:report.update(status='stopped_unresolved',reason=repr(error));save()
    finally:driver.snapshot();report['finished']=stamp();save()

if __name__=='__main__':main()
