"""Bounded i841 middle-frequency known-seed audit; at most six seeds/twelve beta arms."""
import concurrent.futures,csv,datetime,hashlib,json,math,pathlib,re,struct,sys,time
import numpy as np
AREA=pathlib.Path(__file__).resolve().parent;RUN=AREA.parent.parent
COMP=RUN/'outputs/upper_boundary_139_v10_20261023/compensation_20261023_2148'
sys.path.insert(0,str(COMP));import compensate as c
c.AREA=AREA;c.a.AREA=AREA;c.a.driver.ROOT=AREA;c.a.fd.ROOT=AREA
OPTS=list(c.VECOPTS);OPTS[OPTS.index('-pep_tol')+1]='1e-13'
def write(n,v):c.write(n,v)
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def serial(d):return {k:v for k,v in d.items() if k not in ['q','center_q']}
def main_plan():
    O=RUN/'outputs';files=[O/'tail841_all_boundaries_20261024.json',
        O/'upper_boundary_audit_v10_20261023/root_inventory_enriched.json',
        O/'local_peak_tail_beta6_v10_hsmall_progress.json',O/'axial_wall_local_v10_near_homotopy.json']
    docs=[json.loads(p.read_text()) for p in files];jobs=[]
    for label,source in [('low_lobe_upper','tail841_low_lobe_upper'),('high_lobe_lower','tail841_lobe_lower')]:
        z=next(x for x in docs[0]['results'] if x['job']['source']==source and 'beta4' in x['job']['name'])
        jobs.append(dict(label=label,point=z['result']['point'],source_file=str(files[0]),source_sha256=sha(files[0]),
            source_location=z['job']['name'],source_metadata=dict(job=z['job'],verification=z['result']['verification']['name'])))
    gap=[x['point']['omega_r'] for x in jobs];mid=sum(gap)/2
    for label,rid in [('M1high_upper','r00193'),('M2_upper_peak_a','r00084')]:
        root=next(x for x in docs[1]['all_roots'] if x['station']==841 and x['root_id']==rid)
        jobs.append(dict(label=label,point=root['point'],source_file=str(files[1]),source_sha256=sha(files[1]),
            source_location=rid,source_metadata=root))
    record=next(x for x in docs[2]['records'] if x.get('station')==841)
    jobs.append(dict(label='M2_upper_peak_b',point=record['point'],source_file=str(files[2]),source_sha256=sha(files[2]),
        source_location='records station841',source_metadata=record))
    record=next(x for x in docs[3]['records'] if x.get('iloc')==841)
    p=dict(record,iloc=841.,omega_i=0.,beta_i=0.,Re=79100.,outer=0.,inner=-7.,branch=0.,row=1)
    jobs.append(dict(label='nominal_zero_beta_minus7',point=p,source_file=str(files[3]),source_sha256=sha(files[3]),
        source_location='records iloc841',source_metadata=record))
    assert len(jobs)==6
    for j in jobs:j['target_omega']=mid
    plan=dict(created=c.a.stamp(),station=841,omega_gap_edges=sorted(gap),target_omega=mid,beta_domain=[-40,40],
        jobs=jobs,workers=2,max_beta_arms=12,db_target=1,continuation_step=.4,minimum_MAC=.95,
        residual_tol=1e-7,nev=1,retry_eigenpairs=8,station_predictor=False,
        center_and_arm_timeout_seconds=120,fresh_timeout_seconds=30,
        existing_reuse_review='Existing B40/B60 runs use highest-root+0.05 frequency, not this exact midpoint. Existing tail local ridge traces use other omega values and finite local beta windows. Source facts reused; no unverified near-frequency path substituted.',
        qualification='Finite beta sampled paths of six actual known seeds only; missing/failed stays unknown. No full-spectrum or all-beta stability claim.')
    write('audit_plan',plan);return plan
def center(job,origin,refq):
    w=job['target_omega']
    r=c.a.sample(841,origin,w,origin['beta_r']-1e-8,origin['beta_r'],1,'G841_'+job['label']+'_center',timeout=120,
        extra=dict(options=OPTS,seed_shift=1e-6))
    out=dict(run=r,name=r['name'],source_identity_verified=False)
    try:
        meta,q=c.vec(r,origin);out.update(initial_vector=meta,initial_MAC=c.mac(q,refq))
        out['source_identity_verified']=out['initial_MAC']>=1-1e-8
        p=c.endpoint(r,w,origin['beta_r']);meta,q=c.vec(r,p)
        out.update(point=p,center_vector=meta,center_q=q)
    except Exception as e:out['reason']=repr(e)
    out['passed']=r['complete'] and out['source_identity_verified'] and 'center_q' in out
    return out
def arm(job,side,origin,refq):
    end=-40 if side=='left' else 40
    r=c.a.sample(841,origin,job['target_omega'],origin['beta_r'],end,1,'G841_'+job['label']+'_'+side,timeout=120,
        extra=dict(options=OPTS,seed_shift=1e-6))
    out=dict(side=side,endpoint=end,run=r,name=r['name'],source_identity_verified=False)
    try:
        start=c.endpoint(r,job['target_omega'],origin['beta_r']);meta,q=c.vec(r,start)
        out.update(start_vector=meta,start_MAC=c.mac(q,refq),start_alpha_delta=abs(c.point_eigen(start)-c.point_eigen(origin)))
        out['source_identity_verified']=out['start_MAC']>=1-1e-8 and out['start_alpha_delta']<1e-8
        if r['complete']:
            p=c.endpoint(r,job['target_omega'],end);meta,q=c.vec(r,p);out.update(endpoint_actual=p,endpoint_vector=meta)
    except Exception as e:out['reason']=repr(e)
    out['complete']=r['complete'] and out['source_identity_verified'] and 'endpoint_actual' in out
    positives=[p for p in r['points'] if c.a.valid(p) and abs(p['omega_r']-job['target_omega'])<1e-10 and p['growth']>1e-8]
    out['positive_points']=positives
    out['state']='positive_actual_sample' if positives else 'sampled_nonpositive_known_path' if out['complete'] else 'unknown'
    return out
def execute(job):
    started=time.time();out=dict(label=job['label'],job=job,started=c.a.stamp(),state='unknown',arms=[])
    try:
        fresh,info,q=c.fresh(841,job['point'],'G841_'+job['label']+'_fresh')
        p=c.endpoint(fresh,job['point']['omega_r'],job['point']['beta_r'])
        out['fresh']=dict(run=fresh,vector=info,actual_point=p,original_alpha_delta=abs(c.point_eigen(p)-c.point_eigen(job['point'])))
        seedvec=AREA/'runs'/fresh['name']/'out/eigenvectors.petsc';assert sha(seedvec)==info['sha256']
        cen=center(job,p,q);out['center']=serial(cen)
        if not cen['passed']:raise RuntimeError('Center transport or raw-vector source join incomplete')
        positives=[p for p in cen['run']['points'] if c.a.valid(p) and abs(p['omega_r']-job['target_omega'])<1e-10 and p['growth']>1e-8]
        if positives:print('GAP_ACTUAL_POSITIVE '+json.dumps(dict(label=job['label'],stage='center',point=max(positives,key=lambda x:x['growth']))),flush=True)
        for side in ['left','right']:
            a=arm(job,side,cen['point'],cen['center_q']);out['arms'].append(a)
            write('results/'+job['label'],out)
            if a['positive_points']:print('GAP_ACTUAL_POSITIVE '+json.dumps(dict(label=job['label'],stage=side,point=max(a['positive_points'],key=lambda x:x['growth']))),flush=True)
        out['positive_points']=positives+[p for a in out['arms'] for p in a['positive_points']]
        out['state']='positive_actual_sample' if out['positive_points'] else 'sampled_nonpositive_known_path' if all(a['complete'] for a in out['arms']) else 'unknown'
    except Exception as e:out['reason']=repr(e)
    out['seconds']=round(time.time()-started,3);out['finished']=True;write('results/'+job['label'],out)
    print('GAP_SOURCE_DONE '+json.dumps(dict(label=job['label'],state=out['state'],seconds=out['seconds'],reason=out.get('reason'))),flush=True)
    return out
def main():
    AREA.mkdir(exist_ok=True);(AREA/'results').mkdir(exist_ok=True)
    before=c.a.driver.snapshot();assert sha(c.a.SOLVER)==c.a.SHA
    existing=AREA/'audit_plan.json';plan=json.loads(existing.read_text()) if existing.exists() else main_plan()
    results=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        pending={}
        for job in plan['jobs']:
            p=AREA/'results'/(job['label']+'.json')
            if p.exists() and json.loads(p.read_text()).get('finished'):results.append(json.loads(p.read_text()))
            else:pending[pool.submit(execute,job)]=job
        for future in concurrent.futures.as_completed(pending):
            job=pending[future]
            try:r=future.result()
            except Exception as e:r=dict(label=job['label'],job=job,state='unknown',reason=repr(e),finished=True)
            results.append(r)
            write('audit_progress',dict(updated=c.a.stamp(),completed=len(results),planned=6,finished=len(results)==6,results=results))
    after=c.a.driver.snapshot();assert before==after and sha(c.a.SOLVER)==c.a.SHA
    write('audit_integrity',dict(created=c.a.stamp(),before=before,after=after,unchanged=True,solver_sha256=c.a.SHA))
    print('GAP_ALL_DONE '+json.dumps([(x['label'],x['state']) for x in results]),flush=True)
if __name__=='__main__':main()

