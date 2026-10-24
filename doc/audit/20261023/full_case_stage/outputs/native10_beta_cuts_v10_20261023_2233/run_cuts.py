"""Native10 beta-minimum seed: two bounded frequency cuts, never relabel as max."""
import collections,concurrent.futures,datetime,hashlib,json,math,pathlib,sys
import numpy as np
AREA=pathlib.Path(__file__).resolve().parent;RUN=AREA.parent.parent
COMP=RUN/'outputs/upper_boundary_139_v10_20261023/compensation_20261023_2148'
sys.path.insert(0,str(COMP));import compensate as c
c.AREA=AREA;c.a.AREA=AREA;c.a.driver.ROOT=AREA;c.a.fd.ROOT=AREA
SOURCE=RUN/'outputs/matrix_temporal_M1c_native10_secant_20261023.json'
CURVATURE=RUN/'outputs/matrix_temporal_M1c_native10_curvature_20261023.json'
OPTIONS=list(c.VECOPTS);OPTIONS[OPTIONS.index('-pep_tol')+1]='1e-13'
def write(name,data):c.write(name,data)
def snap(path,name):
 raw=path.read_bytes();target=AREA/'inputs'/name;target.parent.mkdir(exist_ok=True);target.write_bytes(raw)
 return dict(source=str(path),snapshot=str(target),sha256=hashlib.sha256(raw).hexdigest())
def read_arrays(path):
 import struct
 data=path.read_bytes();offset=0;arrays=[]
 while offset<len(data):
  cls,n=struct.unpack_from('>ii',data,offset);offset+=8;assert cls==1211214 and 0<n<1000000
  end=offset+16*n;assert end<=len(data);arrays.append(np.frombuffer(data[offset:end],dtype='>c16').astype(np.complex128));offset=end
 return arrays
def valid(r):return c.a.valid(r)
def center(tag,target,origin,reference_q):
 run=c.a.sample(10,origin,target,origin['beta_r']-1e-8,origin['beta_r'],1,'N10_'+tag+'_center',timeout=120,
  extra=dict(options=OPTIONS,seed_shift=1e-6))
 result=dict(tag=tag,target_omega=target,kind='center',run=run,name=run['name'],source_identity_verified=False)
 try:
  info,q=c.vec(run,origin);result['initial_source_vector']=info;result['initial_source_MAC']=c.mac(q,reference_q)
  result['source_identity_verified']=result['initial_source_MAC']>=1-1e-8
  point=c.endpoint(run,target,origin['beta_r']);result['point']=point
  info,q=c.vec(run,point);result['center_vector']=info;result['center_q']=q
 except Exception as e:result['identity_or_endpoint_reason']=repr(e)
 result['eligible_for_arms']=run['complete'] and result['source_identity_verified'] and 'center_q' in result
 return result
def arm(tag,side,origin,reference_q):
 endpoint=-40 if side=='left' else 40
 run=c.a.sample(10,origin,origin['omega_r'],origin['beta_r'],endpoint,1,'N10_'+tag+'_'+side,timeout=120,
  extra=dict(options=OPTIONS,seed_shift=1e-6))
 result=dict(tag=tag,target_omega=origin['omega_r'],kind='arm',side=side,endpoint=endpoint,run=run,name=run['name'],source_identity_verified=False)
 try:
  start=c.endpoint(run,origin['omega_r'],origin['beta_r'])
  result['start_alpha_difference']=abs(c.point_eigen(start)-c.point_eigen(origin))
  info,q=c.vec(run,start);result['start_vector']=info;result['start_MAC']=c.mac(q,reference_q)
  result['source_identity_verified']=result['start_alpha_difference']<1e-7 and result['start_MAC']>=1-1e-8
 except Exception as e:result['identity_or_endpoint_reason']=repr(e)
 return result
def serializable(result):return {k:v for k,v in result.items() if k!='center_q'}
def notice(result):
 run=result['run'];good=[p for p in run['points'] if valid(p)]
 positive=[p for p in good if p['growth']>0]
 endpoint=[p for p in positive if abs(abs(p['beta_r'])-40)<1e-7]
 out=dict(tag=result['tag'],kind=result['kind'],side=result.get('side'),run=result['name'],omega=result['target_omega'],
  complete=run['complete'],source_identity_verified=result['source_identity_verified'],
  valid_samples=len(good),positive_points=positive,positive_endpoints=endpoint,peak=max(good,key=lambda p:p['growth']) if good else None,
  minimum_accepted_MAC=run['metrics']['minimum_accepted_mac'],maximum_accepted_residual=run['metrics']['maximum_accepted_residual'],
  state='positive_exists' if positive else 'sampled_nonpositive' if run['complete'] and result['source_identity_verified'] else 'unknown',
  qualification='Actual positive samples stand independently of a return-to-neutral test. Negative path claims require complete source-consistent tracking. Finite sampled path only.')
 if positive:print('NATIVE10_ACTUAL_POSITIVE '+json.dumps(out),flush=True)
 return out
def main():
 before=c.a.driver.snapshot();assert c.a.driver.sha(c.a.SOLVER)==c.a.SHA
 AREA.mkdir(exist_ok=True);(AREA/'inputs').mkdir(exist_ok=True)
 proof=json.loads(SOURCE.read_text());fresh=proof['native_fresh'];source_dir=pathlib.Path(fresh['config_path']).parent
 snapshots=[snap(SOURCE,SOURCE.name),snap(CURVATURE,CURVATURE.name)]
 for p in [pathlib.Path(fresh['native_vector_file']),pathlib.Path(fresh['config_path']),source_dir/'out/neutral_field.plt',source_dir/'solver.log']:
  if p.exists():snapshots.append(snap(p,'native_source_'+p.name))
 rawfile=pathlib.Path(fresh['native_vector_file']);assert c.sha(rawfile)==fresh['vector_sha256']
 qsource=read_arrays(rawfile)[fresh['raw_index']]
 native_fields=c.a.driver.table(source_dir/'out/neutral_field.plt')
 p=min(native_fields,key=lambda p:abs(p['beta_r']-fresh['fixed_parameters']['beta_r']))
 origin=dict(p,omega_r=fresh['fixed_parameters']['omega_r'],beta_r=fresh['fixed_parameters']['beta_r'],
  alpha_r=fresh['alpha_r'],alpha_i=fresh['alpha_i'],growth=fresh['growth'])
 seedrun,seedinfo,q=c.fresh(10,origin,'N10_independent_actual_seed')
 source_mac=c.mac(q,qsource)
 validation=dict(created=c.a.stamp(),origin=origin,origin_classification='beta_minimum_not_upper_maximum',
  source_matrix_report=str(SOURCE),raw_MAC_to_actual_native_source=source_mac,source_vector_sha256=fresh['vector_sha256'],
  fresh_run=c.brief(seedrun),fresh_vector=seedinfo,passed=source_mac>=1-1e-8,snapshots=snapshots)
 write('source_fresh_join',validation);assert validation['passed'],'Independent actual native seed vector mismatch'
 cuts=[('wplus005',origin['omega_r']+.05),('w40',40.)]
 plan=dict(created=c.a.stamp(),station=10,source_seed=origin,origin_classification='beta_minimum',
  cuts=[dict(tag=t,omega=w,beta_center=origin['beta_r'],beta_domain=[-40,40]) for t,w in cuts],
  beta_db_target=1,continuation_step=.4,minimum_MAC=.95,residual_tol=1e-7,nev=1,retry=8,station_predictor=False,
  center_timeout_seconds=120,arm_timeout_seconds=120,workers=2,solver_sha256=c.a.SHA,
  qualification='Two single-frequency finite beta cuts. Positive endpoint may indicate beta/frequency truncation; never require an interior upper neutral maximum.')
 write('cut_plan',plan);results=[];summaries=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
  pending={pool.submit(center,tag,w,origin,qsource):('center',tag) for tag,w in cuts}
  while pending:
   done,_=concurrent.futures.wait(pending,return_when=concurrent.futures.FIRST_COMPLETED)
   for future in done:
    kind,tag=pending.pop(future)
    try:
     result=future.result();out=serializable(result);write('results/'+result['name'],out);results.append(out)
     summary=notice(result);summaries.append(summary)
     if kind=='center':
      # Center state is saved before submitting either arm.
      write('centers',dict(updated=c.a.stamp(),centers=[r for r in results if r['kind']=='center']))
      if result['eligible_for_arms']:
       for side in ('left','right'):pending[pool.submit(arm,tag,side,result['point'],result['center_q'])]=('arm',tag)
      else:
       for side in ('left','right'):summaries.append(dict(tag=tag,kind='arm',side=side,omega=result['target_omega'],state='unknown',reason='center transport/identity incomplete; arm not attempted'))
    except Exception as e:summaries.append(dict(tag=tag,kind=kind,state='unknown',reason=repr(e)))
    write('cut_progress',dict(created=plan['created'],updated=c.a.stamp(),results=summaries,
      full_results=[dict(file=str(AREA/'results'/(r['name']+'.json')),tag=r['tag'],kind=r['kind'],side=r.get('side')) for r in results],
      pending_tasks=len(pending),finished=not pending,
      qualification=plan['qualification']))
 after=c.a.driver.snapshot();assert before==after;assert c.a.driver.sha(c.a.SOLVER)==c.a.SHA
 write('integrity',dict(before=before,after=after,unchanged=True,solver_sha256=c.a.SHA))
if __name__=='__main__':main()

