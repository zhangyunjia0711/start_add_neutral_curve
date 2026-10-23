"""One bounded, vector-checked compensation per original unknown item."""
import argparse,collections,concurrent.futures,csv,datetime,hashlib,json,math,pathlib,re,struct,sys,time
import numpy as np
AREA=pathlib.Path(__file__).resolve().parent;PARENT=AREA.parent;RUN=PARENT.parent.parent
PREVIOUS=RUN/'outputs/upper_boundary_audit_v10_20261023'
sys.path.insert(0,str(PREVIOUS));import audit_followup as a
a.AREA=AREA;a.driver.ROOT=AREA;a.fd.ROOT=AREA
VECOPTS=list(a.BASE['options'])+['-pep_view_vectors','binary:out/eigenvectors.petsc:default:append']
FLOAT=r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?'
PAIR=re.compile(r'^\s*(\d+)\s+\(\s*('+FLOAT+r')\s*,\s*('+FLOAT+r')\s*\)\s+The error is\s+('+FLOAT+r')')
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def write(name,data):
 p=AREA/(name+'.json');p.parent.mkdir(parents=True,exist_ok=True);t=p.with_suffix('.tmp');t.write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+'\n');t.replace(p)
def point_eigen(p):return complex(p['alpha_r'],p['alpha_i'])
def mac(q,p):return float(abs(np.vdot(q,p))**2/(np.vdot(q,q).real*np.vdot(p,p).real))
def points_for(path):
 if path.get('kind')=='reused_complete_B40':
  raw=json.loads(pathlib.Path(path['evidence']).read_text());return raw.get('field',raw.get('points',[]))
 return [p for c in path.get('chunks',[]) for p in json.loads(pathlib.Path(c['audit_file']).read_text()).get('points',[])]
def brief(run):
 return dict(name=run['name'],state=run['state'],complete=run['complete'],summary=run['summary'],metrics=run['metrics'],
  point_count=len(run['points']),actual_seed=run['actual_seed'],positive_points=run['positive_points'],
  raw_audit=str(AREA/'runs'/run['name']/'chunk_audit.json'))
def vec(run,expected):
 assert run['complete'],'incomplete solver path'
 p=AREA/'runs'/run['name'];data=(p/'out/eigenvectors.petsc').read_bytes();arrays=[];offset=0
 while offset<len(data):
  cls,n=struct.unpack_from('>ii',data,offset);offset+=8
  assert cls==1211214 and 0<n<1000000,(cls,n)
  end=offset+16*n;assert end<=len(data);arrays.append(np.frombuffer(data[offset:end],dtype='>c16').astype(np.complex128));offset=end
 pairs=[];groups=[]
 for line in (p/'solver.log').read_text(errors='replace').splitlines():
  if 'ncov=' in line:groups.append([])
  m=PAIR.search(line)
  if m:
   item=dict(index=int(m[1]),alpha_r=float(m[2].replace('D','E')),alpha_i=float(m[3].replace('D','E')),residual=float(m[4].replace('D','E')))
   assert groups,'eigenpair before solver-group marker'
   groups[-1].append(item);pairs.append(item)
 assert len(pairs)==len(arrays)>0,(len(pairs),len(arrays))
 matches=[i for i,pair in enumerate(pairs) if abs(point_eigen(pair)-point_eigen(expected))<1e-8 and pair['residual']<1e-7]
 assert matches,'No exported actual eigenpair matches requested endpoint'
 i=min(matches,key=lambda j:abs(point_eigen(pairs[j])-point_eigen(expected)))
 delta=abs(point_eigen(pairs[i])-point_eigen(expected));assert delta<1e-8 and pairs[i]['residual']<1e-7,(delta,pairs[i])
 q=arrays[i];assert np.isfinite(q).all() and np.vdot(q,q).real>0
 if any(mac(q,arrays[j])<1-1e-8 for j in matches):raise AssertionError('Ambiguous equal-eigenvalue exported vectors')
 return dict(pair=pairs[i],delta_alpha=delta,vector_index=i,vector_count=len(arrays),solve_group_sizes=[len(g) for g in groups],matched_indices=matches,expected_parameters={k:expected[k] for k in ('iloc','omega_r','beta_r')},path=str(p/'out/eigenvectors.petsc'),sha256=sha(p/'out/eigenvectors.petsc')),q
def endpoint(run,omega,beta):
 good=[p for p in run['points'] if a.valid(p) and abs(p['omega_r']-omega)<1e-10 and abs(p['beta_r']-beta)<1e-10]
 assert good,'no actual requested endpoint'
 return good[-1]
def fresh(number,p,label):
 options=[v for v in VECOPTS];options[options.index('-pep_tol')+1]='1e-13';options[options.index('-pep_ncv')+1]='16'
 r=a.sample(number,p,p['omega_r'],p['beta_r']-1e-8,p['beta_r'],1,label,timeout=30,
   extra=dict(options=options,seed_shift=1e-6))
 actual=endpoint(r,p['omega_r'],p['beta_r']);assert abs(point_eigen(actual)-point_eigen(p))<1e-8
 info,q=vec(r,p);return r,info,q
def beta_line(number,start,target,label):
 return a.sample(number,start,start['omega_r'],start['beta_r'],target,1,label,timeout=120,extra=dict(options=VECOPTS))
def omega_reverse(number,start,target,label):
 count=max(2,math.ceil(abs(start['omega_r']-target)/.4)+1)
 return a.sample(number,start,start['omega_r'],start['beta_r']-1e-8,start['beta_r'],1,label,timeout=120,
  extra=dict(options=VECOPTS,outer_parameter='omega',outer_min=min(target,start['omega_r']),
   outer_max=max(target,start['omega_r']),outer_points=count,reverse_outer=target<start['omega_r']))
def build_plan():
 inv=json.loads((PARENT/'boundary_inventory.json').read_text());jobs=[];missing=[]
 for file in sorted((PARENT/'stations').glob('*.json')):
  d=json.loads(file.read_text())
  if d['state']!='unknown':continue
  if not d.get('rounds'):
   missing.append(dict(station=d['station'],reason=d.get('reason'),state='unknown_no_actual_seed'));continue
  round=d['rounds'][-1];roots={r['root_id']:r for r in round['roots']}
  for failure in d['unknown_paths']:
   rid=failure['root_id'];root=roots[rid];side=failure['side'];reference=None
   job=dict(station=d['station'],x_mm=d['x_mm'],root_id=rid,side=side,cap=round['cap_omega'],failure=failure,
    source_station_file=str(file),source_station_sha256=sha(file),root=root,
    qualification='One bounded trial. Reuse is allowed only after actual-vector MAC and residual/scalar/parameter joins.')
   if side=='both':
    origin=root['point'];candidates=[]
    for other in roots.values():
     if other['root_id']==rid:continue
     for ps in other['paths'].values():
      for p in points_for(ps):
       if a.valid(p):candidates.append((abs(p['beta_r']-origin['beta_r']),other['root_id'],p))
    if not candidates:job.update(strategy='no_alternate_actual_cap_point',executable=False)
    else:
     _,oid,p=min(candidates,key=lambda x:x[0]);job.update(strategy='center_reverse_omega',executable=True,
       alternate_root_id=oid,start=p,reference=origin,original_beta=origin['beta_r'])
   else:
    path=root['paths'][side];reference=path['current_actual'];target=path['endpoint'];candidates=[]
    highest=inv['stations'][str(d['station'])]['highest_positive']['root_id']
    for other in roots.values():
     if other['root_id']==rid:continue
     ps=other['paths'].get(side,{})
     if not ps.get('complete'):continue
     for p in points_for(ps):
      if a.valid(p) and abs(p['beta_r']-target)<1e-7:
       candidates.append((other['root_id']!=highest,-p['growth'],other['root_id'],p))
    if candidates:
     _,_,oid,p=min(candidates,key=lambda x:(x[0],x[1]));job.update(strategy='reverse_endpoint',executable=True,
       alternate_root_id=oid,start=p,reference=reference,target=target)
    else:
     other_side='left' if side=='right' else 'right';own=root['paths'].get(other_side,{})
     start=own.get('current_actual')
     if start is None:
      vals=[]
      for ps in root['paths'].values():vals.extend(points_for(ps))
      vals=[p for p in vals if a.valid(p)]
      start=max(vals,key=lambda p:abs(p['beta_r']-reference['beta_r'])) if vals else reference
     job.update(strategy='own_actual_out_and_back',executable=True,start=start,reference=start,target=target,
       original_failure_join=reference)
   job['name']=f'i{job["station"]:03d}_{rid}_{side}';jobs.append(job)
 plan=dict(created=a.stamp(),inventory_sha256=sha(PARENT/'boundary_inventory.json'),
  original_progress_sha256=sha(PARENT/'boundary_progress.json'),jobs=jobs,missing=missing,
  counts=dict(collections.Counter(j['strategy'] for j in jobs)),workers_max=4,
  timeout_bounds=dict(fresh_seconds=30,beta_line_seconds=120,omega_reverse_seconds=120),
  estimate_minutes=dict(expected=(3,8),conservative_upper=15),
  qualification='No rerun of 122 passed stations; no attempt to fabricate seeds for i1/i2. Parent core and immutable v10 remain unchanged.')
 write('compensation_plan',plan);return plan
def execute(job):
 out=dict(job=job,started=a.stamp(),state='unknown',success=False,stages=[]);t=time.time()
 if not job['executable']:out['reason']=job['strategy'];write('results/'+job['name'],out);return out
 n=job['station'];label='C139_'+job['name'];reference=job['reference']
 try:
  ref,refinfo,refq=fresh(n,reference,label+'_reference');out['stages'].append(brief(ref));out['reference_vector']=refinfo
  if job['strategy']=='reverse_endpoint':
   line=beta_line(n,job['start'],reference['beta_r'],label+'_reverse');out['stages'].append(brief(line))
   actual=endpoint(line,reference['omega_r'],reference['beta_r']);info,q=vec(line,actual)
   out.update(join_alpha=abs(point_eigen(actual)-point_eigen(reference)),join_mac=mac(refq,q),returned_vector=info)
   out['success']=line['complete'] and out['join_alpha']<1e-8 and out['join_mac']>=1-1e-8
  elif job['strategy']=='own_actual_out_and_back':
   outward=beta_line(n,job['start'],job['target'],label+'_outward');out['stages'].append(brief(outward))
   far=endpoint(outward,reference['omega_r'],job['target'])
   returned=beta_line(n,far,reference['beta_r'],label+'_return');out['stages'].append(brief(returned))
   actual=endpoint(returned,reference['omega_r'],reference['beta_r']);info,q=vec(returned,actual)
   out.update(join_alpha=abs(point_eigen(actual)-point_eigen(reference)),join_mac=mac(refq,q),returned_vector=info)
   out['success']=outward['complete'] and returned['complete'] and out['join_alpha']<1e-8 and out['join_mac']>=1-1e-8
  elif job['strategy']=='center_reverse_omega':
   capline=beta_line(n,job['start'],reference['beta_r'],label+'_cap_beta_bridge');out['stages'].append(brief(capline))
   atcap=endpoint(capline,job['cap'],reference['beta_r']);capinfo,capq=vec(capline,atcap)
   capfresh,caprefinfo,caprefq=fresh(n,atcap,label+'_cap_join_reference');out['stages'].append(brief(capfresh))
   out['cap_beta_join_mac']=mac(capq,caprefq);assert out['cap_beta_join_mac']>=1-1e-8
   reverse=omega_reverse(n,atcap,reference['omega_r'],label+'_reverse_omega');out['stages'].append(brief(reverse))
   actual=endpoint(reverse,reference['omega_r'],reference['beta_r']);info,q=vec(reverse,actual)
   out.update(join_alpha=abs(point_eigen(actual)-point_eigen(reference)),join_mac=mac(refq,q),returned_vector=info,
    cap_actual=atcap,cap_vector=capinfo)
   out['success']=reverse['complete'] and out['join_alpha']<1e-8 and out['join_mac']>=1-1e-8
  out['state']='vector_verified_compensation' if out['success'] else 'unknown_pair_mismatch'
 except Exception as error:out['reason']=repr(error)
 positive=[p for stage in out['stages'] for p in stage.get('positive_points',[]) if p['growth']>0]
 # Neutral-reference roundoff is separated from robust cap-frequency competition.
 cap_positive=[p for p in positive if abs(p['omega_r']-job['cap'])<1e-8]
 out['cap_positive_points']=cap_positive;out['seconds']=round(time.time()-t,3)
 write('results/'+job['name'],out)
 print('COMPENSATION_RESULT '+json.dumps(dict(name=job['name'],state=out['state'],success=out['success'],seconds=out['seconds'],
  join_mac=out.get('join_mac'),join_alpha=out.get('join_alpha'),cap_positive=cap_positive,reason=out.get('reason'))),flush=True)
 return out
def main():
 p=argparse.ArgumentParser();p.add_argument('--execute',action='store_true');p.add_argument('--workers',type=int,default=4)
 args=p.parse_args();assert 1<=args.workers<=4;AREA.mkdir(exist_ok=True)
 existing=AREA/'compensation_plan.json'
 plan=json.loads(existing.read_text()) if existing.exists() else build_plan()
 if not args.execute:
  print(json.dumps(dict(jobs=len(plan['jobs']),missing=plan['missing'],counts=plan['counts'],estimate=plan['estimate_minutes'],
   choices=[dict(name=j['name'],strategy=j['strategy'],alternate_root=j.get('alternate_root_id'),start_beta=j.get('start',{}).get('beta_r'),
    reference_beta=j.get('reference',{}).get('beta_r'),cap=j['cap']) for j in plan['jobs']])));return
 before=a.driver.snapshot();assert a.driver.sha(a.SOLVER)==a.SHA;results=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
  futures={pool.submit(execute,j):j for j in plan['jobs']}
  for future in concurrent.futures.as_completed(futures):
   j=futures[future]
   try:r=future.result()
   except Exception as e:r=dict(job=j,state='unknown',success=False,reason=repr(e))
   results.append(r);write('compensation_progress',dict(updated=a.stamp(),results=results,completed=len(results),planned=len(plan['jobs'])))
 after=a.driver.snapshot();assert before==after;assert a.driver.sha(a.SOLVER)==a.SHA
 write('compensation_integrity',dict(before=before,after=after,unchanged=True,solver_sha256=a.SHA))
if __name__=='__main__':main()

