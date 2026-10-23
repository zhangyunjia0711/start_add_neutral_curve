"""Known-seed, single-frequency B40 competition audit over 139 native stations.
Finite sampled paths only. Missing/failed paths remain unknown. A positive
actual point triggers bounded independent FD and, if verified, a new cap round.
"""
import argparse,concurrent.futures,copy,csv,datetime,hashlib,json,math,os,pathlib,sys,threading,time
AREA=pathlib.Path(__file__).resolve().parent;RUN=AREA.parent.parent
PREVIOUS=RUN/'outputs/upper_boundary_audit_v10_20261023'
sys.path.insert(0,str(PREVIOUS));import audit_followup as a
a.AREA=AREA;a.INVENTORY=AREA/'boundary_inventory.json';a.driver.ROOT=AREA;a.fd.ROOT=AREA
def fd_run(spec):return a.driver.run(dict(spec,timeout=min(spec.get('timeout',60),60)))
a.fd.run=fd_run
LOCK=threading.Lock();EVENTS=[];COMPACT={};INVENTORY_HASH=None
def sha(path):return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()
def write(path,data):
 path=pathlib.Path(path);path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix('.tmp')
 tmp.write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+'\n');tmp.replace(path)
def compact(result):
 return {k:result.get(k) for k in ('station','x_mm','state','initial_cap_omega','final_cap_omega','seed_count','correction_count','seconds','positive_count','unknown_paths','evidence')}
def publish(result):
 number=result['station'];path=AREA/'stations'/f'i{number:03d}.json';write(path,result)
 with LOCK:
  COMPACT[number]=dict(compact(result),evidence=str(path),evidence_sha256=sha(path))
  data=dict(schema_version=1,created=STARTED,updated=a.stamp(),inventory=str(a.INVENTORY),inventory_sha256=INVENTORY_HASH,
   solver_sha256=a.SHA,planned_stations=139,completed_stations=sum(r['state'] not in ('pending','running') for r in COMPACT.values()),
   results=[COMPACT[k] for k in sorted(COMPACT)],positive_events=str(AREA/'positive_events.json'),
   qualification='Finite sampled real-beta paths for explicit known roots. Unknown is not stable; not a complete-spectrum proof.')
  write(AREA/'boundary_progress.json',data)
  with (AREA/'boundary_station_summary.csv').open('w',newline='') as f:
   columns=sorted(set(k for row in data['results'] for k in row));writer=csv.DictWriter(f,fieldnames=columns);writer.writeheader();writer.writerows(data['results'])
def event(data):
 with LOCK:
  data=dict(created=a.stamp(),**data);EVENTS.append(data)
  write(AREA/'positive_events.json',dict(updated=a.stamp(),events=EVENTS))
 print('ACTUAL_POSITIVE_EVENT '+json.dumps(data),flush=True)
def semantic(p,q):return abs(p['omega_r']-q['omega_r'])<1e-4 and abs(p['beta_r']-q['beta_r'])<1e-3 and abs(complex(p['alpha_r']-q['alpha_r'],p['alpha_i']-q['alpha_i']))<1e-3
OLD=json.loads((PREVIOUS/'B40_progress.json').read_text())['results']
def old_paths(number,root,cap):
 return [r for r in OLD if r['job']['station']==number and semantic(root['point'],r['job']['root']['point']) and abs(r['job']['target_omega']-cap)<1e-8]
def root_info(root):
 return dict(root_id=root['root_id'],point=root['point'],sources=sorted(set(p['source'] for p in root['provenance'])),
  inventory_reference=root['root_id'],derived=any(p.get('derived_joint_conjugate') for p in root['provenance']))
def brief(chunk):
 spec=chunk['spec'];d=AREA/'runs'/spec['name'];good=[p for p in chunk['points'] if a.valid(p)]
 return dict(name=spec['name'],state=chunk['state'],complete=chunk['complete'],actual_seed=chunk['actual_seed'],
  samples=len(chunk['points']),valid=len(good),peak=max(good,key=lambda p:p['growth']) if good else None,
  positive_points=chunk['positive_points'],metrics=chunk['metrics'],seconds=chunk['summary']['seconds'],
  exit=chunk['summary']['exit'],audit_file=str(d/'chunk_audit.json'),audit_sha256=sha(d/'chunk_audit.json'))
def positive_refine(number,root,cap,point,name,round_number):
 token=hashlib.sha256(json.dumps(point,sort_keys=True).encode()).hexdigest()[:10]
 ev=dict(station=number,root_id=root['root_id'],sources=root_info(root)['sources'],cap_omega=cap,
  actual_seed=point,run=name,status='actual_positive_pending_FD',round=round_number)
 event(ev)
 result=a.fd.refine(number,point,point['omega_r'],point['beta_r'],a.SOLVER,'wall_local',
  label=f'N139_i{number:03d}_r{round_number}_positiveFD_{token}',h=.02,max_iterations=12,
  omega_bounds=(0,400),beta_bounds=(-40,40),allow_intermediate_curvature=True)
 path=AREA/'refinements'/f'i{number:03d}_round{round_number}_{token}.json'
 output=dict(created=a.stamp(),event=ev,result=result);write(path,output)
 verified=result.get('status')=='candidate_converged' and result['verification']['J'][0][0]<0
 if verified:
  with LOCK:
   export_path=AREA/'boundary_new_roots.json'
   export=json.loads(export_path.read_text()) if export_path.exists() else dict(results=[])
   export['results'].append(dict(job=dict(name=f'boundary139_i{number:03d}_round{round_number}_{token}',station=number,
    source='boundary_competition',origin_sources=root_info(root)['sources'],actual_seed=point,
    solver=str(a.SOLVER),solver_sha256=a.SHA,artifact=str(path)),result=result))
   export.update(updated=a.stamp(),qualification='Verified local upper candidates from actual positive points; separate finite-beta cap audit remains required.')
   write(export_path,export)
 event(dict(station=number,root_id=root['root_id'],cap_omega=cap,status='FD_verified_upper_root' if verified else 'FD_unknown',
  point=result.get('point'),reason=result.get('reason'),artifact=str(path),artifact_sha256=sha(path),round=round_number))
 return result if verified else None,dict(path=str(path),sha256=sha(path),verified=verified)
def run_station(number,inv):
 view=inv['stations'][str(number)];start=time.time();existing=AREA/'stations'/f'i{number:03d}.json'
 if existing.exists():
  prior=json.loads(existing.read_text())
  if prior.get('inventory_sha256')==INVENTORY_HASH and prior['state'] not in ('pending','running'):
   publish(prior);return prior
 x=35*GRID[number-1]
 result=dict(station=number,x_mm=x,state='running',inventory_sha256=INVENTORY_HASH,solver_sha256=a.SHA,
  initial_cap_omega=view['target_omega'],final_cap_omega=view['target_omega'],seed_count=len(view['roots']),
  correction_count=0,positive_count=0,unknown_paths=[],rounds=[],corrections=[],
  unqualified_higher_actual_candidates=view.get('unqualified_higher_actual_candidates',[]),
  qualification='Finite known-source sample paths. MAC within processes; joins require scalar alpha consistency only.')
 if not view['upper_derivative_verified']:
  result.update(state='unknown',reason='No actual positive-frequency root with verified negative g_omega',seconds=time.time()-start)
  publish(result);return result
 roots=sorted(copy.deepcopy(view['roots']),key=lambda r:r['point']['omega_r'],reverse=True)
 cap=view['target_omega'];publish(result)
 for round_number in range(3):
  round_record=dict(round=round_number,cap_omega=cap,roots=[],state='running');result['rounds'].append(round_record)
  change=None;positive_unresolved=False
  for root in roots:
   old=old_paths(number,root,cap);rr=dict(**root_info(root),paths={});round_record['roots'].append(rr)
   centers=[p for r in old for p in r.get('field',[]) if a.valid(p) and abs(p['beta_r']-root['point']['beta_r'])<1e-7]
   if centers:
    center=centers[0];rr['center']=dict(state='reused_actual_B40_center',point=center)
   else:
    p=root['point'];c=a.sample(number,p,cap,p['beta_r'],p['beta_r']+1e-8,1,f'N139_i{number:03d}_r{round_number}_{root["root_id"]}_center')
    rr['center']=brief(c)
    good=[p for p in c['points'] if a.valid(p) and abs(p['beta_r']-root['point']['beta_r'])<1e-10]
    if c['positive_points']:
     result['positive_count']+=len(c['positive_points']);point=max(c['positive_points'],key=lambda p:p['growth'])
     change,proof=positive_refine(number,root,cap,point,c['name'],round_number);result['corrections'].append(proof)
     if change:break
     positive_unresolved=True
    if not c['complete'] or not good:
     rr.update(state='unknown',reason='center transport failed; both beta directions uncomputed')
     result['unknown_paths'].append(dict(round=round_number,root_id=root['root_id'],side='both',reason=rr['reason']))
     publish(result);continue
    center=good[0]
   for side,direction in [('left',-1),('right',1)]:
    prior=next((r for r in old if r['job']['side']==side and r['complete']),None)
    if prior:
     rr['paths'][side]=dict(state=prior['state'],complete=True,kind='reused_complete_B40',
      evidence=str(PREVIOUS/'runs'/prior['job']['spec']['name']/'beta_audit.json'))
     if prior['positive_points']:
      result['positive_count']+=len(prior['positive_points']);point=max(prior['positive_points'],key=lambda p:p['growth'])
      change,proof=positive_refine(number,root,cap,point,prior['job']['spec']['name'],round_number);result['corrections'].append(proof)
      if not change:positive_unresolved=True
     if change:break
     continue
    path=dict(state='running',complete=False,endpoint=direction*40,current_actual=center,chunks=[])
    rr['paths'][side]=path;current=dict(center);chunk_number=0
    while direction*(direction*40-current['beta_r'])>1e-7:
     end=current['beta_r']+direction*min(5,abs(direction*40-current['beta_r']))
     c=a.sample(number,current,cap,current['beta_r'],end,1,f'N139_i{number:03d}_r{round_number}_{root["root_id"]}_{side}_c{chunk_number:02d}')
     small=brief(c);path['chunks'].append(small);chunk_number+=1
     good=[p for p in c['points'] if a.valid(p)]
     same=[p for p in good if abs(p['beta_r']-current['beta_r'])<1e-7]
     delta=abs(complex(same[0]['alpha_r']-current['alpha_r'],same[0]['alpha_i']-current['alpha_i'])) if same else None
     small['join_alpha_difference']=delta
     if c['positive_points']:
      result['positive_count']+=len(c['positive_points']);point=max(c['positive_points'],key=lambda p:p['growth'])
      change,proof=positive_refine(number,root,cap,point,c['name'],round_number);result['corrections'].append(proof)
      if change:break
      positive_unresolved=True
     matched=[p for p in good if abs(p['beta_r']-end)<1e-7]
     if not c['complete'] or not matched or delta is None or delta>=1e-7:
      path.update(state='unknown',reason='incomplete/MAC/runtime/scalar-join failure',last_accepted_beta=current['beta_r'])
      result['unknown_paths'].append(dict(round=round_number,root_id=root['root_id'],side=side,
        from_beta=current['beta_r'],remaining_endpoint=direction*40,reason=path['reason']))
      break
     current=matched[0];path['current_actual']=current
     publish(result)
    if change:break
    if abs(current['beta_r']-direction*40)<1e-7:
     path.update(complete=True,state='positive_competitor_exists' if any(c['positive_points'] for c in path['chunks']) else 'sampled_nonpositive')
   if change:break
   rr['state']='sampled_nonpositive' if len(rr['paths'])==2 and all(p['complete'] and p['state']=='sampled_nonpositive' for p in rr['paths'].values()) else 'unknown'
   publish(result)
  if change:
   newpoint=change['point'];newcap=newpoint['omega_r']+.05
   round_record['state']='cap_corrected_after_actual_positive'
   if newcap<=cap+1e-7 or round_number>=2:
    result.update(state='positive_competitor_unresolved',reason='FD upper root did not raise cap or two correction rounds exhausted');break
   root=dict(root_id=f'fd_i{number}_round{round_number}',station=number,point=newpoint,
    g_omega_values=[change['verification']['J'][0][0]],g_beta_beta_values=[change['verification']['J'][1][1]],
    provenance=[dict(source='actual_positive_FD',derived_joint_conjugate=False,artifact=result['corrections'][-1])])
   roots=[root]+[r for r in roots if not semantic(r['point'],newpoint)]
   cap=newcap;result.update(final_cap_omega=cap,correction_count=round_number+1,seed_count=len(roots))
   publish(result);continue
  complete=all(r.get('state')=='sampled_nonpositive' for r in round_record['roots']) and len(round_record['roots'])==len(roots)
  round_record['state']='sampled_nonpositive_known_seeds' if complete and not positive_unresolved else 'unknown'
  result['state']='positive_competitor_unresolved' if positive_unresolved else round_record['state']
  break
 else:result.update(state='unknown',reason='correction rounds exhausted')
 result['seconds']=round(time.time()-start,3);publish(result)
 print('STATION_END '+json.dumps(compact(result)),flush=True);return result
def monitor(stop,filename):
 samples=[];peaks={}
 while not stop.is_set():
  owned=[]
  for p in pathlib.Path('/proc').iterdir():
   if not p.name.isdigit():continue
   try:
    cwd=(p/'cwd').resolve()
    if AREA/'runs' not in cwd.parents:continue
    cmd=(p/'cmdline').read_bytes()
    if str(a.SOLVER).encode() not in cmd:continue
    rss=next(int(l.split()[1])*1024 for l in (p/'status').read_text().splitlines() if l.startswith('VmRSS:'))
    owned.append(dict(pid=int(p.name),rss_bytes=rss,cwd=str(cwd)));peaks[p.name]=max(peaks.get(p.name,0),rss)
   except (OSError,StopIteration,ValueError):pass
  available=next(int(l.split()[1])*1024 for l in pathlib.Path('/proc/meminfo').read_text().splitlines() if l.startswith('MemAvailable:'))
  samples.append(dict(time=a.stamp(),mem_available_bytes=available,loadavg=os.getloadavg(),owned=owned))
  write(AREA/filename,dict(created=STARTED,updated=a.stamp(),samples=samples,peak_rss_by_pid=peaks,
   max_process_rss_bytes=max(peaks.values(),default=0),minimum_available_bytes=min(s['mem_available_bytes'] for s in samples),
   qualification='Observed every 2 seconds; short-lived process peaks can be missed.'))
  stop.wait(2)
def main():
 global STARTED,INVENTORY_HASH,GRID
 p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=4);p.add_argument('--limit',type=int);p.add_argument('--tag',required=True)
 args=p.parse_args();assert 1<=args.workers<=8
 STARTED=a.stamp();INVENTORY_HASH=sha(a.INVENTORY);inv=json.loads(a.INVENTORY.read_text())
 GRID=json.loads((AREA/'inputs/seed_catalog.json').read_text())['grids']['main_841']['x_nondim']
 before=a.driver.snapshot();assert a.driver.sha(a.SOLVER)==a.SHA
 plan=json.loads((AREA/'coverage_plan.json').read_text())
 order=[r['station'] for r in plan['first_20_priority']]
 # The complete ordering is in priority.csv, not just the preview.
 with (AREA/'priority.csv').open() as f:order=[int(r['station']) for r in csv.DictReader(f)]
 extra=[n for n in inv['representative_stations'] if n not in order]
 selected=order[:args.limit] if args.limit else order+extra
 for n in inv['representative_stations']:
  COMPACT[n]=dict(station=n,x_mm=35*GRID[n-1],state='pending',initial_cap_omega=inv['stations'][str(n)]['target_omega'],
   final_cap_omega=None,seed_count=len(inv['stations'][str(n)]['roots']))
 if (AREA/'positive_events.json').exists():EVENTS.extend(json.loads((AREA/'positive_events.json').read_text())['events'])
 write(AREA/(args.tag+'_execution_plan.json'),dict(created=STARTED,workers=args.workers,selected_stations=selected,
  inventory_sha256=INVENTORY_HASH,solver_sha256=a.SHA,arguments=vars(args)))
 print('RUNNER_STARTED '+json.dumps(dict(tag=args.tag,workers=args.workers,stations=selected)),flush=True)
 sys.stdout=(AREA/(args.tag+'_controller.log')).open('a',buffering=1)
 stop=threading.Event();thread=threading.Thread(target=monitor,args=(stop,args.tag+'_resources.json'),daemon=True);thread.start()
 try:
  with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
   pending={pool.submit(run_station,n,inv):n for n in selected}
   for future in concurrent.futures.as_completed(pending):
    number=pending[future]
    try:future.result()
    except Exception as error:
     result=dict(station=number,x_mm=35*GRID[number-1],state='unknown',reason=repr(error),inventory_sha256=INVENTORY_HASH,
       initial_cap_omega=inv['stations'][str(number)]['target_omega'],final_cap_omega=None,seed_count=len(inv['stations'][str(number)]['roots']))
     publish(result);print('STATION_EXCEPTION '+json.dumps(result),flush=True)
 finally:stop.set();thread.join(timeout=5)
 after=a.driver.snapshot();assert before==after;assert a.driver.sha(a.SOLVER)==a.SHA
 write(AREA/(args.tag+'_integrity.json'),dict(before=before,after=after,unchanged=True,solver_sha256=a.SHA))
if __name__=='__main__':main()

