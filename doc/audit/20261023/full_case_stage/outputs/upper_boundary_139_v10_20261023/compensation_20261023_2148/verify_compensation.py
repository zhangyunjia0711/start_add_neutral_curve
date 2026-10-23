"""Verify stored returned vectors, then add diagnostic start-point vector joins.
No physical route is re-run. A distinct returned pair never fills the gap.
"""
import collections,csv,hashlib,json,pathlib
import compensate as c
def readraw(stage):return json.loads(pathlib.Path(stage['raw_audit']).read_text())
def csvwrite(name,rows,cols):
 with (c.AREA/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=cols,extrasaction='ignore');w.writeheader();w.writerows(rows)
def main():
 progress=json.loads((c.AREA/'compensation_progress.json').read_text());verified=[];summary=[];cap_positive=[]
 for original in progress['results']:
  out=dict(original);job=out['job'];ref=job['reference'];stages=out['stages'];out['success']=False
  out['verification_note']='Stored actual returned eigenpair vector is selected by its own field lambda, then compared with original reference; no route solver repeated.'
  try:
   assert len(stages)>=2,'No completed alternate route'
   rawref=readraw(stages[0]);refinfo,refq=c.vec(rawref,ref)
   raw=readraw(stages[-1]);assert raw['complete'],'Bounded route timeout/incomplete'
   actual=c.endpoint(raw,ref['omega_r'],ref['beta_r'])
   info,q=c.vec(raw,actual)
   out.update(returned_actual=actual,returned_vector=info,join_alpha=abs(c.point_eigen(actual)-c.point_eigen(ref)),join_mac=c.mac(refq,q))
   matched=out['join_alpha']<1e-8 and out['join_mac']>=1-1e-8
   out['end_pair_matched']=matched;out['start_joins']=[]
   if matched:
    for index,stage in enumerate(stages[1:],1):
     if '_reference_' in stage['name']:continue
     route=readraw(stage);assert route['complete'],'Incomplete component route'
     seed=route['actual_seed'];at=c.endpoint(route,seed['omega_r'],seed['beta_r'])
     delta=abs(c.point_eigen(at)-c.point_eigen(seed));assert delta<1e-8,('Actual start scalar join mismatch',delta)
     start_info,start_q=c.vec(route,at)
     diagnostic,diagnostic_info,dq=c.fresh(job['station'],at,'P139_'+job['name']+f'_stage{index}_start')
     overlap=c.mac(start_q,dq);assert overlap>=1-1e-8,('Actual start raw-vector mismatch',overlap)
     out['start_joins'].append(dict(stage=stage['name'],delta_alpha=delta,mac=overlap,path_start_vector=start_info,
      independent_reference=diagnostic_info,diagnostic=c.brief(diagnostic)))
    assert all(s['complete'] for s in stages),'Incomplete stage'
    out['success']=True;out['state']='vector_verified_compensation'
   else:out['state']='unknown_distinct_returned_eigenpair'
  except Exception as error:
   out['post_verification_reason']=repr(error)
   out['state']='unknown_bounded_timeout' if any(s.get('summary',{}).get('exit')==124 for s in stages) else 'unknown_verification_failed'
  positives=[p for s in stages for p in s.get('positive_points',[]) if abs(p['omega_r']-job['cap'])<1e-8 and p['growth']>0]
  out['cap_positive_points']=positives;cap_positive.extend(dict(station=job['station'],job=job['name'],point=p) for p in positives)
  if positives:out.update(success=False,state='positive_cap_competitor')
  c.write('verified_results/'+job['name'],out);verified.append(out)
  summary.append(dict(station=job['station'],root_id=job['root_id'],side=job['side'],strategy=job['strategy'],
   state=out['state'],success=out['success'],seconds=out['seconds'],join_alpha=out.get('join_alpha'),join_mac=out.get('join_mac'),
   returned_alpha_r=out.get('returned_actual',{}).get('alpha_r'),returned_alpha_i=out.get('returned_actual',{}).get('alpha_i'),
   reference_alpha_r=ref['alpha_r'],reference_alpha_i=ref['alpha_i'],start_join_count=len(out.get('start_joins',[])),
   reason=out.get('post_verification_reason',out.get('reason')),evidence=str(c.AREA/'verified_results'/(job['name']+'.json'))))
 csvwrite('compensation_summary.csv',summary,['station','root_id','side','strategy','state','success','seconds','join_alpha','join_mac',
  'returned_alpha_r','returned_alpha_i','reference_alpha_r','reference_alpha_i','start_join_count','reason','evidence'])
 # Qualified overlay preserves the original frozen 139-station audit.
 original=json.loads((c.PARENT/'boundary_progress.json').read_text());overlay=[];improved=[]
 for station in original['results']:
  n=station['station'];initial=json.loads((c.PARENT/'stations'/f'i{n:03d}.json').read_text())
  failures=initial.get('unknown_paths',[])
  attempts=[v for v in verified if v['job']['station']==n]
  solved={(v['job']['root_id'],v['job']['side']) for v in attempts if v['success']}
  remaining=[p for p in failures if (p['root_id'],p['side']) not in solved]
  qualified=station['state']=='sampled_nonpositive_known_seeds' or bool(failures) and not remaining
  state='sampled_nonpositive_known_seeds' if qualified else 'unknown'
  if any(v['cap_positive_points'] for v in attempts):state='positive_cap_competitor'
  added=state=='sampled_nonpositive_known_seeds' and station['state']=='unknown'
  if added:improved.append(n)
  overlay.append(dict(station=n,x_mm=station['x_mm'],before_state=station['state'],after_state=state,
   qualified_after_vector_compensation=state=='sampled_nonpositive_known_seeds',improved=added,
   remaining_failure_count=len(remaining),remaining_failures=remaining,
   evidence=[str(c.AREA/'verified_results'/(v['job']['name']+'.json')) for v in attempts],
   no_seed_reason=initial.get('reason') if not initial.get('rounds') else None))
 csvwrite('boundary_audit_after_compensation.csv',overlay,['station','x_mm','before_state','after_state','qualified_after_vector_compensation','improved','remaining_failure_count','remaining_failures','no_seed_reason','evidence'])
 report=dict(created=c.a.stamp(),finished=len(verified)==19,before=dict(collections.Counter(s['state'] for s in original['results'])),
  after=dict(collections.Counter(s['after_state'] for s in overlay)),improved_stations=improved,
  remaining_unknown_stations=[r['station'] for r in overlay if r['after_state']=='unknown'],
  unknown_item_counts=dict(collections.Counter(v['state'] for v in verified)),
  vector_compensations=summary,states=overlay,cap_positive_points=cap_positive,
  original_progress_sha256=c.sha(c.PARENT/'boundary_progress.json'),compensation_plan_sha256=c.sha(c.AREA/'compensation_plan.json'),
  solver_sha256=c.a.SHA,integrity=str(c.AREA/'compensation_integrity.json'),
  qualification='Finite B40 sampled paths for explicit known actual roots. Unknown remains unknown; these are not full-spectrum or continuous-beta guarantees. Original audit files remain unchanged.')
 c.write('compensation_report',report)
 c.write('boundary_audit_after_compensation',dict(results=overlay,report=str(c.AREA/'compensation_report.json'),qualification=report['qualification']))
 print(json.dumps({k:report[k] for k in ['before','after','improved_stations','remaining_unknown_stations','unknown_item_counts','cap_positive_points']},ensure_ascii=False))
if __name__=='__main__':main()

