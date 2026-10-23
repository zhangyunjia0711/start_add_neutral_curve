"""Reduce real boundary-audit fields without treating missing paths as stable."""
import collections,csv,datetime,hashlib,json,pathlib,re,sys
AREA=pathlib.Path(__file__).resolve().parent;RUN=AREA.parent.parent
sys.path.insert(0,str(RUN/'work'));import full_driver
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def dump(name,data): (AREA/name).write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def csvout(name,rows,columns):
 with (AREA/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=columns,extrasaction='ignore');w.writeheader();w.writerows(rows)
def fields(path):
 d=json.loads(pathlib.Path(path).read_text());return d.get('field',d.get('points',[]))
def main():
 inv=json.loads((AREA/'boundary_inventory.json').read_text());progress=json.loads((AREA/'boundary_progress.json').read_text())
 current={r['station']:r for r in progress['results']};stations=[];samples=[];path_rows=[];joins=[];full=[]
 for number in inv['representative_stations']:
  summary=current[number];file=AREA/'stations'/f'i{number:03d}.json'
  if not file.exists():
   stations.append(dict(station=number,x_mm=summary['x_mm'],state='pending',unknown=True));continue
  data=json.loads(file.read_text());full.append(data);final_round=len(data.get('rounds',[]))-1;final_points=[]
  for round_rec in data.get('rounds',[]):
   for root in round_rec['roots']:
    sources=';'.join(root['sources'])
    for side,path in root['paths'].items():
     sets=[]
     if path.get('kind')=='reused_complete_B40':sets=[dict(name='reused',evidence=path['evidence'])]
     else:sets=[dict(name=c['name'],evidence=c['audit_file']) for c in path.get('chunks',[])]
     for item in sets:
      values=fields(item['evidence'])
      for point in values:
       samples.append(dict(station=number,x_mm=data['x_mm'],round=round_rec['round'],is_final_round=round_rec['round']==final_round,
        root_id=root['root_id'],sources=sources,side=side,run=item['name'],evidence=item['evidence'],**point))
      if round_rec['round']==final_round:final_points.extend(p for p in values if p['status']==0)
     for c in path.get('chunks',[]):
      if c.get('join_alpha_difference') is not None:
       joins.append(dict(station=number,run=c['name'],delta_alpha=c['join_alpha_difference'],passed=c['join_alpha_difference']<1e-7))
     path_rows.append(dict(station=number,x_mm=data['x_mm'],round=round_rec['round'],cap_omega=round_rec['cap_omega'],
      root_id=root['root_id'],sources=sources,side=side,state=path['state'],complete=path['complete'],
      kind=path.get('kind','chunked'),current_beta=path.get('current_actual',{}).get('beta_r'),
      target_endpoint=path.get('endpoint'),evidence=str(file)))
    if not root['paths']:
     path_rows.append(dict(station=number,x_mm=data['x_mm'],round=round_rec['round'],cap_omega=round_rec['cap_omega'],
      root_id=root['root_id'],sources=sources,side='both',state='unknown',complete=False,kind='center_transport_failed',evidence=str(file)))
  peak=max(final_points,key=lambda p:p['growth']) if final_points else {}
  stations.append(dict(station=number,x_mm=data['x_mm'],state=data['state'],
   highest_known_neutral_omega=(data.get('final_cap_omega')-.05) if data.get('final_cap_omega') is not None else None,
   initial_cap_omega=data.get('initial_cap_omega'),final_cap_omega=data.get('final_cap_omega'),seed_count=data.get('seed_count'),
   correction_count=data.get('correction_count'),positive_count=data.get('positive_count'),seconds=data.get('seconds'),
   final_cap_maximum_sampled_growth=peak.get('growth'),final_cap_peak_beta=peak.get('beta_r'),
   final_round_qualified=data['state']=='sampled_nonpositive_known_seeds',unknown=data['state'] in ('unknown','positive_competitor_unresolved'),
   evidence=str(file),evidence_sha256=sha(file)))
 csvout('boundary_audit_stations.csv',stations,['station','x_mm','state','highest_known_neutral_omega','initial_cap_omega','final_cap_omega',
  'seed_count','correction_count','positive_count','seconds','final_cap_maximum_sampled_growth','final_cap_peak_beta','final_round_qualified','unknown','evidence','evidence_sha256'])
 csvout('boundary_audit_paths.csv',path_rows,['station','x_mm','round','cap_omega','root_id','sources','side','state','complete','kind','current_beta','target_endpoint','evidence'])
 csvout('boundary_audit_samples.csv',samples,['station','x_mm','round','is_final_round','root_id','sources','side','run','evidence']+full_driver.COLUMNS+['row'])
 resources=[]
 for name in ('pilot4a_resources.json','full8_resources.json'):
  p=AREA/name
  if not p.exists():continue
  q=json.loads(p.read_text());resources.append(dict(path=str(p),sha256=sha(p),samples=len(q['samples']),
   max_process_rss_MiB=q['max_process_rss_bytes']/2**20,min_available_GiB=q['minimum_available_bytes']/2**30,
   maximum_observed_load1=max(s['loadavg'][0] for s in q['samples'])))
 events=json.loads((AREA/'positive_events.json').read_text()) if (AREA/'positive_events.json').exists() else dict(events=[])
 newroots=json.loads((AREA/'boundary_new_roots.json').read_text())
 done=all(r['state'] not in ('pending','running') for r in stations)
 output=dict(schema_version=1,created=datetime.datetime.now().astimezone().isoformat(),finished=done,stations=stations,
  states=dict(collections.Counter(r['state'] for r in stations)),input_inventory_sha256=sha(AREA/'boundary_inventory.json'),
  solver_sha256='634b45c5d30e076f5b86ccb136e0734c1905beedde108477dceb8a75adc3853c',
  exported_samples=len(samples),unknown_paths=[r for r in path_rows if r['state']=='unknown'],
  positive_events=events,new_roots=newroots,scalar_joins=joins,all_recorded_scalar_joins_passed=all(j['passed'] for j in joins),
  resources=resources,integrity=str(AREA/'full8_integrity.json') if (AREA/'full8_integrity.json').exists() else None,
  qualification='All claims apply only to these finite sampled real-beta paths and known actual seeds. Failures/missing roots remain unknown; continuous beta and full spectrum are not certified.')
 dump('boundary_audit_report.json',output)
 files=[]
 for p in sorted(AREA.iterdir()):
  if p.is_file() and p.name!='boundary_audit_manifest.json':files.append(dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size))
 dump('boundary_audit_manifest.json',dict(created=output['created'],finished=done,files=files,
  snapshots=[dict(path=str(p),sha256=sha(p)) for p in sorted((AREA/'inputs').rglob('*')) if p.is_file()],
  raw_evidence='stations/ and runs/ retain all numerical and failed attempts; controller_error_before_solver/ is explicitly non-numerical controller failure.'))
 print(json.dumps(dict(finished=done,states=output['states'],samples=len(samples),unknown_paths=len(output['unknown_paths']),new_roots=len(newroots['results']))),flush=True)
if __name__=='__main__':main()

