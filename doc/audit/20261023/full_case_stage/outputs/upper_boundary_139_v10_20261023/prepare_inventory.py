"""Snapshot the 139 native-station known-neutral inventory before execution."""
import collections,csv,datetime,hashlib,importlib.util,json,pathlib,sys
AREA=pathlib.Path(__file__).resolve().parent;RUN=AREA.parent.parent
PREVIOUS=RUN/'outputs/upper_boundary_audit_v10_20261023'
sys.path.insert(0,str(RUN/'work'))
from centered_survey import coarse_stations
gridpath=RUN/'outputs/seed_catalog.json'
gridraw=gridpath.read_bytes();catalog=json.loads(gridraw);grid=catalog['grids']['main_841']['x_nondim']
stations=coarse_stations(grid);assert len(stations)==139
AREA.joinpath('inputs').mkdir(parents=True,exist_ok=True)
(AREA/'inputs/seed_catalog.json').write_bytes(gridraw)
module_spec=importlib.util.spec_from_file_location('neutral_collector139',PREVIOUS/'collect_roots.py')
collector=importlib.util.module_from_spec(module_spec);module_spec.loader.exec_module(collector)
collector.AREA=AREA;collector.RUN=RUN;collector.REPRESENTATIVES=stations
inventory=collector.collect()
# Enrichment implementation is retained locally and operates only on hashed
# explicit origin artifacts, matching their actual coordinates before copying.
enrich_source=(PREVIOUS/'enrich_verification.py').read_text()
local=AREA/'enrich_verification.py'
local.write_text(enrich_source)
exec(compile(enrich_source,str(local),'exec'),dict(__file__=str(local),__name__='inventory_enrichment'))
inventory=json.loads((AREA/'root_inventory_enriched.json').read_text())
old=json.loads((PREVIOUS/'B40_progress.json').read_text())
roots=inventory['all_roots'];coverage=[];reuse=[]
for number in stations:
 view=inventory['stations'][str(number)]
 # Select highest verified upper neutral candidate, and separately retain any
 # higher unqualified actual candidate as an explicit limitation.
 positive=[r for r in view['roots'] if r['point']['omega_r']>0]
 eligible=[r for r in positive if r['g_omega_values'] and all(v<0 for v in r['g_omega_values'])]
 highest=max(eligible,key=lambda r:r['point']['omega_r']) if eligible else None
 view.update(highest_positive=highest,upper_derivative_verified=highest is not None,
  target_omega=highest['point']['omega_r']+.05 if highest else None)
 higher=[r['root_id'] for r in positive if highest and r['point']['omega_r']>highest['point']['omega_r']+1e-4]
 view['unqualified_higher_actual_candidates']=higher
 hits=[]
 if highest:
  cap=view['target_omega']
  for root in view['roots']:
   p=root['point']
   for result in old['results']:
    q=result['job']['root']['point']
    same=(result['job']['station']==number and abs(result['job']['target_omega']-cap)<1e-8
     and abs(p['omega_r']-q['omega_r'])<1e-4 and abs(p['beta_r']-q['beta_r'])<1e-3
     and abs(complex(p['alpha_r']-q['alpha_r'],p['alpha_i']-q['alpha_i']))<1e-3)
    if same:
     hits.append(dict(station=number,new_root_id=root['root_id'],old_root_id=result['job']['root']['root_id'],
       side=result['job']['side'],state=result['state'],complete=result['complete'],evidence=str(PREVIOUS/'runs'/result['job']['spec']['name']/'beta_audit.json'),
       evidence_sha256=hashlib.sha256((PREVIOUS/'runs'/result['job']['spec']['name']/'beta_audit.json').read_bytes()).hexdigest()))
 view['reusable_B40_paths']=hits;reuse.extend(hits)
 coverage.append(dict(station=number,x_mm=35*grid[number-1],known_roots=len(view['roots']),positive_roots=len(positive),
  sources=';'.join(sorted(set(p['source'] for r in view['roots'] for p in r['provenance']))),
  verified_upper_roots=len(eligible),highest_omega=highest['point']['omega_r'] if highest else None,
  highest_beta=highest['point']['beta_r'] if highest else None,cap_omega=view['target_omega'],
  state='ready_known_seeds' if highest else 'unknown_no_verified_upper_root',
  unqualified_higher_actual_candidates=';'.join(higher),reusable_complete_paths=sum(h['complete'] for h in hits),
  reusable_positive_paths=sum(h['state']=='positive_competitor_exists' for h in hits)))
inventory['representative_stations']=stations
inventory['grid_snapshot']=dict(path=str(AREA/'inputs/seed_catalog.json'),sha256=hashlib.sha256(gridraw).hexdigest())
inventory['selection_rule']='Highest positive-frequency candidate with recorded verified negative g_omega. Higher unqualified actual candidates remain separately flagged, not suppressed.'
(AREA/'boundary_inventory.json').write_text(json.dumps(inventory,ensure_ascii=False,indent=2)+'\n')
# Greedy farthest physical-distance sampling prioritizes the biggest gaps
# relative to fully reused audits. Preserve end points and then fill gaps.
ready=[r for r in coverage if r['state']=='ready_known_seeds']
covered=[r['x_mm'] for r in ready if r['reusable_complete_paths']>=2*r['known_roots']]
pending=[r for r in ready if r['reusable_complete_paths']<2*r['known_roots']]
priority=[]
while pending:
 item=max(pending,key=lambda r:min((abs(r['x_mm']-x) for x in covered),default=float('inf')))
 item=dict(item,priority=len(priority)+1,distance_to_previous_audit_mm=min((abs(item['x_mm']-x) for x in covered),default=None))
 priority.append(item);covered.append(item['x_mm']);pending=[r for r in pending if r['station']!=item['station']]
for row in coverage:
 row['priority']=next((p['priority'] for p in priority if p['station']==row['station']),None)
with (AREA/'coverage.csv').open('w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=list(coverage[0]));writer.writeheader();writer.writerows(coverage)
with (AREA/'priority.csv').open('w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=list(priority[0]));writer.writeheader();writer.writerows(priority)
# Estimate from completed original B40 real solve time, not a claimed ETA.
durations=[r['summary']['seconds'] for r in old['results'] if r['complete']]
per_path=sum(durations)/len(durations)
unreused=sum(max(0,2*r['known_roots']-r['reusable_complete_paths']) for r in ready)
estimate=dict(observed_complete_paths=len(durations),observed_mean_seconds_per_path=per_path,
 unreused_directional_paths=unreused,mean_based_cpu_worker_seconds=unreused*per_path,
 estimated_minutes_one_worker=unreused*per_path/60,estimated_minutes_two_workers=unreused*per_path/120,
 qualification='Observed prior complete-path time only; new center/chunk setup, FD correction, timeouts and contention can increase elapsed time. Not a guarantee.')
summary=dict(created=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),
 representatives=139,trace_snapshots=len(inventory['inputs']),unique_actual_roots=len(roots),
 ready_stations=len(ready),unknown_stations=[r for r in coverage if r['state']!='ready_known_seeds'],
 source_counts=dict(collections.Counter(p['source'] for r in roots for p in r['provenance'])),
 reused_complete_directional_paths=sum(h['complete'] for h in reuse),reuse=reuse,
 higher_unqualified=[r for r in coverage if r['unqualified_higher_actual_candidates']],estimate=estimate,
 first_20_priority=priority[:20],coverage=coverage,
 qualification='Known actual roots only; finite sampled real-beta audit cannot certify full spectrum. Missing root/verification/path is unknown.')
(AREA/'coverage_plan.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
print('COVERAGE_SUMMARY '+json.dumps({k:v for k,v in summary.items() if k not in ('coverage','reuse','source_counts','unknown_stations')},ensure_ascii=False),flush=True)
print('UNKNOWN_STATIONS '+json.dumps([r['station'] for r in summary['unknown_stations']]),flush=True)

