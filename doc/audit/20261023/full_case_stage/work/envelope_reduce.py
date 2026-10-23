#!/usr/bin/env python3
"""Read-only six-start survey reduction; NEVER runs a solver.

CLI: python3 work/envelope_reduce.py --phase pilot --coordinates global
     python3 work/envelope_reduce.py --phase vectorpilot --coordinates wall_local
Repeat --phase to include several phases without merging their rows. Optional
--root and --output-prefix select the input task folder and output file prefix.
JSON plus eight CSV tables use preliminary schema envelope_reduce/1.0.
known_G is an actual-sample lower bound, never a certified continuous-beta or
full-spectrum maximum. Omega segments describe sampled nodes, not filled regions.
Root-refinement inputs: neutral_brackets[*].left/right.seed and referenced
beta_peak_candidates contain actual omega_r/beta_r/alpha_r/alpha_i and provenance.
"""
from __future__ import annotations
import argparse, configparser, csv, datetime as dt, hashlib, json, math, re
from pathlib import Path
from collections import defaultdict
ROOT=Path(__file__).resolve().parent.parent
SCHEMA='envelope_reduce/1.1'
SOURCES=tuple(f'{family}{suffix}' for family in ('M1','M2','M3') for suffix in ('p','c'))
COLUMNS=('outer','inner','iloc','x','Re','alpha_r','alpha_i','beta_r','beta_i','omega_r','omega_i','growth','branch','status')
PARAM_TOL=1e-8
SIGN_TOL=1e-8
PEAK_TOL=1e-12

def stamp():return dt.datetime.now(dt.timezone(dt.timedelta(hours=8))).isoformat()
def finite(v):return isinstance(v,(int,float)) and math.isfinite(v)
def key(v):return format(float(v),'.12g')
def unique(vs):return sorted({key(v):float(v) for v in vs if finite(v)}.values())
def clean(v):
    if isinstance(v,float) and not math.isfinite(v):return None
    if isinstance(v,dict):return {k:clean(x) for k,x in v.items() if not k.startswith('_')}
    if isinstance(v,(list,tuple)):return [clean(x) for x in v]
    return v

def read_json(path,issues):
    try:return json.loads(path.read_text(encoding='utf-8-sig'))
    except (ValueError,OSError) as e:
        issues.append(dict(path=str(path),issue=repr(e)));return None

def digest(path):
    if not path.is_file():return None
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()

def match(v,grid):
    if not grid:return float(v)
    nearest=min(grid,key=lambda w:abs(w-v))
    return nearest if abs(nearest-v)<=PARAM_TOL*max(1.,abs(v),abs(nearest)) else float(v)

def axis(spec,parameter):
    if parameter=='omega':
        if spec.get('outer_parameter','omega')!='omega':return []
        lo,hi,n=(spec.get('outer_min',spec.get('omega_min')),spec.get('outer_max',spec.get('omega_max')),spec.get('outer_points',spec.get('omega_points')))
    else:
        if spec.get('inner_parameter','beta')!='beta':return []
        lo,hi,n=(spec.get('inner_min',spec.get('beta_min')),spec.get('inner_max',spec.get('beta_max')),spec.get('inner_points',spec.get('beta_points')))
    if not all(finite(v) for v in (lo,hi,n)) or int(n)<1:return []
    n=int(n);return unique([float(lo)+j*(float(hi)-float(lo))/max(1,n-1) for j in range(n)])

def source_label(spec):
    family=spec.get('family')
    if family in ('M1','M2','M3','M3near','M3_local_a','M3_local_b'):
        return family+('c' if spec.get('conjugate',False) else 'p')
    m=re.search(r'_(M[123])([pc])_i\d+',spec.get('name',''))
    return m[1]+m[2] if m else None

def table(path):
    result,bad,zone=[],[],0
    if not path.exists():return result,bad
    for line_number,line in enumerate(path.read_text(errors='replace').splitlines(),1):
        m=re.search(r'ZONE\s+T="row_(\d+)"',line,re.I)
        if m:zone=int(m[1]);continue
        text=line.strip()
        if not text or text.upper().startswith(('VARIABLES','ZONE','TITLE','#')):continue
        try:
            values=[float(x.replace('D','E').replace('d','e')) for x in text.split()]
            if len(values)!=14:raise ValueError('expected 14 columns, got '+str(len(values)))
            result.append(dict(zip(COLUMNS,values),field_row=zone,field_line=line_number))
        except ValueError as e:bad.append(dict(line=line_number,reason=str(e),raw=line))
    return result,bad

def load_run(root,phase,plan,progress,coordinates,issues):
    name=plan['name'];directory=root/'runs'/name
    summary_path=directory/'summary.json';summary=read_json(summary_path,issues)
    summary=summary if isinstance(summary,dict) else None
    spec=dict(plan)
    if summary:spec.update(summary.get('spec',{}))
    config=configparser.ConfigParser()
    try:config.read(directory/'config.ini')
    except configparser.Error as e:issues.append(dict(path=str(directory/'config.ini'),issue=repr(e)))
    evidence={}
    if plan.get('coordinates'):evidence['plan']=plan['coordinates']
    if summary and summary.get('spec',{}).get('coordinates'):evidence['summary_spec']=summary['spec']['coordinates']
    if config.has_option('lst','profile_coordinates'):evidence['actual_config']=config.get('lst','profile_coordinates').strip()
    if (directory/'wall_local_profile_metadata.csv').exists():evidence['wall_local_metadata']='wall_local'
    models=set(evidence.values())
    actual_model=next(iter(models)) if len(models)==1 else 'conflict' if models else 'global'
    if not evidence:evidence['legacy_default']='global'
    selected=actual_model==coordinates;quality=[]
    if summary is None:quality.append('summary_missing_or_unreadable')
    else:
        if summary.get('spec',{}).get('name',name)!=name:quality.append('summary_name_mismatch')
        if summary.get('petsc_errors',0):quality.append('PETSc_or_SLEPc_error')
        if summary.get('unused_options',0):quality.append('unused_solver_options')
        if summary.get('coordinate_metadata_missing',False):quality.append('required_coordinate_metadata_missing')
        if summary.get('exit')==124:quality.append('run_timeout')
    if not selected:quality.append('coordinate_model_filtered_or_conflicting')
    logfile=directory/'solver.log';log=logfile.read_text(errors='replace') if logfile.exists() else ''
    if re.search(r'SIGSEGV|segmentation fault|Fortran runtime error|floating point exception',log,re.I):quality.append('fatal_runtime_error_in_run')
    if re.search(r'PETSC ERROR|SLEPC ERROR',log,re.I):quality.append('PETSc_or_SLEPc_error')
    if re.search(r'Option left:|options you set that were not used',log,re.I):quality.append('unused_solver_options')
    tracking={k:config.get('neutral',k,fallback=None) for k in ('adaptive_steps','vector_tracking','minimum_overlap','minimum_step','continuation_step','eigen_error_tol','track_initial_target')}
    path=directory/'out/neutral_field.plt';field,bad=table(path) if selected else ([],[])
    if bad:quality.append('malformed_field_lines')
    if summary and selected and summary.get('field_records') is not None and summary['field_records']!=len(field):quality.append('summary_field_record_count_mismatch')
    source=source_label(spec);station=int(spec['station']);run_id=phase+':'+name;samples=[]
    for values in field:
        reasons=list(dict.fromkeys(quality));raw=values['status']
        status=int(round(raw)) if finite(raw) else None
        if status!=0:reasons.append('solver_field_status_nonzero_or_nonfinite')
        if not all(finite(values[k]) for k in ('iloc','omega_r','beta_r','alpha_r','alpha_i','growth')):reasons.append('nonfinite_required_field')
        if finite(values['iloc']) and abs(values['iloc']-station)>1e-7:reasons.append('station_mismatch')
        samples.append(dict(preliminary=True,sample_id=f"{run_id}:L{values['field_line']}",phase=phase,
            coordinates=actual_model,station=station,source=source,run_id=run_id,run_name=name,field_path=str(path),
            **values,status_code=status,accepted_for_known_growth=not reasons,exclusion_reasons=reasons))
    return dict(preliminary=True,run_id=run_id,name=name,phase=phase,station=station,source=source,
        coordinates=actual_model,coordinate_evidence=evidence,selected_coordinate_model=selected,
        summary_path=str(summary_path),field_path=str(path),summary_sha256=digest(summary_path),field_sha256=digest(path),
        solver_sha256=summary.get('solver_sha256') if summary else spec.get('solver_sha256'),
        tracking_configuration=tracking,seed_provenance=spec.get('seed_provenance'),
        exit=summary.get('exit') if summary else None,terminal=summary.get('terminal') if summary else 'missing summary',
        partial_or_failed=not summary or summary.get('terminal') not in ('converged','no bracketed root','sampled'),
        quality_issues=list(dict.fromkeys(quality)),malformed_field_lines=bad,actual_samples=len(samples),
        eligible_samples=sum(v['accepted_for_known_growth'] for v in samples),progress_status=progress.get('status'),
        progress_exception=progress.get('exception'),continuity_flags=progress.get('continuity_flags',[]),
        _spec=spec,_omegas=axis(spec,'omega') or unique([v.get('omega') for v in progress.get('rows',[])]),_betas=axis(spec,'beta'),_samples=samples)

def load_evidence(root,phases,coordinates):
    issues,manifests,runs=[],[],[]
    for phase in phases:
        pp=root/'outputs'/f'{phase}_progress.json';pl=root/'outputs'/f'{phase}_plan.json'
        progress=read_json(pp,issues) or {};plan=read_json(pl,issues) or []
        results={v['name']:v for v in progress.get('results',[]) if 'name' in v}
        # Centered surveys split each source into frequency/beta arms. Preserve
        # every arm and use the explicitly declared family set for this phase.
        plan_items=plan.get('jobs',[]) if isinstance(plan,dict) else plan
        families=plan.get('families') if isinstance(plan,dict) else None
        expected_sources=tuple(f+s for f in families for s in ('p','c')) if families else SOURCES
        plans={v['name']:v for v in plan_items if isinstance(v,dict) and 'name' in v}
        for name,result in results.items():plans.setdefault(name,{k:result[k] for k in ('name','family','conjugate','station') if k in result})
        stations=unique(progress.get('stations',[])+[v.get('station') for v in plans.values()])
        for station in stations:
            for source in expected_sources:
                if not any(v.get('station')==station and source_label(v)==source for v in plans.values()):
                    name=f'{phase}_{source}_i{int(station):03d}'
                    plans[name]=dict(name=name,station=int(station),family=source[:-1],conjugate=source.endswith('c'))
        manifests.append(dict(phase=phase,progress_path=str(pp),plan_path=str(pl),progress_sha256=digest(pp),plan_sha256=digest(pl),
            progress_completed=progress.get('completed'),progress_planned=progress.get('planned'),run_slots=len(plans),expected_sources=list(expected_sources)))
        for name,spec in sorted(plans.items()):
            if 'station' not in spec or source_label(spec) is None:issues.append(dict(run=name,issue='unidentified station/source'));continue
            record=load_run(root,phase,spec,results.get(name,{}),coordinates,issues)
            record['expected_sources']=list(expected_sources)
            runs.append(record)
    return runs,manifests,issues

def best(samples):
    eligible=[v for v in samples if v['accepted_for_known_growth']]
    return max(eligible,key=lambda v:v['growth']) if eligible else None

def seed(sample):return {k:sample[k] for k in ('omega_r','beta_r','alpha_r','alpha_i')} if sample else None

def state(g,complete):
    if g is None:return 'unknown'
    if g>SIGN_TOL:return 'positive_exists'
    if not complete:return 'unknown_nonpositive_observed'
    return 'sampled_zero' if abs(g)<=SIGN_TOL else 'sampled_negative'

def point(sample):
    if sample is None:return None
    return dict(sample_id=sample['sample_id'],run_id=sample['run_id'],source=sample['source'],growth=sample['growth'],
        beta=sample['beta_r'],seed=seed(sample),field_path=sample['field_path'],field_line=sample['field_line'])

def local_peaks(row,cells,by_id):
    candidates=[];expected=[j for j,c in enumerate(cells) if c['expected']]
    edges={expected[0],expected[-1]} if expected else set()
    for j,cell in enumerate(cells):
        if not cell['best_sample_id']:continue
        sample=by_id[cell['best_sample_id']]
        neighbors=[cells[i] for i in (j-1,j+1) if 0<=i<len(cells)]
        known=[by_id[c['best_sample_id']] for c in neighbors if c['best_sample_id']]
        missing=any(not c['best_sample_id'] for c in neighbors)
        local=all(sample['growth']>=v['growth']-PEAK_TOL for v in known)
        if not local and j not in edges:continue
        kind='domain_endpoint' if j in edges else 'unknown_neighbor_peak_candidate' if missing else 'interior_local_maximum'
        candidates.append(dict(preliminary=True,peak_id=f"{row['row_id']}:P{j}",source_row_id=row['row_id'],
            phase=row['phase'],coordinates=row['coordinates'],station=row['station'],omega=row['omega'],source=row['source'],
            beta=cell['beta'],growth=sample['growth'],alpha_r=sample['alpha_r'],alpha_i=sample['alpha_i'],kind=kind,
            is_local_max_on_available_neighbors=local,is_domain_endpoint=j in edges,unknown_neighbor=missing,
            all_requested_cells_accepted=row['all_requested_cells_accepted'],sample_id=sample['sample_id'],seed=seed(sample),
            field_path=sample['field_path'],field_line=sample['field_line'],continuous_beta_max_certified=False))
    return candidates

def make_rows(runs,coordinates):
    selected=[r for r in runs if r['selected_coordinate_model']]
    samples=[v for r in selected for v in r['_samples']];by_id={v['sample_id']:v for v in samples}
    station_groups=defaultdict(list)
    for run in selected:station_groups[run['phase'],run['station']].append(run)
    rows,peaks,unknown=[],[],[]
    for (phase,station),station_runs in sorted(station_groups.items()):
        expected_sources=station_runs[0].get('expected_sources',list(SOURCES))
        omega_grid=unique([w for r in station_runs for w in r['_omegas']])
        omega_grid=unique(omega_grid+[match(v['omega_r'],omega_grid) for r in station_runs for v in r['_samples'] if finite(v['omega_r'])])
        station_betas=unique([b for r in station_runs for b in r['_betas']])
        for omega in omega_grid:
            for source in expected_sources:
                source_runs=[r for r in station_runs if r['source']==source]
                planned=[r for r in source_runs if any(key(w)==key(omega) for w in r['_omegas'])]
                beta_grid=unique([b for r in planned for b in r['_betas']]);grid_known=bool(beta_grid)
                if not beta_grid:beta_grid=unique([b for r in source_runs for b in r['_betas']]) or station_betas
                values=[v for r in source_runs for v in r['_samples'] if finite(v['omega_r']) and key(match(v['omega_r'],omega_grid))==key(omega)]
                buckets=defaultdict(list);unlocated=[]
                for value in values:
                    if finite(value['beta_r']):buckets[key(match(value['beta_r'],beta_grid))].append(value)
                    else:unlocated.append(value['sample_id'])
                beta_all=unique(beta_grid+[float(b) for b in buckets]);cells=[]
                for beta in beta_all:
                    actual=buckets.get(key(beta),[]);winner=best(actual)
                    all_good=bool(actual) and all(v['accepted_for_known_growth'] for v in actual)
                    expected=any(key(b)==key(beta) for b in beta_grid)
                    status='accepted' if all_good else 'accepted_with_failed_or_unverified_repeat' if winner else 'failed_or_unverified' if actual else 'missing'
                    cells.append(dict(beta=beta,expected=expected,status=status,sample_ids=[v['sample_id'] for v in actual],
                        best_sample_id=winner['sample_id'] if winner else None,known_growth=winner['growth'] if winner else None,all_evaluations_accepted=all_good))
                complete=bool(grid_known and cells and all(c['all_evaluations_accepted'] for c in cells) and not unlocated)
                winner=best(values);growth=winner['growth'] if winner else None
                rid=f'{phase}:{coordinates}:i{station:03d}:{source}:w{key(omega)}'
                row=dict(preliminary=True,row_id=rid,phase=phase,coordinates=coordinates,station=station,omega=omega,source=source,
                    run_ids=[r['run_id'] for r in source_runs],actual_sample_count=len(values),eligible_sample_count=sum(v['accepted_for_known_growth'] for v in values),
                    requested_frequency=bool(planned),requested_beta_grid_known=grid_known,expected_beta_count=sum(c['expected'] for c in cells),
                    unknown_beta_count=sum(not c['all_evaluations_accepted'] for c in cells)+(not cells),unknown_expected_beta_count=sum(c['expected'] and not c['all_evaluations_accepted'] for c in cells),all_requested_cells_accepted=complete,
                    sampled_grid_max_is_capped=complete,continuous_beta_max_certified=False,spectral_complete=False,known_G=growth,
                    beta_star=winner['beta_r'] if winner else None,winner=point(winner),evidence_state=state(growth,complete),
                    raw_status_counts={str(code):sum(v['status_code']==code for v in values) for code in sorted({v['status_code'] for v in values},key=str)},
                    cells=cells,unlocated_sample_ids=unlocated,expected_sources=expected_sources)
                for cell in cells:
                    if not cell['all_evaluations_accepted']:
                        unknown.append(dict(preliminary=True,phase=phase,coordinates=coordinates,station=station,omega=omega,source=source,
                            source_row_id=rid,beta=cell['beta'],expected=cell['expected'],status=cell['status'],known_growth=cell['known_growth'],sample_ids=cell['sample_ids'],requested_grid_known=grid_known))
                if not cells:
                    unknown.append(dict(preliminary=True,phase=phase,coordinates=coordinates,station=station,omega=omega,source=source,
                        source_row_id=rid,beta=None,expected=False,status='beta_grid_unknown',known_growth=None,sample_ids=[],requested_grid_known=False))
                candidate=local_peaks(row,cells,by_id);row['beta_peak_candidate_ids']=[p['peak_id'] for p in candidate]
                rows.append(row);peaks.extend(candidate)
    joint_groups=defaultdict(list);joints=[]
    for row in rows:joint_groups[row['phase'],row['station'],key(row['omega'])].append(row)
    for (phase,station,_),source_rows in sorted(joint_groups.items(),key=lambda pair:(pair[0][0],pair[0][1],float(pair[0][2]))):
        omega=source_rows[0]['omega'];available=[r for r in source_rows if r['known_G'] is not None]
        winner=max(available,key=lambda r:r['known_G']) if available else None;growth=winner['known_G'] if winner else None
        expected_sources=source_rows[0].get('expected_sources',list(SOURCES))
        complete=len(source_rows)==len(expected_sources) and all(r['all_requested_cells_accepted'] for r in source_rows)
        joints.append(dict(preliminary=True,row_id=f'{phase}:{coordinates}:i{station:03d}:joint:w{key(omega)}',phase=phase,
            coordinates=coordinates,station=station,omega=omega,source_row_ids=[r['row_id'] for r in source_rows],expected_sources=expected_sources,
            eligible_sources=[r['source'] for r in available],unknown_sources=[r['source'] for r in source_rows if not r['all_requested_cells_accepted']],
            known_G=growth,beta_star=winner['beta_star'] if winner else None,winner_source=winner['source'] if winner else None,
            winner=winner['winner'] if winner else None,all_requested_cells_accepted=complete,sampled_grid_max_is_capped=complete,
            continuous_beta_max_certified=False,spectral_complete=False,evidence_state=state(growth,complete),
            unknown_cell_count=sum(r['unknown_beta_count'] for r in source_rows),beta_peak_candidate_ids=[p for r in source_rows for p in r['beta_peak_candidate_ids']]))
    return samples,rows,peaks,joints,unknown

def sign(v):return None if v is None else 1 if v>SIGN_TOL else -1 if v< -SIGN_TOL else 0

def omega_structure(rows,scope):
    groups=defaultdict(list)
    for row in rows:groups[row['phase'],row['station'],row.get('source','joint')].append(row)
    segments,brackets=[],[]
    for (phase,station,source),group in sorted(groups.items()):
        group.sort(key=lambda r:r['omega']);part=[]
        def flush():
            if part:
                segments.append(dict(preliminary=True,scope=scope,phase=phase,coordinates=part[0]['coordinates'],station=station,source=source,
                    evidence_state=part[0]['evidence_state'],omega_start=part[0]['omega'],omega_end=part[-1]['omega'],
                    omega_nodes=[r['omega'] for r in part],row_ids=[r['row_id'] for r in part],interpolation_assumed=False,
                    interpretation='consecutive sampled nodes only; no continuous interval is filled'))
        for row in group:
            if part and row['evidence_state']!=part[-1]['evidence_state']:flush();part=[]
            part.append(row)
        flush()
        def endpoint(row):return dict(row_id=row['row_id'],omega=row['omega'],known_G=row['known_G'],beta_star=row['beta_star'],
            winner=row['winner'],seed=row['winner']['seed'] if row['winner'] else None,evidence_state=row['evidence_state'],
            all_requested_cells_accepted=row['all_requested_cells_accepted'],beta_peak_candidate_ids=row['beta_peak_candidate_ids'])
        for row in group:
            if sign(row['known_G'])==0:
                brackets.append(dict(preliminary=True,bracket_id=row['row_id']+':zero',scope=scope,phase=phase,coordinates=row['coordinates'],
                    station=station,source=source,kind='near_zero_sample_candidate',left=endpoint(row),right=endpoint(row),
                    all_endpoint_requested_cells_accepted=row['all_requested_cells_accepted'],continuous_envelope_root_certified=False,refinement_required=True))
        for left,right in zip(group,group[1:]):
            a,b=sign(left['known_G']),sign(right['known_G'])
            if a is None or b is None or a*b>=0:continue
            complete=left['all_requested_cells_accepted'] and right['all_requested_cells_accepted']
            brackets.append(dict(preliminary=True,bracket_id=left['row_id']+':to:'+key(right['omega']),scope=scope,phase=phase,
                coordinates=left['coordinates'],station=station,source=source,kind='adjacent_observed_sign_change',left=endpoint(left),right=endpoint(right),
                all_endpoint_requested_cells_accepted=complete,negative_endpoint_not_capped=not (left if a<0 else right)['all_requested_cells_accepted'],
                winner_changed=left['winner']['source']!=right['winner']['source'],continuous_envelope_root_certified=False,refinement_required=True,
                warning='Discrete observed maxima only; unknown cells or unsampled beta competitors may remove this candidate bracket.'))
    return segments,brackets

def write_csv(path,records,omit=()):
    fields=list(dict.fromkeys(k for r in records for k in r if k not in omit)) or ['preliminary','schema_version']
    with path.open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');writer.writeheader()
        for r in records:
            writer.writerow({k:json.dumps(clean(v),separators=(',',':'),allow_nan=False) if isinstance(v,(dict,list,tuple)) else clean(v) for k,v in r.items() if k not in omit})

def validate(output):
    ids={v['sample_id'] for v in output['samples']}
    assert len(ids)==len(output['samples']),'sample IDs are not unique'
    for r in output['source_rows']+output['joint_rows']:
        assert not r['continuous_beta_max_certified'] and not r['spectral_complete']
        assert r['known_G'] is None or r['winner']['sample_id'] in ids
        if r['evidence_state'] in ('sampled_negative','sampled_zero'):assert r['all_requested_cells_accepted'],'unknown cells treated as stable'
        if not r['all_requested_cells_accepted']:assert not r['sampled_grid_max_is_capped']
    for p in output['beta_peak_candidates']:assert p['sample_id'] in ids and p['preliminary']
    assert all(v['coordinates']==output['coordinates'] for v in output['samples'])
    return dict(sample_ids_unique=True,candidate_samples_exist=True,coordinate_models_separated=True,
        unknown_cells_never_imply_stability=True,no_continuous_or_spectral_completeness_claim=True)

def self_check():
    # Two distinct interior beta maxima; a sampled stable hole and unknown hole.
    betas=[-2.,-1.,0.,1.,2.];omegas=[0.,1.,2.,3.,4.]
    shape=[0.,.5,0.,.6,-.2];level=[0.,-2.,0.,None,-2.];runs=[]
    for source in SOURCES:
        values=[]
        for w in omegas:
            for j,b in enumerate(betas):
                good=level[int(w)] is not None
                g=level[int(w)]+shape[j] if good else math.nan
                values.append(dict(preliminary=True,sample_id=f'test:{source}:{w}:{b}',phase='test',coordinates='global',
                    station=1,source=source,run_id='test:'+source,run_name=source,field_path='synthetic',field_line=len(values)+1,
                    omega_r=w,beta_r=b,alpha_r=10+b,alpha_i=-g,growth=g,status_code=0 if good else 2,
                    accepted_for_known_growth=good,exclusion_reasons=[] if good else ['synthetic failure']))
        runs.append(dict(phase='test',station=1,source=source,run_id='test:'+source,selected_coordinate_model=True,
            _omegas=omegas,_betas=betas,_samples=values))
    # An excluded-model run must contribute neither samples nor maxima.
    runs.append(dict(phase='test',station=1,source='M1p',run_id='excluded',selected_coordinate_model=False,
        _omegas=omegas,_betas=betas,_samples=[]))
    samples,rows,peaks,joints,unknown=make_rows(runs,'global')
    expected=['positive_exists','sampled_negative','positive_exists','unknown','sampled_negative']
    assert [r['evidence_state'] for r in joints]==expected
    segments,brackets=omega_structure(joints,'joint')
    assert [v['evidence_state'] for v in segments]==expected,'stable/unknown holes were filled'
    assert [(b['left']['omega'],b['right']['omega']) for b in brackets]==[(0.,1.),(1.,2.)],'bracket crossed unknown row'
    first=[v for v in peaks if v['source']=='M1p' and v['omega']==0]
    assert {v['beta'] for v in first}=={-2.,-1.,1.,2.},'a beta local maximum or endpoint was lost'
    assert len(samples)==150 and len(unknown)==30
    # A single failed competitor cannot turn an incomplete negative row stable.
    bad=runs[0]['_samples'][5];bad.update(accepted_for_known_growth=False,status_code=2,growth=math.nan)
    _,_,_,new_joint,_=make_rows(runs,'global')
    assert new_joint[1]['evidence_state']=='unknown_nonpositive_observed'
    assert not new_joint[1]['sampled_grid_max_is_capped']
    return dict(two_beta_peaks_and_endpoints_preserved=True,stable_hole_preserved=True,unknown_hole_preserved=True,
        no_bracket_across_unknown=True,failed_competitor_not_stable=True,excluded_model_omitted=True)

def main():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--self-test',action='store_true',help='Pure-data regression; no file reads or solver calls')
    p.add_argument('--root',type=Path,default=ROOT)
    p.add_argument('--phase',action='append',help='Exact phase; repeat to include separate phases. Default all *_progress.json phases.')
    p.add_argument('--coordinates',choices=('global','wall_local'),required=True)
    p.add_argument('--output-prefix',type=Path,help='Default outputs/envelope_reduce_PHASE_COORDINATES')
    args=p.parse_args()
    if args.self_test:print(json.dumps(self_check()));return 0
    root=args.root.resolve()
    phases=list(dict.fromkeys(args.phase or [f.name[:-len('_progress.json')] for f in sorted((root/'outputs').glob('*_progress.json'))]))
    if not phases or any(not re.fullmatch(r'[A-Za-z0-9_.-]+',v) for v in phases):p.error('provide a valid phase or outputs/*_progress.json')
    runs,manifests,issues=load_evidence(root,phases,args.coordinates)
    samples,rows,peaks,joints,unknown=make_rows(runs,args.coordinates)
    ss,sb=omega_structure(rows,'source');js,jb=omega_structure(joints,'joint')
    output=dict(schema_version=SCHEMA,preliminary=True,created=stamp(),root=str(root),phases=phases,coordinates=args.coordinates,
        expected_sources_by_phase={m['phase']:m['expected_sources'] for m in manifests},inputs=manifests,
        limitations=['Declared search starts are provenance labels, not complete spectral enumeration.',
            'known_G is a lower bound from eligible actual samples, never a certified continuous-beta upper bound.',
            'Positive eligible growth is numerical existence evidence; missing/failed or nonpositive-incomplete cells are not stable evidence.',
            'All_requested_cells_accepted describes only discrete requested cells, not unsampled beta/frequency regions.',
            'Omega segments retain separate sampled-node sign and unknown holes; no continuous interval is filled.',
            'Partial-run successful samples require a readable summary and no global PETSc/options/runtime/model integrity issue.',
            'All raw success/failure samples, partial runs and excluded-model metadata are retained; inputs are never modified.'],
        tolerances=dict(parameter_match_relative=PARAM_TOL,growth_sign=SIGN_TOL,peak_comparison=PEAK_TOL),runs=runs,input_issues=issues,
        samples=samples,source_rows=rows,beta_peak_candidates=peaks,joint_rows=joints,omega_segments=js+ss,neutral_brackets=jb+sb,unknown_cells=unknown)
    output['validation']=validate(output)
    output['summary']=dict(run_slots=len(runs),selected_model_runs=sum(r['selected_coordinate_model'] for r in runs),
        excluded_model_runs=sum(not r['selected_coordinate_model'] for r in runs),actual_samples=len(samples),
        eligible_samples=sum(v['accepted_for_known_growth'] for v in samples),source_rows=len(rows),joint_rows=len(joints),
        beta_peak_candidates=len(peaks),unknown_cells=len(unknown),positive_joint_nodes=sum(r['evidence_state']=='positive_exists' for r in joints),
        sampled_negative_joint_nodes=sum(r['evidence_state']=='sampled_negative' for r in joints),
        unknown_joint_nodes=sum(r['evidence_state'].startswith('unknown') for r in joints),neutral_candidates=len(jb+sb),formal_complete=False)
    prefix=args.output_prefix or root/'outputs'/('envelope_reduce_'+'_'.join(phases)+'_'+args.coordinates)
    if not prefix.is_absolute():prefix=root/prefix
    prefix.parent.mkdir(parents=True,exist_ok=True);path=prefix.with_suffix('.json');tmp=path.with_suffix('.json.tmp')
    tmp.write_text(json.dumps(clean(output),indent=2,allow_nan=False)+'\n');tmp.replace(path)
    tables=(('samples',samples,()),('source_rows',rows,('cells',)),('beta_peaks',peaks,()),('joint',joints,()),
        ('omega_segments',js+ss,()),('neutral_brackets',jb+sb,()),('unknown_cells',unknown,()),('runs',runs,()))
    for suffix,records,omit in tables:write_csv(prefix.parent/(prefix.name+'_'+suffix+'.csv'),clean(records),omit)
    print(json.dumps(dict(json=str(path),**output['summary'],validation=output['validation']),allow_nan=False),flush=True)
    return 0

if __name__=='__main__':raise SystemExit(main())



