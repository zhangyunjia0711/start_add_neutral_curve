#!/usr/bin/env python3
"""Actual local-wall-coordinate seeds on the main 841-station grid.

Each direction starts at a freshly solved seed and traverses in physical order.
An unresolved station stops that direction; later samples are never propagated
across it. Family names describe initial-guess provenance, not modal identities.
"""
import argparse
import concurrent.futures
import configparser
import csv
import copy
import datetime
import json
import math
from pathlib import Path
import re
import threading
import time

import full_driver as driver

ROOT = driver.ROOT
OUTPUT = ROOT / 'outputs/axial_wall_local.json'
OPTIONS = ['-pep_tol', '1e-12', '-pep_max_it', '1000', '-pep_nev', '1',
           '-pep_ncv', '8', '-st_pc_factor_mat_solver_type', 'superlu']
LOCK = threading.Lock()
STATE = {}
ARGS = None
X = []
EPS = 1e-8
NUMBER = r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?'
PAIR = re.compile(r'^\s*\d+\s+\(\s*(' + NUMBER + r')\s*,\s*(' + NUMBER
                  + r')\s*\)\s+The error is\s+(' + NUMBER + r')', re.M)


def stamp():
    return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()


def seed_dict(value):
    return {key: value[key] for key in ['omega_r', 'beta_r', 'alpha_r', 'alpha_i']}


def initial_guesses(include_nose_generic=False,include_nearwall=False):
    catalogue = json.loads((ROOT / 'outputs/seed_catalog.json').read_text())
    baseline = {v['mode']: v for v in catalogue['baseline_converged_seeds']}
    result = {}
    for family, station in [('M1', 317), ('M2', 312)]:
        value = copy.deepcopy(baseline[family])
        value.update(family=family, iloc=station,
                     guess_kind='archived global-model flat-profile guess; fresh local solve required')
        result[family] = value
    local = json.loads((ROOT / 'outputs/wall_local_actual_seeds_v5.json').read_text())
    candidates = sorted([v for v in local if v['iloc'] == 119 and v['status'] == 0],
                        key=lambda v: v['alpha_r'])
    if len(candidates) != 2:
        raise RuntimeError('Exactly two separately documented i119 local-model roots are required')
    for suffix, value in zip(['a', 'b'], candidates):
        result['M3_local_' + suffix] = dict(copy.deepcopy(value), family='M3_local_' + suffix,
            guess_kind='actual v5 local-wall root; freshly re-solve with immutable v6')
    if include_nose_generic:
        value = next(v for v in local if v['iloc'] == 1 and v['status'] == 0)
        result['nose_generic'] = dict(copy.deepcopy(value), family='nose_generic',
            guess_kind='actual v5 i1 generic-shift local root; independent unclassified search start')
    if include_nearwall:
        value=json.loads((ROOT/'outputs/local_mode_audit_nearwall_seed.json').read_text())
        assert value['family']=='M3near' and value['coordinates']=='wall_local' and value['status']==0
        result['M3near']=dict(copy.deepcopy(value),
            guess_kind='independently audited actual local near-wall candidate; fresh solve required')
    return result, catalogue['grids']['main_841']['x_nondim']


def persist(family, local, records):
    with LOCK:
        existing = {(v['family'], v['iloc']): v for v in STATE['records']}
        for value in records:
            key = value['family'], value['iloc']
            if key not in existing or existing[key]['status'] != 0 or value['status'] == 0:
                existing[key] = value
        STATE['records'] = sorted(existing.values(), key=lambda v: (v['family'], v['iloc']))
        STATE['families'][family] = copy.deepcopy(local)
        STATE['updated'] = stamp()
        STATE['coverage'] = {}
        for label in STATE['initial_guesses']:
            values = [v for v in STATE['records'] if v['family'] == label]
            accepted = {v['iloc'] for v in values if v['status'] == 0}
            STATE['coverage'][label] = dict(accepted=len(accepted),
                failed_or_unverified=[v['iloc'] for v in values if v['status'] != 0],
                not_yet_covered=[i for i in range(1, 842) if i not in accepted])
        temporary = OUTPUT.with_suffix('.json.tmp')
        temporary.write_text(json.dumps(STATE, indent=2, allow_nan=False) + '\n')
        temporary.replace(OUTPUT)


def export_coverage():
    lookup={(v['family'],v['iloc']):v for v in STATE['records']}
    columns=['family','iloc','x_mm','coordinates','status','alpha_r','alpha_i','growth',
             'omega_r','beta_r','source_path','source_solver_sha256','reason']
    target=OUTPUT.with_name(OUTPUT.stem+'_coverage.csv')
    with target.open('w',newline='') as out:
        writer=csv.DictWriter(out,fieldnames=columns);writer.writeheader()
        for family in STATE['initial_guesses']:
            for station in range(1,842):
                value=lookup.get((family,station))
                row={key:value.get(key) if value else None for key in columns}
                row.update(family=family,iloc=station,x_mm=35*X[station-1],coordinates='wall_local',
                    status=value['status'] if value else 'unknown',
                    reason='; '.join(value['reasons']) if value else 'no accepted actual same-station axial seed')
                writer.writerow(row)
    return str(target)


def execute_block(family, lo, hi, anchor, direction, jump_limit=2, prefix='axial_local'):
    name = f'{prefix}_{family}_{direction}_i{lo:03d}_{hi:03d}_{ARGS.solver_sha[:8]}'
    spec = dict(name=name, solver=str(ARGS.solver), solver_sha256=ARGS.solver_sha,
        coordinates='wall_local', mode='M1', family=family, station=lo,
        seed_station=anchor['iloc'], seed=seed_dict(anchor), seed_shift=1e-4,
        task=33, outer_parameter='iloc', inner_parameter='omega',
        outer_min=lo, outer_max=hi, outer_points=hi-lo+1,
        inner_min=anchor['nominal_omega'], inner_max=anchor['nominal_omega']+EPS,
        inner_points=2, reverse_outer=direction == 'up', anchor_row_seed=True,
        sample_only=True, adaptive=True, vector_tracking=True, minimum_overlap=.95,
        minimum_step=1e-5, retry_eigenpairs=8, step=.4, jump=jump_limit,
        options=OPTIONS, timeout=1200)
    if ARGS.station_predictor:spec['station_predictor']=True
    summary = driver.run(spec)
    directory = ROOT / 'runs' / name
    actual_config = configparser.ConfigParser()
    actual_config.read(directory / 'config.ini')
    neutral = actual_config['neutral']
    istrue = lambda v: v.strip().lower() in ('t','true','yes','1')
    configuration_ok = all(istrue(neutral.get(k,'')) == bool(spec[k]) for k in
        ['sample_only','vector_tracking','anchor_row_seed','reverse_outer'])
    configuration_ok = configuration_ok and istrue(neutral.get('adaptive_steps',''))
    configuration_ok = configuration_ok and int(neutral.get('retry_eigenpairs','0')) == 8
    configuration_ok = configuration_ok and abs(float(neutral.get('minimum_overlap','0'))-.95) < 1e-12
    if ARGS.station_predictor:configuration_ok=configuration_ok and istrue(neutral.get('station_predictor',''))
    field = driver.table(directory / 'out/neutral_field.plt')
    text = (directory / 'solver.log').read_text(errors='replace')
    pairs = [(complex(float(m[1].replace('D', 'E')), float(m[2].replace('D', 'E'))),
              float(m[3].replace('D', 'E'))) for m in PAIR.finditer(text)]
    fatal = bool(not configuration_ok or summary.get('petsc_errors') or summary.get('unused_options')
                 or summary.get('coordinate_metadata_missing') or summary['exit'] == 124
                 or re.search(r'SIGSEGV|segmentation fault|Fortran runtime error', text, re.I))
    records = []
    traversal = list(range(lo, hi+1))
    if direction == 'up': traversal.reverse()
    previous, first_gap = anchor, None
    for station in traversal:
        values = [v for v in field if int(round(v['iloc'])) == station]
        good = [v for v in values if v['status'] == 0 and all(math.isfinite(v[k]) for k in
            ['x', 'omega_r', 'beta_r', 'alpha_r', 'alpha_i', 'growth'])
            and abs(v['omega_r']-anchor['nominal_omega']) <= EPS+2e-12
            and abs(v['beta_r']-anchor['beta_r']) <= 1e-10]
        chosen = min(good, key=lambda v: abs(v['omega_r']-anchor['nominal_omega'])) if good else None
        reasons = []
        if len(values) != 2: reasons.append('expected two actual tiny-frequency samples')
        if any(v['status'] != 0 for v in values): reasons.append('failed sample at station')
        if chosen is None: reasons.append('no valid local-wall fixed-parameter sample')
        if fatal: reasons.append('block contains PETSc/options/runtime/coordinate-metadata error')
        matching, residual = None, None
        if chosen is not None:
            alpha = complex(chosen['alpha_r'], chosen['alpha_i'])
            if pairs:
                matching, residual = min(((abs(z-alpha), err) for z, err in pairs), key=lambda p:p[0])
            if matching is None or matching > 1e-8: reasons.append('no matching printed PEP residual')
            elif residual > 1e-7: reasons.append('selected original-polynomial residual exceeds 1e-7')
            if abs(chosen['x']-X[station-1]) > 1e-10: reasons.append('x differs from main grid')
        value = dict(family=family, coordinates='wall_local', iloc=station, x=X[station-1],
            x_mm=35*X[station-1], nominal_omega=anchor['nominal_omega'],
            omega_r=chosen['omega_r'] if chosen else anchor['nominal_omega'],
            beta_r=anchor['beta_r'], alpha_r=chosen['alpha_r'] if chosen else None,
            alpha_i=chosen['alpha_i'] if chosen else None,
            growth=chosen['growth'] if chosen else None, status=0 if not reasons else 2,
            reasons=reasons, direction=direction, source_path=str(directory/'out/neutral_field.plt'),
            source_row=chosen['row'] if chosen else None, source_solver_sha256=ARGS.solver_sha,
            printed_pep_error=residual, printed_pep_matching_distance=matching,
            template='M1 main_841', physical_mode_identity='unclassified continuation of named search start')
        if first_gap is not None:
            value['status'] = 3
            value['reasons'].append('beyond unresolved station; no family continuity assumed')
        elif value['status'] == 0:
            if previous['iloc'] != station:
                jump = abs(complex(value['alpha_r']-previous['alpha_r'], value['alpha_i']-previous['alpha_i']))
                value['alpha_jump_from_adjacent'] = jump
                if abs(previous['iloc']-station) != 1 or jump > jump_limit:
                    value['status'] = 3
                    value['reasons'].append('nonadjacent anchor or alpha jump exceeds '+str(jump_limit))
            if value['status'] == 0: previous = value
        if value['status'] != 0 and first_gap is None: first_gap = station
        records.append(value)
    result = dict(name=name, range=[lo,hi], direction=direction, terminal=summary['terminal'],
        exit=summary['exit'], seconds=summary['seconds'], solves=summary.get('solves'),
        rejected=summary.get('rejected'), failed_field_records=summary.get('failed_field_records'),
        first_unresolved_station=first_gap, accepted_stations=[v['iloc'] for v in records if v['status']==0],
        summary_path=str(directory/'summary.json'))
    result['actual_configuration_verified'] = configuration_ok
    overlaps=[(float(m[1].replace('D','E')), m[2].upper()=='T') for m in re.finditer(
        r'vector_overlap=\s*('+NUMBER+r')\s+accepted=\s*([TF])',
        (directory/'out/neutral_run.log').read_text())]
    result['overlap_checks']=[dict(overlap=o,accepted=accepted) for o,accepted in overlaps]
    result['jump_limit']=jump_limit
    return result, records, previous


def gap_probes():
    """Independent larger-jump diagnostics; never rewrite the axial seed JSON."""
    previous=json.loads(OUTPUT.read_text())
    assert previous['solver_sha256']==ARGS.solver_sha
    planned=[]
    selected=ARGS.families or ['M1','M2']
    locations={'M1':(56,55),'M2':(20,19),'M3near':(56,55)}
    for family in selected:
        old,new=locations[family]
        anchor=next(v for v in previous['records'] if v['family']==family and v['iloc']==old and v['status']==0)
        for jump_limit in [4,8]:planned.append((family,old,new,jump_limit,anchor))
    results=[]
    def probe(item):
        family,old,new,jump_limit,anchor=item
        prefix=f'gapcheck_jump{jump_limit}'
        block,records,_=execute_block(family,new,new,anchor,'out',jump_limit,prefix)
        result=dict(family=family,from_station=old,to_station=new,jump_limit=jump_limit,
                    original_anchor=anchor,forward_block=block,forward_records=records)
        actual=records[0]
        if actual['status']==0:
            result['actual_delta_alpha']=abs(complex(actual['alpha_r']-anchor['alpha_r'],actual['alpha_i']-anchor['alpha_i']))
            back,back_records,_=execute_block(family,old,old,actual,'back',jump_limit,prefix)
            result.update(back_block=back,back_records=back_records)
            if back_records[0]['status']==0:
                returned=back_records[0]
                result['roundtrip_alpha_difference']=abs(complex(returned['alpha_r']-anchor['alpha_r'],returned['alpha_i']-anchor['alpha_i']))
        return result
    with concurrent.futures.ThreadPoolExecutor(max_workers=min(2,ARGS.workers)) as pool:
        for future in concurrent.futures.as_completed([pool.submit(probe,item) for item in planned]):
            results.append(future.result())
    output=dict(created=stamp(),solver=str(ARGS.solver),solver_sha256=ARGS.solver_sha,
        qualification='Independent jump-limit diagnosis; original axial JSON is unchanged',
        minimum_overlap=.95,eigen_error_tol=1e-7,results=results)
    version=re.search(r'_v(\d+)',ARGS.solver.name)[1]
    target=ROOT/f'outputs/axial_wall_local_gap_probes_{"_".join(selected)}_v{version}.json'
    target.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print('GAP PROBES '+json.dumps(dict(output=str(target),results=[{k:v.get(k) for k in
        ['family','jump_limit','actual_delta_alpha','roundtrip_alpha_difference']} for v in results])),flush=True)


def family_worker(family, guess):
    with LOCK: local = copy.deepcopy(STATE['families'].get(family, {}))
    if 'seed' not in local:
        anchor = dict(guess, nominal_omega=guess['omega_r'])
        block, records, _ = execute_block(family, guess['iloc'], guess['iloc'], anchor, 'seed',ARGS.jump,
                                         'axial_local_j'+str(ARGS.jump))
        seed = records[0]
        if seed['status'] == 0:
            distance = abs(complex(seed['alpha_r']-guess['alpha_r'],seed['alpha_i']-guess['alpha_i']))
            seed['bootstrap_guess_alpha_distance'] = distance
            if family.startswith('M3') and distance > 1e-3:
                seed['status']=3; seed['reasons'].append('fresh local root differs from recorded local bootstrap')
            if family == 'M2' and distance > 1:
                seed['status']=3; seed['reasons'].append('mapped M2 initial root requires separate review')
        local = dict(seed=seed, guess_provenance=guess, blocks=[block], directions={})
        persist(family, local, [seed])
        print('AXIAL LOCAL SEED '+json.dumps(seed), flush=True)
    if local['seed']['status'] != 0 or ARGS.phase == 'bootstrap': return
    for direction, sign in [('down',1),('up',-1)]:
        state = local['directions'].setdefault(direction, dict(anchor=copy.deepcopy(local['seed']),
            next_station=local['seed']['iloc']+sign, stopped=False, completed_blocks=0))
        while not state['stopped'] and 1 <= state['next_station'] <= 841:
            if ARGS.phase == 'diagnostic' and state['completed_blocks'] >= 1: break
            start = state['next_station']
            lo,hi = ((start,min(841,start+ARGS.block_size-1)) if sign==1 else
                     (max(1,start-ARGS.block_size+1),start))
            block, records, endpoint = execute_block(family,lo,hi,state['anchor'],direction,ARGS.jump,
                                                      'axial_local_j'+str(ARGS.jump))
            local['blocks'].append(block); state['completed_blocks'] += 1
            state['anchor'] = endpoint
            if block['first_unresolved_station'] is not None:
                state['stopped']=True
                state['stop_reason']='unresolved station; no propagation across gap'
            else: state['next_station'] = hi+1 if sign==1 else lo-1
            persist(family,local,records)
            print('AXIAL LOCAL BLOCK '+json.dumps(dict(family=family,**block)), flush=True)
        state['finished_native_boundary'] = not 1 <= state['next_station'] <= 841
    persist(family,local,[])


def main():
    global ARGS,STATE,X,OUTPUT
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver',type=Path,required=True);parser.add_argument('--solver-sha',required=True)
    parser.add_argument('--phase',choices=['bootstrap','diagnostic','full','gap-probes'],default='diagnostic')
    parser.add_argument('--workers',type=int,default=4);parser.add_argument('--block-size',type=int,default=20)
    parser.add_argument('--include-nose-generic',action='store_true');parser.add_argument('--families',nargs='*')
    parser.add_argument('--include-nearwall',action='store_true')
    parser.add_argument('--output',type=Path,default=OUTPUT)
    parser.add_argument('--resume-from',type=Path,
        help='Clone a historical actual-seed state into a NEW output; retain per-record solver provenance')
    parser.add_argument('--jump',type=float,default=2)
    parser.add_argument('--station-predictor',action='store_true',help='Reviewed bounded station prediction, v10 or later')
    ARGS=parser.parse_args();ARGS.solver=ARGS.solver.resolve()
    assert ARGS.solver.is_file() and driver.sha(ARGS.solver)==ARGS.solver_sha
    version=re.search(r'_v(\d+)',ARGS.solver.name)
    assert version and int(version[1])>=6,'New traversal API requires reviewed v6-or-later build'
    if ARGS.station_predictor:assert int(version[1])>=10
    OUTPUT=ARGS.output.resolve()
    assert OUTPUT.parent==ROOT/'outputs','Output must stay within this task outputs directory'
    assert ARGS.jump>0
    driver.snapshot();guesses,X=initial_guesses(ARGS.include_nose_generic,ARGS.include_nearwall)
    if ARGS.phase=='gap-probes':
        gap_probes();driver.snapshot();return
    if ARGS.resume_from and not OUTPUT.exists():
        inherited=json.loads(ARGS.resume_from.read_text())
        STATE=copy.deepcopy(inherited)
        STATE.update(created=stamp(),solver=str(ARGS.solver),solver_sha256=ARGS.solver_sha,
            inherited_seed_file=str(ARGS.resume_from.resolve()),
            inherited_seed_file_sha256=driver.sha(ARGS.resume_from),
            inherited_executor_solver_sha256=inherited['solver_sha256'],
            inheritance_qualification='Historical actual records retain their source_solver_sha256; current executor freshly solves every new continuation point',
            active_jump_limit=ARGS.jump)
        assert STATE['initial_guesses']==guesses
        for family,local in STATE['families'].items():
            for direction,state in local['directions'].items():
                if state['stopped'] and state['anchor']['status']==0:
                    state['historical_stop_reason']=state.get('stop_reason')
                    state['stopped']=False
                    state['next_station']=state['anchor']['iloc']+(1 if direction=='down' else -1)
                    state['retry_under_jump_limit']=ARGS.jump
        OUTPUT.write_text(json.dumps(STATE,indent=2,allow_nan=False)+'\n')
    if OUTPUT.exists():
        STATE=json.loads(OUTPUT.read_text());assert STATE['solver_sha256']==ARGS.solver_sha
        assert STATE['initial_guesses']==guesses,'initial search-start set changed; use distinct output'
    else:
        STATE=dict(created=stamp(),solver=str(ARGS.solver),solver_sha256=ARGS.solver_sha,
            initial_guesses=guesses,records=[],families={},coordinates='wall_local',
            grid='M1 main_841',tracking=dict(minimum_overlap=.95,retry_eigenpairs=8,normal_nev=1),
            qualification='Finite named search starts; no whole-spectrum or physical mode identification claim')
    families=ARGS.families or list(guesses)
    assert all(f in guesses for f in families)
    with concurrent.futures.ThreadPoolExecutor(max_workers=ARGS.workers) as pool:
        jobs=[pool.submit(family_worker,f,guesses[f]) for f in families]
        for job in concurrent.futures.as_completed(jobs):job.result()
    driver.snapshot()
    print('AXIAL LOCAL COVERAGE CSV '+export_coverage(),flush=True)
    print('AXIAL LOCAL COVERAGE '+json.dumps(STATE['coverage']),flush=True)


if __name__=='__main__':main()
