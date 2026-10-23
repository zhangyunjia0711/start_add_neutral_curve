#!/usr/bin/env python3
"""Bounded adjacent-station and return verification for the M3near nose start."""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace
import axial_wall_local as axial

ROOT = axial.ROOT
SOLVER = ROOT/'outputs/bin/START_LST_wall_local_v10'
SHA = '634b45c5d30e076f5b86ccb136e0734c1905beedde108477dceb8a75adc3853c'

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--station', type=int, required=True)
    p.add_argument('--jumps', type=float, nargs='+', required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    assert not args.output.exists(), 'Do not overwrite prior evidence'
    assert args.output.parent.resolve() == ROOT/'outputs'
    axial.ARGS = SimpleNamespace(solver=SOLVER, solver_sha=SHA, station_predictor=True)
    _, axial.X = axial.initial_guesses(include_nearwall=True)
    assert axial.driver.sha(SOLVER) == SHA
    source_hash = axial.driver.sha(args.source)
    prior = json.loads(args.source.read_text())
    anchor = next(v for v in prior['records'] if v['family']=='M3near'
                  and v['iloc']==args.station and v['status']==0)
    result = dict(created=axial.stamp(), source=str(args.source),source_sha256=source_hash,
        solver=str(SOLVER), solver_sha256=SHA, original_anchor=anchor,
        minimum_overlap=.95, residual_limit=1e-7, roundtrip_limit=1e-7,
        qualification='Fixed omega0 beta-7 adjacent continuation; physical modal identity unclassified',
        probes=[])
    def save():
        args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    axial.driver.snapshot()
    for jump in args.jumps:
        prefix=f'nose_near_probe_from{args.station}_j{jump:g}'
        block, records, _ = axial.execute_block('M3near',args.station-1,args.station-1,
                                               anchor,'out',jump,prefix)
        item = dict(jump_limit=jump, forward_block=block, forward_records=records)
        actual=records[0]
        if actual['status']==0:
            item['delta_alpha']=abs(complex(actual['alpha_r']-anchor['alpha_r'],
                                            actual['alpha_i']-anchor['alpha_i']))
            back, returned, _=axial.execute_block('M3near',args.station,args.station,
                                                 actual,'back',jump,prefix)
            item.update(back_block=back,back_records=returned)
            if returned[0]['status']==0:
                item['roundtrip_alpha_difference']=abs(complex(returned[0]['alpha_r']-anchor['alpha_r'],
                                                          returned[0]['alpha_i']-anchor['alpha_i']))
        checks=[q for k in ['forward_block','back_block'] for q in item.get(k,{}).get('overlap_checks',[])]
        item['min_actual_mac']=min([q['overlap'] for q in checks],default=None)
        item['all_actual_mac_checks_accepted']=bool(checks) and all(q['accepted'] for q in checks)
        item['passed']=bool(actual['status']==0 and item.get('back_records',[{}])[0].get('status')==0
            and item['all_actual_mac_checks_accepted'] and item['min_actual_mac']>=.95
            and item['roundtrip_alpha_difference']<=1e-7)
        result['probes'].append(item)
        save()
        print('NOSE PROBE '+json.dumps({k:item.get(k) for k in
            ['jump_limit','delta_alpha','roundtrip_alpha_difference','min_actual_mac','passed']}),flush=True)
    passed=[v['forward_records'][0] for v in result['probes'] if v['passed']]
    result['cross_jump_root_spread']=max((abs(complex(v['alpha_r']-w['alpha_r'],v['alpha_i']-w['alpha_i']))
        for v in passed for w in passed),default=None)
    result['all_requested_probes_passed_same_root']=bool(len(passed)==len(args.jumps)
        and result['cross_jump_root_spread']<=1e-7)
    assert axial.driver.sha(args.source)==source_hash
    assert axial.driver.sha(SOLVER)==SHA
    axial.driver.snapshot()
    result['completed']=axial.stamp();save()
    print('NOSE PROBE RESULT '+str(args.output),flush=True)

if __name__=='__main__':main()
