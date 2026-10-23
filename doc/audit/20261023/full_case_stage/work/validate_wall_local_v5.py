#!/usr/bin/env python3
"""Actual frozen-binary checks for the optional wall-local spatial LST profiles.

Run only after the coordinating agent freezes v5. This script never builds or
changes source. Original-global nose guesses are explicitly guesses; no modal
identity is inferred from an initial guess or a finite eigenvalue.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import full_driver as driver
from seed_catalog import load, get_candidates

ROOT = Path(__file__).resolve().parent.parent


def metadata(directory, station, catalogue):
    path = directory/'wall_local_profile_metadata.csv'
    if not path.exists():
        return dict(valid=False, issue='missing metadata', source_path=str(path))
    rows=list(csv.DictReader(line for line in path.read_text().splitlines() if not line.startswith('#')))
    selected=[row for row in rows if int(row['iloc'])==station]
    if len(selected)!=1:
        return dict(valid=False, issue='expected one metadata station', source_path=str(path), rows=rows)
    row=selected[0]
    numeric={key:float(value) for key,value in row.items()
             if key not in ('iloc','identity_flat','normal_derivatives_finite')}
    expected_x=catalogue['grids']['main_841']['x_nondim'][station-1]
    t=[numeric['t_x'],numeric['t_y']];n=[numeric['n_x'],numeric['n_y']]
    valid=all(math.isfinite(v) for v in numeric.values())
    valid &= abs(numeric['x_global_N']-expected_x)<1e-12
    valid &= abs(sum(v*v for v in t)-1)<1e-12
    valid &= abs(sum(v*v for v in n)-1)<1e-12
    valid &= abs(sum(a*b for a,b in zip(t,n)))<1e-12
    valid &= numeric['min_delta_n']>0 and numeric['n_min']==0
    valid &= row['normal_derivatives_finite']=='T'
    return dict(valid=bool(valid),source_path=str(path),values=row)


def execute(solver, sha, name, station, seed, coordinates):
    spec=dict(name=name,station=station,mode='M1',coordinates=coordinates,
              seed={key:seed[key] for key in ('omega_r','beta_r','alpha_r','alpha_i')},
              solver=str(solver),solver_sha256=sha,task=33,sample_only=True,
              omega_min=seed['omega_r'],omega_max=seed['omega_r'],omega_points=1,
              beta_min=seed['beta_r']-1e-4,beta_max=seed['beta_r']+1e-4,beta_points=3,
              step=.02,seed_shift=1e-3,eigen_tol=1e-7,jump=2,timeout=180,
              guess_note=seed.get('guess_note','archived initial guess; no exact modal mapping assumed'),
              options=['-pep_tol','1e-12','-pep_max_it','1000','-pep_nev','8','-pep_ncv','24',
                       '-st_pc_factor_mat_solver_type','superlu'])
    summary=driver.run(spec)
    directory=ROOT/'runs'/name
    field=[row for row in driver.table(directory/'out/neutral_field.plt') if row['status']==0]
    center=[row for row in field if abs(row['beta_r']-seed['beta_r'])<1e-10]
    usable=bool(center) and summary['petsc_errors']==0 and summary['unused_options']==0
    usable &= bool(center) and all(math.isfinite(value) for value in center[0].values())
    return dict(name=name,summary=summary,field=field,center=center[0] if center else None,
                finite_center_without_solver_errors=usable)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver',type=Path,required=True)
    parser.add_argument('--tag',default='v5')
    args=parser.parse_args()
    solver=args.solver.resolve()
    expected_sha=driver.sha(solver)
    catalogue=load()
    report=dict(solver=str(solver),solver_sha256=expected_sha,flat=[],nose=[],symmetry=[],
                notes=['Actual solver checks; no new builds.',
                       'Finite local roots alone do not establish a physical modal identity.',
                       'New wall-local nose guesses are not claimed as exact mapped global roots.'])
    out=ROOT/'outputs'/('wall_local_validation_'+args.tag+'.json')
    def save():
        out.write_text(json.dumps(report,indent=2)+'\n')
        seeds=[]
        for station_report in report['nose']:
            for attempt in station_report['attempts']:
                if attempt['finite_center_without_solver_errors']:
                    source=str(ROOT/'runs'/attempt['name']/'out/neutral_field.plt')
                    seeds.append(dict(attempt['center'],iloc=int(attempt['center']['iloc']),status=0,
                                      coordinates='wall_local',family='M3',
                                      family_role='bootstrap/probe source label only, not physical mode identification',
                                      sourcepath=source,source_path=source,source_solver_sha256=expected_sha,
                                      guess_note=attempt['summary']['spec']['guess_note'],
                                      physical_mode_identity='unclassified new local-model root'))
        (ROOT/'outputs'/('wall_local_actual_seeds_'+args.tag+'.json')).write_text(json.dumps(seeds,indent=2)+'\n')

    for station in (242,317,841):
        seed=get_candidates(catalogue,station,'fixed_M1',omega=8.17,beta=-.98903,limit=1)[0]
        runs=[execute(solver,expected_sha,f'check_{args.tag}_flat_i{station}_{coordinates}',
                      station,seed,coordinates) for coordinates in ('global','wall_local')]
        a,b=runs
        check=metadata(ROOT/'runs'/b['name'],station,catalogue)
        comparison=dict(station=station,runs=runs,metadata=check,passed=False)
        if a['finite_center_without_solver_errors'] and b['finite_center_without_solver_errors']:
            ar=a['center'];br=b['center']
            difference=abs(complex(ar['alpha_r'],ar['alpha_i'])-complex(br['alpha_r'],br['alpha_i']))
            comparison.update(alpha_difference=difference,
                              passed=bool(difference<1e-8 and check['valid'] and
                                          check['values']['identity_flat']=='T'))
        report['flat'].append(comparison);save()

    accepted_nose={}
    for station in (119,1,2,21,61,241):
        # Local-curved archive guess is merely an initial guess for the new
        # zero-curvature model. A separate generic guess is retained on failure.
        attempts=[]
        guesses=[dict(omega_r=0.0,beta_r=-7.0,alpha_r=.5,alpha_i=-.5,
                      guess_note='generic small local alpha initial guess'),
                 dict(omega_r=-2.387183,beta_r=-6.6972,alpha_r=16.42,alpha_i=-2.24,
                      guess_note='curvature-model archived value used only as a zero-curvature local initial guess'),
                 dict(omega_r=0.0,beta_r=-7.0,alpha_r=5.0,alpha_i=-.5,
                      guess_note='generic larger local alpha initial guess')]
        if station==119:
            old=get_candidates(catalogue,119,'fixed_M3',omega=0.0,beta=-7.0,limit=1)[0]
            xwall=catalogue['grids']['main_841']['x_nondim'][118]
            t_x=math.sqrt(1-xwall*xwall)
            projected=dict(omega_r=old['omega_r'],beta_r=old['beta_r'],
                           alpha_r=t_x*old['alpha_r'],alpha_i=t_x*old['alpha_i'],
                           guess_note='t_x times old global alpha is a geometric guess, not an exact modal transform')
            unchanged=dict(omega_r=old['omega_r'],beta_r=old['beta_r'],
                           alpha_r=old['alpha_r'],alpha_i=old['alpha_i'],
                           guess_note='unchanged old global alpha used only as an initial guess in new local model')
            guesses=[projected,unchanged]+guesses
        for number,seed in enumerate(guesses,1):
            trial=execute(solver,expected_sha,f'check_{args.tag}_nose_i{station}_guess{number}',
                          station,seed,'wall_local')
            trial['metadata']=metadata(ROOT/'runs'/trial['name'],station,catalogue)
            attempts.append(trial)
            if trial['finite_center_without_solver_errors']:
                if station not in accepted_nose:accepted_nose[station]=trial
                if station!=119 or number>=2:break
        report['nose'].append(dict(station=station,attempts=attempts,
                                   geometry_passed=all(a['metadata']['valid'] for a in attempts),
                                   finite_root_found=station in accepted_nose))
        save()

    original=[report['flat'][1]['runs'][1]]
    if 119 in accepted_nose:original.append(accepted_nose[119])
    for first in original:
        if not first['finite_center_without_solver_errors']:continue
        row=first['center'];station=int(row['iloc'])
        seed=dict(omega_r=-row['omega_r'],beta_r=-row['beta_r'],
                  alpha_r=-row['alpha_r'],alpha_i=row['alpha_i'])
        second=execute(solver,expected_sha,f'check_{args.tag}_symmetry_i{station}',station,seed,'wall_local')
        item=dict(station=station,original_run=first['name'],transformed_run=second,passed=False)
        if second['finite_center_without_solver_errors']:
            other=second['center']
            error=abs(complex(other['alpha_r'],other['alpha_i'])-
                      complex(-row['alpha_r'],row['alpha_i']))
            item.update(alpha_symmetry_error=error,growth_difference=abs(other['growth']-row['growth']),
                        passed=bool(error<1e-7))
        report['symmetry'].append(item);save()
    report['solver_unchanged']=driver.sha(solver)==expected_sha
    report['passed']=bool(report['solver_unchanged'] and all(x['passed'] for x in report['flat']) and
                          all(x['geometry_passed'] and x['finite_root_found'] for x in report['nose']) and
                          len(report['symmetry'])==2 and all(x['passed'] for x in report['symmetry']))
    save()
    print(json.dumps(dict(report=str(out),passed=report['passed'],
                          flat_pass=[x['passed'] for x in report['flat']],
                          nose_geometry_pass=[x['geometry_passed'] for x in report['nose']],
                          nose_root_found=[x['finite_root_found'] for x in report['nose']],
                          symmetry_pass=[x['passed'] for x in report['symmetry']]),indent=2),flush=True)
    raise SystemExit(0 if report['passed'] else 1)


if __name__=='__main__':main()
