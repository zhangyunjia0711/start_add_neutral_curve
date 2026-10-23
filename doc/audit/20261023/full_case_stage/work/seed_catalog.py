#!/usr/bin/env python3
"""Build a provenance-preserving catalogue of *existing* START solutions.

This never calls the solver.  A family denotes the origin of a solve, not a
verified physical mode.  M2 values retain their native uniform-grid station;
the nearest 841-grid station is only a place to try that value as a guess.
"""
import argparse
import bisect
import collections
import configparser
import datetime
import hashlib
import json
import math
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parent.parent
CASE = ROOT.parent
FIXED = CASE / 'neutral_xomega_20261022'
TRIAL = CASE / 'neutral_beta_envelope_20261022'
OUTPUT = ROOT / 'outputs' / 'seed_catalog.json'
COLUMNS = ['outer', 'inner', 'iloc', 'x', 'Re', 'alpha_r', 'alpha_i',
           'beta_r', 'beta_i', 'omega_r', 'omega_i', 'growth', 'branch', 'status']
FAMILIES = {
    'fixed_M1': {'mode': 'M1', 'native_grid': 'main_841'},
    'fixed_M2': {'mode': 'M2', 'native_grid': 'uniform_661'},
    'fixed_M3': {'mode': 'M3', 'native_grid': 'main_841'},
    'trial_M1': {'mode': 'M1', 'native_grid': 'main_841'},
}


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()


def read_grid(path):
    import h5py
    import numpy as np
    with h5py.File(path, 'r') as stream:
        values = np.array(stream['grid'])
    assert np.isfinite(values).all()
    assert np.count_nonzero(values[..., 1]) == 0
    x = values[0, :, 0, 0].tolist()
    assert all(b > a for a, b in zip(x, x[1:]))
    return dict(source_path=str(path), sha256=sha(path),
                shape=list(values.shape), x_nondim=x,
                x_mm=[35.0 * value for value in x],
                imag_nonzero_count=0, coordinate='global wall Cartesian x/Lref')


def nearest_station(x_values, target):
    right = bisect.bisect_left(x_values, target)
    choices = [i for i in (right - 1, right) if 0 <= i < len(x_values)]
    return min(choices, key=lambda i: (abs(x_values[i] - target), i)) + 1


def table(path):
    zone = 0
    with path.open() as stream:
        for line_number, line in enumerate(stream, 1):
            match = re.search(r'ZONE\s+T="row_(\d+)"', line)
            if match:
                zone = int(match[1])
                continue
            try:
                values = [float(v.replace('D', 'E')) for v in line.split()]
            except ValueError:
                continue
            if len(values) != len(COLUMNS):
                continue
            yield line_number, zone, dict(zip(COLUMNS, values))


def terminal_for_field(path):
    log = path.with_name('neutral_run.log')
    if not log.exists():
        return dict(result='MISSING', source_path=str(log))
    text = log.read_text()
    matches = list(re.finditer(r'evaluations=\s*(\d+)\s+solves=\s*(\d+)\s+rejected=\s*(\d+)\s+result=([^\n]+)', text))
    if not matches:
        # The final aggregate is derived from several real runs.
        return dict(result='derived aggregate; see provenance.json', source_path=str(log))
    match = matches[-1]
    return dict(result=match[4].strip(), evaluations=int(match[1]),
                solves=int(match[2]), rejected=int(match[3]), source_path=str(log))


def complex_pair(text):
    return [float(v.strip().replace('D', 'E').replace('d', 'e')) for v in text.split(',')]


def original_seed(mode, grids):
    path = CASE / mode / 'config.ini'
    cfg = configparser.ConfigParser()
    cfg.read(path)
    station = cfg.getint('domain', 'iloc')
    native = 'uniform_661' if mode == 'M2' else 'main_841'
    x = grids[native]['x_nondim'][station - 1]
    omega, beta, alpha = (complex_pair(cfg['lst'][k]) for k in ('omega', 'beta', 'alpha'))
    nearest = nearest_station(grids['main_841']['x_nondim'], x)
    return dict(id='original_config_' + mode, mode=mode, native_grid=native,
                native_iloc=station, iloc=station, x=x, x_mm=35*x,
                omega_r=omega[0], omega_i=omega[1], beta_r=beta[0], beta_i=beta[1],
                alpha_r=alpha[0], alpha_i=alpha[1], status='initial guess; not a converged solution',
                value_kind='original configuration initial guess', source_path=str(path),
                source_sha256=sha(path), main_nearest_iloc=nearest,
                main_nearest_delta_x_mm=35*(grids['main_841']['x_nondim'][nearest-1]-x),
                cross_grid_mapping_is_initial_guess_only=(mode == 'M2'))


def conjugate_seed(seed):
    result = dict(seed)
    result.update(id='derived_conjugate_of_' + seed['id'], derived_from=seed['id'],
                  value_kind='algebraic common conjugate transform; not separately solved',
                  omega_r=-seed['omega_r'], omega_i=seed['omega_i'],
                  beta_r=-seed['beta_r'], beta_i=seed['beta_i'],
                  alpha_r=-seed['alpha_r'], alpha_i=seed['alpha_i'],
                  status='derived initial guess; requires actual solve')
    return result


def build():
    grids = {
        'main_841': read_grid(CASE/'M1/in/grid_vertical2d_N.h5'),
        'uniform_661': read_grid(CASE/'M2/in/grid_vertical2d_N_uniform.h5'),
    }
    sources = [(FIXED / mode / 'final/out/neutral_field.plt', 'fixed_' + mode)
               for mode in ('M1', 'M2', 'M3')]
    # Only this bounded task directory is traversed, exactly one run level.
    sources += [(path, 'trial_M1') for path in sorted((TRIAL/'runs').glob('*/out/neutral_field.plt'))]
    records, source_meta = [], []
    groups = {family: collections.defaultdict(list) for family in FAMILIES}
    for path, family in sources:
        assert path.is_file(), str(path)
        info = FAMILIES[family]
        native = info['native_grid']
        counts = collections.Counter()
        terminal = terminal_for_field(path)
        for line_number, zone, row in table(path):
            counts['numeric_rows'] += 1
            if row['status'] != 0:
                counts['nonzero_status_excluded'] += 1
                continue
            if not all(math.isfinite(v) for v in row.values()):
                counts['nonfinite_excluded'] += 1
                continue
            native_station = int(row['iloc'])
            assert row['iloc'] == native_station
            assert 1 <= native_station <= len(grids[native]['x_nondim'])
            exact_x = grids[native]['x_nondim'][native_station - 1]
            assert abs(row['x'] - exact_x) < 1e-7, (path, line_number, row['x'], exact_x)
            closest = nearest_station(grids['main_841']['x_nondim'], exact_x)
            record_id = len(records)
            rec = dict(row, id=record_id, family=family, mode=info['mode'],
                       native_grid=native, native_iloc=native_station, iloc=native_station,
                       x_mm=35*row['x'], native_grid_x=exact_x,
                       source_path=str(path), source_line=line_number, source_row=zone,
                       source_terminal=terminal['result'],
                       main_nearest_iloc=closest,
                       main_nearest_delta_x_mm=35*(grids['main_841']['x_nondim'][closest-1]-exact_x),
                       value_kind='actual status-zero field solution; useful as initial guess',
                       cross_grid_mapping_is_initial_guess_only=(native != 'main_841'))
            records.append(rec)
            groups[family][str(native_station)].append(record_id)
            counts['accepted_status_zero_rows'] += 1
        source_meta.append(dict(path=str(path), family=family, sha256=sha(path),
                                bytes=path.stat().st_size, counts=dict(counts), terminal=terminal))

    group_meta = {}
    summary = {}
    for family, station_map in groups.items():
        coords = grids[FAMILIES[family]['native_grid']]['x_nondim']
        group_meta[family] = {
            station: dict(native_iloc=int(station), x=coords[int(station)-1],
                          x_mm=35*coords[int(station)-1], record_ids=ids)
            for station, ids in sorted(station_map.items(), key=lambda pair: int(pair[0]))
        }
        ids = [record_id for values in station_map.values() for record_id in values]
        rows = [records[record_id] for record_id in ids]
        summary[family] = dict(
            record_count=len(rows), native_station_count=len(station_map),
            native_station_range=[min(int(v) for v in station_map), max(int(v) for v in station_map)],
            source_file_count=sum(v['family'] == family for v in source_meta),
            actual_ranges={key: [min(row[key] for row in rows), max(row[key] for row in rows)]
                           for key in ('x', 'x_mm', 'omega_r', 'beta_r', 'alpha_r', 'alpha_i')},
            family_is_provenance_not_verified_physical_mode=True)

    nearest = {}
    for main_i, main_x in enumerate(grids['main_841']['x_nondim'], 1):
        nearest[str(main_i)] = {}
        for family, available in group_meta.items():
            chosen = min(available.values(), key=lambda group: (abs(group['x']-main_x), group['native_iloc']))
            nearest[str(main_i)][family] = dict(
                station_group=str(chosen['native_iloc']),
                source_native_iloc=chosen['native_iloc'], source_x=chosen['x'],
                source_x_mm=chosen['x_mm'], target_main_x=main_x,
                signed_distance_mm=35*(chosen['x']-main_x),
                absolute_distance_mm=35*abs(chosen['x']-main_x),
                candidate_record_count=len(chosen['record_ids']),
                mapping_is_initial_guess_only=True)

    original = [original_seed(mode, grids) for mode in ('M1', 'M2', 'M3')]
    baseline_path = FIXED/'baseline_summary.json'
    baseline = []
    for item in json.loads(baseline_path.read_text()):
        seed = dict(next(s for s in original if s['mode'] == item['mode']))
        state = item['state']
        seed.update(id='baseline_converged_' + item['mode'], omega_r=state[5], omega_i=state[6],
                    alpha_r=state[7], alpha_i=state[8], beta_r=state[9], beta_i=state[10],
                    value_kind='previous actually converged baseline; rounded archived output',
                    status='previous solve exit ' + str(item['exit']),
                    source_path=str(baseline_path), source_sha256=sha(baseline_path),
                    original_solver_directory=item['directory'])
        baseline.append(seed)

    result = dict(
        created_at=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),
        schema_version=1, root=str(ROOT),
        notes=[
            'Read-only extraction: this script does not run a solver or alter original inputs.',
            'Every field record retains its actual tuple, native station, and absolute source path/line.',
            'No frequency, beta, alpha, or mode identity is interpolated or blended.',
            'Family labels describe provenance, not independently verified physical mode identity.',
            'A status-zero sample from an unsuccessful overall run remains labelled with that terminal state.',
            'M2 native values are not claimed as solutions on the main 841-grid.',
            'Nearest per-family records can be far from the target; distances are explicit.',
            'Derived conjugate seeds use (-omega,-beta,-conj(alpha)); no separate solve is claimed.',
            'Only fixed final field tables and trial runs at one explicitly bounded directory depth were read.',
        ],
        grids=grids, families=FAMILIES, family_summary=summary,
        source_files=source_meta, records=records, station_groups=group_meta,
        nearest_by_main_station=nearest,
        original_config_seeds=original, baseline_converged_seeds=baseline,
        derived_conjugate_seeds=[conjugate_seed(seed) for seed in original+baseline])
    OUTPUT.parent.mkdir(exist_ok=True)
    OUTPUT.write_text(json.dumps(result, ensure_ascii=False, separators=(',', ':'))+'\n')
    print(json.dumps(dict(output=str(OUTPUT), bytes=OUTPUT.stat().st_size,
                          records=len(records), family_summary=summary,
                          main_station_count=len(nearest), original_config_seeds=len(original),
                          baseline_seeds=len(baseline), derived_conjugate_seeds=6), ensure_ascii=False, indent=2))
    return result


def load(path=OUTPUT):
    return json.loads(Path(path).read_text())


def get_candidates(catalogue, main_station, family, omega=None, beta=None, limit=8):
    """Return separate real source records, never parameter-interpolated seeds.

    First selects the nearest available native station in this provenance
    family.  If omega/beta are supplied, orders its records by that requested
    parameter distance, then by highest observed growth.  This does not prove
    physical mode identity or validate propagation to the requested station.
    """
    target = catalogue['nearest_by_main_station'][str(main_station)][family]
    group = catalogue['station_groups'][family][target['station_group']]
    records = [catalogue['records'][record_id] for record_id in group['record_ids']]
    def key(row):
        distance = (abs(row['omega_r']-omega) if omega is not None else 0.0)
        distance += (abs(row['beta_r']-beta) if beta is not None else 0.0)
        return distance, -row['growth'], row['id']
    unique, seen = [], set()
    for row in sorted(records, key=key):
        signature = tuple(row[key] for key in ('omega_r', 'omega_i', 'beta_r', 'beta_i', 'alpha_r', 'alpha_i'))
        if signature in seen:
            continue
        seen.add(signature)
        unique.append(dict(row, requested_main_iloc=int(main_station),
                           source_to_target_distance_mm=target['absolute_distance_mm'],
                           source_to_target_signed_distance_mm=target['signed_distance_mm']))
        if len(unique) >= limit:
            break
    return unique


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--query-main', type=int)
    parser.add_argument('--family', choices=list(FAMILIES), default='fixed_M1')
    parser.add_argument('--omega', type=float)
    parser.add_argument('--beta', type=float)
    parser.add_argument('--limit', type=int, default=8)
    args = parser.parse_args()
    if args.query_main is None:
        build()
    else:
        print(json.dumps(get_candidates(load(), args.query_main, args.family,
                                        args.omega, args.beta, args.limit), ensure_ascii=False, indent=2))
