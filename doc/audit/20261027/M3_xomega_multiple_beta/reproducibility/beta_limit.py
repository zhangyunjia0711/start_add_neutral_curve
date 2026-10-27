"""Supplementary positive-beta cutoff of the optimized tracked representative."""
import json
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares, minimize_scalar
from common import ROOT, CASES, sha
from solver import Pencil, qep, mac, original
from physical import Pair
W = 6.0727867021443425 / 50


def solve(pencil, point, q):
    b0, w0 = point['beta_r'], point['omega_r']
    def evaluate(x):
        b, omega = x
        pencil.omega = float(omega)
        pencil.C = [list(row) for row in pencil.baseC]
        pencil.C[0][0] = pencil.C[0][0] + omega * pencil.G
        a, v, r = pencil.solve_alpha(float(b), point['alpha_r'], q, True)
        return np.array([a.imag, 10 * r['alpha_omega_i']]), a, v, r
    x0 = np.array([b0, w0])
    opt = least_squares(lambda x: evaluate(x)[0], x0, jac='3-point',
                        bounds=(x0 - [5., 10.], x0 + [5., 10.]),
                        xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=100)
    F, a, v, r = evaluate(opt.x)
    assert max(abs(F)) < 1e-8 and r['right_residual'] < 1e-7 and mac(q, v) > .95, (opt.message, opt.x, F, r)
    b, omega = opt.x
    return dict(alpha_r=float(a.real), alpha_i=float(a.imag), beta_r=float(b), omega_r=float(omega),
                f_kHz=float(omega / W), right_residual=r['right_residual'], raw_MAC=mac(q, v)), v, []


def main():
    case = 'rough'
    source = ROOT.parent / 'neutral_m3_frequency_limits_20261027/outputs/rough_f250kHz_neutral.json'
    d = json.loads(source.read_text())
    seed = max([p for p in d['points'] if p['kind'] == 0], key=lambda p: p['beta_r'])
    i0 = int(seed['xi'])
    initial = (seed, np.load(seed['vector_file']))
    records, failures = [], []
    for direction in [1, -1]:
        point, q = initial
        stations = range(i0, CASES[case]['stations'] + 1) if direction > 0 else range(i0 - 1, 1, -1)
        for i in stations:
            try:
                p = Pencil(case, i, 0.)
                point, q, derivative = solve(p, point, q)
                vf = ROOT / 'outputs' / f'rough_positive_beta_limit_i{i:03d}.npy'
                np.save(vf, q)
                records.append(dict(point, station=i, xi=float(i), x_mm=CASES[case]['x_mm'][i - 1],
                                    vector_file=str(vf), derivative_check=derivative))
            except Exception as error:
                failures.append(dict(station=i, direction=direction, reason=repr(error)))
                break
    records.sort(key=lambda r: r['station'])
    (ROOT / 'outputs/rough_positive_beta_limit_native.json').write_text(json.dumps(dict(records=records, failures=failures), indent=2) + '\n')
    print('BETA_NATIVE', len(records), failures, flush=True)
    peak = max(records, key=lambda p: p['beta_r'])
    by = {p['station']: p for p in records}
    pairs, cache, samples = {}, {}, []

    def evaluate(x):
        key = float(x).hex()
        if key in cache:
            return cache[key]
        i = min(int(np.floor(x)), CASES[case]['stations'] - 1)
        t = float(x - i)
        assert i in by and i + 1 in by
        if i not in pairs:
            pairs[i] = Pair(case, i, 0., 'rough_positive_beta_limit')
        a, b = by[i], by[i + 1]
        guess = {k: (1 - t) * a[k] + t * b[k] for k in ['alpha_r', 'beta_r', 'omega_r']}
        row, q, _ = solve(pairs[i].pencil(t), guess, np.load((a if t < .5 else b)['vector_file']))
        row.update(xi=float(x), x_mm=(1 - t) * CASES[case]['x_mm'][i - 1] + t * CASES[case]['x_mm'][i],
                   generation=str(pairs[i].build(t) / 'generation.json'))
        samples.append({k: row[k] for k in ['xi', 'x_mm', 'beta_r', 'omega_r', 'alpha_r']})
        cache[key] = row, q, pairs[i].pencil(t)
        return cache[key]

    choices = []
    for i in range(max(min(by), peak['station'] - 2), min(max(by), peak['station'] + 2)):
        choices += [evaluate(float(i)), evaluate(float(i + 1))]
        opt = minimize_scalar(lambda t: -evaluate(i + float(t))[0]['beta_r'], bounds=(1e-6, 1 - 1e-6),
                              method='bounded', options=dict(xatol=2e-6, maxiter=35))
        choices.append(evaluate(i + float(opt.x)))
    row, q, p = max(choices, key=lambda r: r[0]['beta_r'])
    point, v, receipt = original(Path(row['generation']).parent, 5, row, row['beta_r'], 'rough_positive_beta_limit_independent')
    assert abs(point['growth']) < 1e-8 and mac(q, v) > 1 - 1e-7
    qfile = ROOT / 'outputs/rough_positive_beta_limit_q.npy'
    np.save(qfile, q)
    out = dict(case=case, status='tracked_native_branch_and_local_profile_peak', source=str(source),
               source_sha256=sha(source), records=records, failures=failures, physical_samples=samples,
               critical=row, vector_file=str(qfile), original_check=receipt,
               scope='Supplementary bound on this tracked positive-beta merger branch; no all-spectrum stability claim')
    (ROOT / 'outputs/rough_positive_beta_limit.json').write_text(json.dumps(out, indent=2) + '\n')
    print('BETA_LIMIT', row, 'native coverage', len(records), failures, flush=True)


if __name__ == '__main__':
    main()
