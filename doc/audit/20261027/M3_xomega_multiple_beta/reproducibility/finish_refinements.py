"""Resolve display density and the optimized representative's positive-beta bound."""
import concurrent.futures
import copy
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from common import ROOT, CASES, sha
from solver import qep, mac, original, Pencil
from physical import Pair
from physical_trace import PhysicalFamily
from beta_limit import solve


def dense():
    source = ROOT / 'outputs/rough_beta-30_temporal_c0_final.json'
    d = json.loads(source.read_text())
    case, beta, tag = 'rough', -30., 'rough_beta-30_dense'
    native = [p for p in d['points'] if p['kind'] == 0]
    tips = sorted(d['physical_tips'], key=lambda p: p['xi'])
    assert len(tips) == 2
    family = PhysicalFamily(case, beta, tag)
    branches = [[], []]
    for station in sorted({p['xi'] for p in native}):
        roots = sorted([p for p in native if p['xi'] == station], key=lambda p: p['omega_r'])
        assert len(roots) == 2
        for k in range(2):
            branches[k].append(roots[k])
    xs = np.linspace(tips[0]['xi'], tips[1]['xi'], 29)[1:-1]
    checks = []
    for j, x in enumerate(xs):
        if abs(x - round(x)) < 1e-9:
            continue
        i, theta = int(np.floor(x)), float(x - np.floor(x))
        A, G = family.physical(i, theta)
        roots = []
        for k in range(2):
            ref = min(branches[k], key=lambda p: abs(p['xi'] - x))
            alpha, omega, q = ref['alpha_r'], ref['omega_r'], np.load(ref['vector_file'])
            reference_q = q.copy()
            for it in range(30):
                K = A[0] + alpha * A[1] + alpha * alpha * A[2]
                w, v, r, p = qep([K, G, G * 0], omega, q, True)
                den = np.vdot(p, G @ v)
                wa = -np.vdot(p, (A[1] + 2 * alpha * A[2]) @ v) / den
                if abs(w.imag) < 1e-11 and np.linalg.norm((K + w.real * G) @ v) < 1e-8:
                    break
                alpha += float(np.clip(-w.imag / wa.imag, -1., 1.))
                omega, q = w.real, v
            else:
                raise RuntimeError('Dense physical neutral root failed')
            assert mac(reference_q, v) > .95
            qfile = ROOT / 'outputs' / f'{tag}_x{j:02d}_side{k}.npy'
            np.save(qfile, v)
            row = dict(case=case, kind=3, xi=float(x), x_mm=(1 - theta) * CASES[case]['x_mm'][i - 1] + theta * CASES[case]['x_mm'][i],
                       alpha_r=float(alpha), alpha_i=0., beta_r=beta, beta_i=0., omega_r=float(w.real), omega_i=float(w.imag),
                       growth=0., right_residual=float(np.linalg.norm((K + w.real * G) @ v)), raw_MAC=mac(reference_q, v),
                       vector_file=str(qfile), generation=str(family.pairs[i].build(theta) / 'generation.json'))
            row['x_flat_mm'] = row['x_mm'] - CASES[case]['flat_start_mm']
            roots.append(row)
            if j in [3, 13, 23]:
                fresh, qcheck, receipt = original(family.pairs[i].build(theta), 5, row, beta, tag + f'_x{j:02d}_side{k}_check')
                assert abs(fresh['growth']) < 1e-8 and mac(qcheck, v) > 1 - 1e-7
                checks.append(dict(xi=float(x), side=k, receipt=receipt, growth=fresh['growth'], raw_MAC=mac(qcheck, v)))
        roots.sort(key=lambda p: p['omega_r'])
        assert roots[1]['omega_r'] - roots[0]['omega_r'] > 1e-6
        for k in range(2):
            branches[k].append(roots[k])
    points = [tips[0]] + sorted(branches[0], key=lambda p: p['xi']) + [tips[1]] + sorted(branches[1], key=lambda p: -p['xi'])
    adjacent = [mac(np.load(a['vector_file']), np.load(b['vector_file'])) for a, b in zip(points, points[1:] + points[:1])]
    assert min(adjacent) > .95
    points.append(dict(points[0], kind=2, verified_closure=True))
    d.update(points=points, dense_source=str(source), dense_source_sha256=sha(source), added_physical_profile_points=sum(p['kind'] == 3 for p in points),
             dense_original_checks=checks, minimum_adjacent_raw_MAC=min(adjacent),
             display_order='Increasing xi on lower omega, decreasing xi on upper omega; original independently verified closure retained')
    d['pre_density_sampled_omega_range'] = d['omega_range']
    d['omega_range'] = [min(p['omega_r'] for p in points), max(p['omega_r'] for p in points)]
    target = ROOT / 'outputs/rough_beta-30_dense_final.json'
    target.write_text(json.dumps(d, indent=2) + '\n')
    print('DENSE_FINAL', len(points), min(adjacent), flush=True)


def beta_bound():
    source = ROOT / 'outputs/rough_positive_beta_limit_upstream.json'
    d = json.loads(source.read_text())
    peak = max(d['records'], key=lambda p: p['beta_r'])
    pairs, cache, samples = {}, {}, []

    def evaluate(x):
        key = float(x).hex()
        if key in cache:
            return cache[key]
        i = int(np.floor(x))
        theta = float(x - i)
        if i not in pairs:
            pairs[i] = Pair('rough', i, 0., 'rough_beta_upper_refined')
        p = pairs[i].pencil(theta)
        row, q, _ = solve(p, peak, np.load(peak['vector_file']))
        row.update(xi=float(x), x_mm=(1 - theta) * CASES['rough']['x_mm'][i - 1] + theta * CASES['rough']['x_mm'][i],
                   generation=str(pairs[i].build(theta) / 'generation.json'))
        samples.append({k: row[k] for k in ['xi', 'x_mm', 'beta_r', 'omega_r']})
        cache[key] = row, q, p
        return cache[key]

    choices = [evaluate(float(peak['xi']))]
    for low, high in [(peak['xi'] - .25, peak['xi']), (peak['xi'], peak['xi'] + .25)]:
        choices += [evaluate(float(low)), evaluate(float(high))]
        opt = minimize_scalar(lambda x: -evaluate(float(x))[0]['beta_r'], bounds=(low + 1e-6, high - 1e-6),
                              method='bounded', options=dict(xatol=2e-6, maxiter=35))
        choices.append(evaluate(float(opt.x)))
    row, q, base = max(choices, key=lambda r: r[0]['beta_r'])
    fresh, v, receipt = original(Path(row['generation']).parent, 5, row, row['beta_r'], 'rough_beta_upper_refined_independent')
    assert abs(fresh['growth']) < 1e-8 and mac(q, v) > 1 - 1e-7
    qfile = ROOT / 'outputs/rough_beta_upper_refined_q.npy'
    np.save(qfile, q)
    tests = []
    for beta in [row['beta_r'] - .01, row['beta_r'] + .01, 20.]:
        p = copy.copy(base)
        evaluations = []
        p.omega = row['omega_r']
        p.C = [list(c) for c in p.baseC]
        p.C[0][0] = p.C[0][0] + p.omega * p.G
        a0, q0, r0 = p.solve_alpha(row['beta_r'], row['alpha_r'], q, True)
        b0 = row['beta_r']
        while abs(b0 - beta) > 1e-12:
            db = np.sign(beta - b0) * min(.1, abs(beta - b0))
            target = a0 + db * complex(r0['alpha_beta_r'], r0['alpha_beta_i'])
            a0, q0, r0 = p.solve_alpha(b0 + db, target, q0, True)
            b0 += db
        states = [(row['omega_r'], a0, q0, r0)]
        def alpha_i(w):
            ww, a, v, r = min(states, key=lambda s: abs(s[0] - w))
            while abs(ww - w) > 1e-12:
                dw = np.sign(w - ww) * min(.1, abs(w - ww))
                p.omega = float(ww + dw)
                p.C = [list(c) for c in p.baseC]
                p.C[0][0] = p.C[0][0] + p.omega * p.G
                target = a + dw * complex(r['alpha_omega_r'], r['alpha_omega_i'])
                a, v, r = p.solve_alpha(float(beta), target, v, True)
                ww += dw
                states.append((ww, a, v, r))
            evaluations.append(dict(omega_r=float(w), growth=float(-a.imag), alpha_r=float(a.real), MAC=mac(q,v)))
            return float(a.imag)
        opt = minimize_scalar(alpha_i, bounds=(row['omega_r'] - 2, row['omega_r'] + 2), method='bounded', options=dict(xatol=1e-8))
        value = alpha_i(float(opt.x))
        assert value < 0 if beta < row['beta_r'] else value > 0
        test = dict(beta=float(beta), omega_at_local_maximum=float(opt.x), maximum_growth=-value, evaluations=evaluations)
        if beta == 20.:
            state = dict(alpha_r=evaluations[-1]['alpha_r'], alpha_i=value, beta_r=beta, omega_r=float(opt.x))
            actual, vector, check = original(Path(row['generation']).parent, 5, state, beta, 'rough_beta20_local_damped_independent')
            assert actual['growth'] < 0 and abs(actual['growth'] + value) < 1e-7
            test['original_check'] = check
        tests.append(test)
    out = dict(critical=row, vector_file=str(qfile), original_check=receipt, physical_samples=samples,
               beta_side_tests=tests, upstream_source=str(source),
               scope='Local-profile refinement of the largest resolved positive-beta merger. Tracked branch, not all eigenmodes; unresolved far-domain projections retained.')
    (ROOT / 'outputs/rough_positive_beta_limit_refined.json').write_text(json.dumps(out, indent=2) + '\n')
    print('REFINED_BETA_BOUND', row['beta_r'], row['x_mm'], tests[-1]['maximum_growth'], flush=True)


if __name__ == '__main__':
    import sys
    if '--beta-only' in sys.argv:
        beta_bound()
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=2) as pool:
            a, b = pool.submit(dense), pool.submit(beta_bound)
            a.result()
            b.result()
