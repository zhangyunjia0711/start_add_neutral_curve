"""Native neutral points and independently reassembled physical-profile x folds."""
import copy
import json
import numpy as np
from common import ROOT, CASES, sha
from solver import Pencil, original, mac
from physical import Pair


def physical_tip(case, beta, tag, middle, before, after, order, ordinal):
    guess = middle['point']
    xi, omega, alpha = guess['xi'], guess['omega_r'], guess['alpha_r']
    q = np.load(middle['vector_file'])
    pairs, trace = {}, []

    def evaluate(x, omega, alpha, q):
        assert 1 <= x <= CASES[case]['stations']
        i = min(int(np.floor(x)), CASES[case]['stations'] - 1)
        theta = float(x - i)
        if i not in pairs:
            pairs[i] = Pair(case, i, 0., tag + f'_tip{ordinal}')
        p = copy.copy(pairs[i].pencil(theta))
        p.omega = float(omega)
        p.C = [list(row) for row in p.baseC]
        p.C[0][0] = p.C[0][0] + omega * p.G
        a, v, r = p.solve_alpha(beta, alpha, q, True)
        return np.array([a.imag, r['alpha_omega_i']]), a, v, r, i, theta, p

    F, a, q, r, i, theta, p = evaluate(xi, omega, alpha, q)
    for it in range(30):
        trace.append(dict(xi=float(xi), omega_r=float(omega), alpha_r=float(a.real), F=F.tolist(), point=r))
        if max(abs(F)) < 2e-10:
            break
        hx, hw = 1e-4, 1e-3
        Fx, *_ = evaluate(xi + hx, omega, a, q)
        Fw, *_ = evaluate(xi, omega + hw, a, q)
        J = np.column_stack([(Fx - F) / hx, (Fw - F) / hw])
        delta = np.linalg.solve(J, -F)
        delta /= max(abs(delta[0]) / .5, abs(delta[1]) / 1., 1.)
        oldnorm = np.linalg.norm(F * [1., 5.])
        accepted = False
        damping = 1.
        for attempt in range(14):
            xx, ww = xi + damping * delta[0], omega + damping * delta[1]
            try:
                FF, aa, v, rr, ii, tt, pp = evaluate(xx, ww, a, q)
            except (AssertionError, RuntimeError):
                damping /= 2
                continue
            if np.linalg.norm(FF * [1., 5.]) < oldnorm:
                xi, omega, F, a, q, r, i, theta, p = xx, ww, FF, aa, v, rr, ii, tt, pp
                accepted = True
                break
            damping /= 2
        assert accepted, ('Physical x fold stalled', tag, ordinal, trace[-1])
    assert max(abs(F)) < 1e-8 and abs(xi - guess['xi']) < 2., (tag, ordinal, F, xi, guess['xi'])
    pair = pairs[i]
    checks = []
    from scipy.sparse.linalg import norm
    for t, station in [(0., i), (1., i + 1)]:
        actual, native = pair.pencil(t), Pencil(case, station, 0.)
        errors = [float(norm(x - y)) / max(1., float(norm(y))) for x, y in zip(actual.spatial(beta), native.spatial(beta))]
        assert max(errors) < 1e-9
        checks.append(dict(station=station, errors=errors))
    state = dict(alpha_r=float(a.real), alpha_i=float(a.imag), beta_r=beta, omega_r=float(omega))
    fresh, v, receipt = original(pair.build(theta), 5, state, beta, tag + f'_tip{ordinal}_independent')
    assert abs(fresh['growth']) < 1e-8 and mac(q, v) > 1 - 1e-7
    qfile = ROOT / 'outputs' / (tag + f'_tip{ordinal}_q.npy')
    np.save(qfile, v)
    generation = json.loads((pair.build(theta) / 'generation.json').read_text())
    return dict(case=case, kind=1, xi=float(xi), x_mm=generation['x_mm'],
                x_flat_mm=generation['x_mm'] - CASES[case]['flat_start_mm'],
                alpha_r=fresh['alpha_r'], alpha_i=fresh['alpha_i'], beta_r=beta, beta_i=0.,
                omega_r=float(omega), omega_i=0., growth=fresh['growth'], right_residual=r['right_residual'],
                left_residual=r['left_residual'], raw_MAC=mac(q, v), guide_order=float(order),
                vector_file=str(qfile), original_receipt=receipt, generation=str(pair.build(theta) / 'generation.json'),
                endpoint_checks=checks, turning_condition='Im(alpha)=Im(alpha_omega)=0 at fixed real beta',
                trace=trace, guide_guess=guess)


def finalize(path):
    from pathlib import Path
    path = Path(path)
    g = json.loads(path.read_text())
    assert g['status'] == 'closed_verified' and np.dot(g['records'][0]['tangent'], g['records'][-1]['tangent']) > .7
    case, beta, tag = g['case'], g['beta'], path.stem.replace('_guide', '')
    out = ROOT / 'outputs' / (tag + '_final.json')
    if out.exists():
        d = json.loads(out.read_text())
        assert d['status'] == 'closed_verified'
        return d
    records = g['records'][:-1]
    native, tips = [], []
    from cross_native import native_records
    crossings = native_records(g, tag)
    for j, record in enumerate(crossings):
        r = record['point']
        if abs(r['xi'] - round(r['xi'])) > 1e-8:
            continue
        row = dict(r, case=case, kind=0, xi=float(round(r['xi'])), guide_order=float(record['order']),
                   vector_file=record['vector_file'], left_vector_file=record['left_vector_file'])
        row['x_mm'] = CASES[case]['x_mm'][int(row['xi']) - 1]
        row['x_flat_mm'] = row['x_mm'] - CASES[case]['flat_start_mm']
        if native and abs(native[-1]['xi'] - row['xi']) < 1e-8 and abs(native[-1]['alpha_r'] - row['alpha_r']) < 1e-7:
            continue
        p = Pencil(case, int(row['xi']), row['omega_r'])
        a, q, rr = p.solve_alpha(beta, row['alpha_r'], np.load(row['vector_file']), True)
        assert abs(a.imag) < 1e-8 and mac(q, np.load(row['vector_file'])) > 1 - 1e-7
        row.update(alpha_r=float(a.real), alpha_i=float(a.imag), growth=float(-a.imag),
                   right_residual=rr['right_residual'], left_residual=rr['left_residual'], alpha_omega_i=rr['alpha_omega_i'])
        native.append(row)
    for j, m in enumerate(records):
        before, after = records[(j - 1) % len(records)], records[(j + 1) % len(records)]
        if (m['point']['xi'] - before['point']['xi']) * (after['point']['xi'] - m['point']['xi']) < -1e-12:
            tips.append(physical_tip(case, beta, tag, m, before, after, j, len(tips)))
    assert len(tips) >= 2 and len(tips) % 2 == 0, ('X extrema missing', tag, len(tips))
    points = sorted(native + tips, key=lambda p: p['guide_order'])
    assert points
    counts = {i: sum(abs(p['xi'] - i) < 1e-8 for p in native) for i in {round(p['xi']) for p in native}}
    assert set(counts) == set(range(min(counts), max(counts) + 1)) and all(n >= 2 and n % 2 == 0 for n in counts.values()), counts
    checks = []
    for j in sorted(set(np.linspace(0, len(native) - 1, 10, dtype=int).tolist())):
        row = native[j]
        fresh, q, receipt = original(case, int(row['xi']), row, beta, tag + f'_nativecheck{j:04d}')
        delta = abs(complex(fresh['alpha_r'], fresh['alpha_i']) - complex(row['alpha_r'], row['alpha_i']))
        overlap = mac(q, np.load(row['vector_file']))
        assert delta < 1e-7 and overlap > 1 - 1e-7 and abs(fresh['growth']) < 1e-8
        sides = []
        for offset in [-1e-4, 1e-4]:
            p = Pencil(case, int(row['xi']), row['omega_r'] + offset)
            a, _, rr = p.solve_alpha(beta, row['alpha_r'], q, True)
            sides.append(dict(offset=offset, growth=float(-a.imag), residual=rr['right_residual']))
        assert sides[0]['growth'] * sides[1]['growth'] < 0, (tag, row, sides)
        checks.append(dict(native_index=j, alpha_difference=delta, raw_MAC=overlap, receipt=receipt, frequency_side_checks=sides))
    adjacent = [mac(np.load(a['vector_file']), np.load(b['vector_file'])) for a, b in zip(points, points[1:] + points[:1])]
    assert min(adjacent) > .9, ('Adjacent vector mismatch', tag, min(adjacent))
    points.append(dict(points[0], kind=2, verified_closure=True))
    d = dict(case=case, beta=beta, status='closed_verified', guide=str(path), guide_sha256=sha(path),
             points=points, native_points=len(native), native_stations=len(counts), native_counts=counts,
             physical_tips=tips, original_checks=checks, minimum_adjacent_raw_MAC=min(adjacent),
             x_range_mm=[min(p['x_mm'] for p in points), max(p['x_mm'] for p in points)],
             omega_range=[min(p['omega_r'] for p in points), max(p['omega_r'] for p in points)])
    out.write_text(json.dumps(d, indent=2) + '\n')
    print('FINAL', tag, len(points), len(tips), d['x_range_mm'], d['omega_range'], flush=True)
    return d
