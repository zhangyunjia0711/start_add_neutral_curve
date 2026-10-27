"""Re-solve every original-station crossing; do not assume predictor nodes hit it."""
import math
import numpy as np
from common import ROOT
from family import Family, SCALE
from solver import mac


def native_records(g, tag):
    family = Family(g['case'], g['beta'])
    # Closure was independently corrected. Canonicalize its final occurrence
    # to the exact start before enumerating integer station intersections.
    records = g['records'][:-1] + [g['records'][0]] if g['status'] == 'closed_verified' else g['records']
    result = []
    directory = ROOT / 'work' / (tag + '_native_crossings')
    directory.mkdir(exist_ok=True)
    for j, (left, right) in enumerate(zip(records, records[1:])):
        a, b = left['point'], right['point']
        lo, hi = min(a['xi'], b['xi']), max(a['xi'], b['xi'])
        if hi - lo < 1e-10:
            continue
        for station in range(math.ceil(lo - 1e-8), math.floor(hi + 1e-8) + 1):
            theta = (station - a['xi']) / (b['xi'] - a['xi'])
            if theta < -1e-6 or theta > 1 + 1e-6:
                continue
            theta = float(np.clip(theta, 0., 1.))
            alpha = (1 - theta) * a['alpha_r'] + theta * b['alpha_r']
            omega = (1 - theta) * a['omega_r'] + theta * b['omega_r']
            reference = left if theta < .5 else right
            q, p = np.load(reference['vector_file']), np.load(reference['left_vector_file'])
            pred = np.array([station, alpha, omega]) / SCALE
            z, r, v, w, _, tr = family.correct(pred, np.array([1., 0., 0.]), q, p)
            assert abs(r['alpha_r'] - alpha) < max(.3, 1.5 * abs(a['alpha_r'] - b['alpha_r']))
            assert abs(r['omega_r'] - omega) < max(.1, 1.5 * abs(a['omega_r'] - b['omega_r']))
            assert mac(q, v) > .95
            order = j + theta
            same = next((old for old in result if abs(old['point']['xi'] - station) < 1e-8
                         and abs(old['point']['alpha_r'] - r['alpha_r']) < 1e-6
                         and abs(old['point']['omega_r'] - r['omega_r']) < 1e-6), None)
            if same is not None:
                assert mac(np.load(same['vector_file']), v) > 1 - 1e-7
                gap = order - same['order']
                assert gap < 3 or (same['order'] < 2 and order > len(records) - 3), ('Interior retracing of a native root', tag, station, same['order'], order)
                continue
            k = len(result)
            qf, pf = directory / f'q{k:04d}.npy', directory / f'p{k:04d}.npy'
            np.save(qf, v)
            np.save(pf, w)
            result.append(dict(index=j, order=order, point=r, vector_file=str(qf), left_vector_file=str(pf),
                               crossing_sources=[left['index'], right['index']]))
    counts = {i: sum(abs(r['point']['xi'] - i) < 1e-8 for r in result)
              for i in {round(r['point']['xi']) for r in result}}
    odd = [i for i, count in counts.items() if count % 2]
    if odd:
        from temporal import Family as FrequencyFamily, S
        frequency_family = FrequencyFamily(g['case'], g['beta'])
        diagnostics = []
        for station in odd:
            refs = [r for r in records if .1 < abs(r['point']['xi'] - station) < 2.1]
            for ref in refs:
                old = ref['point']
                alpha, omega = old['alpha_r'], old['omega_r']
                q = np.load(ref['vector_file'])
                try:
                    for it in range(30):
                        r, v, p = frequency_family.evaluate(np.array([station, alpha]) / S, q, omega)
                        if abs(r['omega_i']) < 1e-11 and r['original_real_frequency_polynomial_residual'] < 1e-8:
                            break
                        derivative = r['omega_alpha_i']
                        assert abs(derivative) > 1e-12
                        alpha += float(np.clip(-r['omega_i'] / derivative, -.5, .5))
                        q, omega = v, r['omega_r']
                    else:
                        raise RuntimeError('Odd-station frequency root iteration limit')
                    assert abs(r['alpha_r'] - old['alpha_r']) < 5 and mac(np.load(ref['vector_file']), v) > .9
                    if any(abs(x['point']['xi'] - station) < 1e-8 and abs(x['point']['alpha_r'] - r['alpha_r']) < 1e-6 and abs(x['point']['omega_r'] - r['omega_r']) < 1e-6 for x in result):
                        continue
                    z = np.array([station / 50., r['alpha_r'] / 20., r['omega_r'] / 10.])
                    choices = []
                    for j, (a, b) in enumerate(zip(records, records[1:])):
                        za = np.array([a['point']['xi'] / 50., a['point']['alpha_r'] / 20., a['point']['omega_r'] / 10.])
                        zb = np.array([b['point']['xi'] / 50., b['point']['alpha_r'] / 20., b['point']['omega_r'] / 10.])
                        delta = zb - za
                        if np.dot(delta, delta) < 1e-24:
                            continue
                        t = float(np.clip(np.dot(z - za, delta) / np.dot(delta, delta), 0., 1.))
                        choices.append((np.linalg.norm(z - za - t * delta), j + t))
                    distance, order = min(choices)
                    assert distance < .05, ('Repaired root far from guide', distance)
                    k = len(result)
                    qf, pf = directory / f'q{k:04d}.npy', directory / f'p{k:04d}.npy'
                    np.save(qf, v)
                    np.save(pf, p)
                    r['repaired_native_frequency_root'] = True
                    result.append(dict(index=ref['index'], order=order, point=r, vector_file=str(qf),
                                       left_vector_file=str(pf), repair_reference=ref['index'], guide_distance=float(distance)))
                    diagnostics.append(dict(station=station, status='additional_root', alpha=r['alpha_r'], omega=r['omega_r']))
                except Exception as error:
                    diagnostics.append(dict(station=station, status='unresolved_attempt', reason=repr(error)))
            assert sum(abs(r['point']['xi'] - station) < 1e-8 for r in result) % 2 == 0, ('Odd root count persists', tag, station, diagnostics)
        import json
        (ROOT / 'outputs' / (tag + '_native_repair.json')).write_text(json.dumps(diagnostics, indent=2) + '\n')
    return sorted(result, key=lambda r: r['order'])
