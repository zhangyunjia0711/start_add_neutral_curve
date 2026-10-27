"""Fixed-beta 3D neutral continuation, crossing every original station exactly."""
import json
import numpy as np
from common import ROOT
from family import Family, SCALE
from solver import mac


def trace(seed, tag, initial_step=.025, max_step=.08):
    family = Family(seed['case'], seed['beta'])
    z, q, p = np.array(seed['z']), np.load(seed['vector_file']), np.load(seed['left_vector_file'])
    e, J, t, q, p, r = family.evaluate(z, q, p)
    assert abs(e) < 1e-10 and r['right_residual'] < 1e-8
    if t[0] < 0:
        t = -t
    z0, q0, p0, t0 = z.copy(), q.copy(), p.copy(), t.copy()
    out, wave = ROOT / 'outputs' / (tag + '_guide.json'), ROOT / 'work' / (tag + '_guide')
    wave.mkdir(exist_ok=True)
    assert not out.exists(), 'Preserve existing trace'
    state = dict(case=seed['case'], beta=seed['beta'], seed=seed, scale=SCALE.tolist(),
                 status='running', records=[], failures=[], model='Guide: native matrices and piecewise native-operator interpolation')

    def write():
        tmp = out.with_suffix('.tmp')
        tmp.write_text(json.dumps(state, indent=2) + '\n')
        tmp.replace(out)

    def store(z, r, q, p, t, step, it):
        j = len(state['records'])
        qf, pf = wave / f'q{j:04d}.npy', wave / f'p{j:04d}.npy'
        np.save(qf, q)
        np.save(pf, p)
        state['records'].append(dict(index=j, z=z.tolist(), point=r, tangent=t.tolist(),
            vector_file=str(qf), left_vector_file=str(pf), step=step, iterations=it))
        if j % 20 == 0:
            write()

    step, arc, half = initial_step, 0., 0
    store(z, r, q, p, t, step, 0)
    for loop in range(8000):
        if abs(z[0] * SCALE[0] - round(z[0] * SCALE[0])) < 1e-9:
            e, J, tt, q, p, r = family.evaluate(z, q, p, side=np.sign(t[0]))
            if tt[0] * t[0] < 0:
                tt = -tt
            t = tt
        oldz, oldq = z.copy(), q.copy()
        oldsection = float(np.dot(z - z0, t0))
        accepted = False
        reached_edge = False
        for retry in range(18):
            try:
                pred = z + step * t
                xi, px = z[0] * SCALE[0], pred[0] * SCALE[0]
                boundary = int(np.floor(xi + 1e-8)) + 1 if t[0] > 0 else int(np.ceil(xi - 1e-8)) - 1
                if not 1 <= boundary <= family.count:
                    reached_edge = True
                    break
                if (px - boundary) * t[0] >= 0:
                    pred = z + (boundary / SCALE[0] - z[0]) / t[0] * t
                    zz, rr, v, w, tt, tr = family.correct(pred, np.array([1., 0., 0.]), q, p, 30)
                    if tt[0] * t[0] < 0:
                        tt = -tt
                else:
                    zz, rr, v, w, tt, tr = family.correct(pred, t, q, p)
                    if np.dot(tt, t) < 0:
                        tt = -tt
                    assert np.dot(tt, t) > .35, 'Sharp neutral turn'
                overlap = mac(oldq, v)
                assert overlap > .95 and np.linalg.norm(zz - oldz) < 2 * step, ('Continuity gate', overlap)
                assert np.linalg.norm(zz - pred) < .5 * step, 'Predictor correction too large'
                z, r, q, p, t = zz, rr, v, w, tt
                arc += np.linalg.norm(z - oldz)
                store(z, r, q, p, t, step, len(tr))
                accepted = True
                break
            except Exception as error:
                state['failures'].append(dict(index=len(state['records']) - 1, step=step, reason=repr(error)))
                step /= 2
                if step < 1e-6:
                    break
        if reached_edge:
            state.setdefault('domain_ends', []).append(dict(point=r, record_index=len(state['records']) - 1))
            if half:
                state['status'] = 'open_at_domain_boundaries'
                break
            half = 1
            state['reverse_half_start_index'] = len(state['records'])
            z, q, p, t = z0.copy(), q0.copy(), p0.copy(), -t0.copy()
            e, J, _, q, p, r = family.evaluate(z, q, p)
            step = initial_step
            store(z, r, q, p, t, step, 0)
            continue
        if not accepted:
            state['status'] = 'unresolved_minimum_step'
            break
        section = float(np.dot(z - z0, t0))
        distance = float(np.linalg.norm(z - z0))
        if not half and len(state['records']) > 30 and distance < 1e-7 and mac(q0, q) > 1 - 1e-7 and np.dot(t, t0) > .7:
            state.update(status='closed_verified', closure=dict(passed=True, scaled_distance=distance,
                         raw_MAC=mac(q0, q), tangent_dot=float(np.dot(t, t0)), method='Actual corrected native-section return'))
            break
        if not half and len(state['records']) > 30 and arc > .5 and oldsection < 0 <= section and distance < max(2 * step, .03) and np.dot(t, t0) > .7:
            normal = np.array([1., 0., 0.])
            normal -= np.dot(normal, t0) * t0
            normal /= np.linalg.norm(normal)
            pred = oldz + (-oldsection) / (section - oldsection) * (z - oldz) + 1e-3 * normal
            pred -= np.dot(pred - z0, t0) * t0
            try:
                zz, rr, v, w, tt, tr = family.correct(pred, t0, q, p)
                closure = dict(scaled_distance=float(np.linalg.norm(zz - z0)), raw_MAC=mac(q0, v),
                               initial_distance=float(np.linalg.norm(pred - z0)), trace=tr)
                closure['passed'] = closure['scaled_distance'] < 1e-7 and closure['raw_MAC'] > 1 - 1e-7
                state['closure'] = closure
                if closure['passed']:
                    store(zz, rr, v, w, t0, 0, len(tr))
                    state['status'] = 'closed_verified'
                    break
            except Exception as error:
                state['failures'].append(dict(kind='closure', reason=repr(error)))
        if len(state['records']) % 100 == 0:
            print('TRACE', tag, len(state['records']), r['xi'], r['omega_r'], flush=True)
        if abs(r['xi'] - round(r['xi'])) < 1e-8:
            for previous in state['records'][1:-30]:
                v = previous['point']
                if abs(v['xi'] - r['xi']) < 1e-8 and abs(v['alpha_r'] - r['alpha_r']) < 1e-6 and abs(v['omega_r'] - r['omega_r']) < 1e-6 and np.dot(previous['tangent'], t) > .99:
                    state.update(status='unresolved_nonseed_cycle', cycle_previous_index=previous['index'])
                    break
            if state['status'] == 'unresolved_nonseed_cycle':
                break
        if len(tr) <= 4 and overlap > .99:
            step = min(max_step, step * 1.15)
    else:
        state['status'] = 'unresolved_point_limit'
    state['finished'] = True
    write()
    print('TRACE_FINAL', tag, state['status'], len(state['records']), flush=True)
    return state, out


def covered(seed, guide):
    s = seed['point']
    matches = [r for r in guide['records'] if abs(r['point']['xi'] - s['xi']) < 1e-7
               and abs(r['point']['alpha_r'] - s['alpha_r']) < 1e-6
               and abs(r['point']['omega_r'] - s['omega_r']) < 1e-6]
    if not matches:
        return False
    return any(mac(np.load(seed['vector_file']), np.load(r['vector_file'])) > 1 - 1e-7 for r in matches)
