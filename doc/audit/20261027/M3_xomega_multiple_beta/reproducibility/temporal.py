"""Exact frequency pencil continuation, matching the verified beta=-7 method."""
import json
import numpy as np
from common import ROOT
from family import Family as Matrices
from solver import qep, mac

S = np.array([50., 20.])


class Family(Matrices):
    def evaluate(self, z, q, target, side=1):
        xi, alpha = z * S
        A, G, dA, dG, i, theta = self.pencil(float(xi), side)
        K = A[0] + alpha * A[1] + alpha * alpha * A[2]
        omega, v, r, p = qep([K, G, G * 0], target, q, True)
        den = np.vdot(p, G @ v)
        assert abs(den) > 1e-14, 'Temporal projection degeneracy'
        wa = -np.vdot(p, (A[1] + 2 * alpha * A[2]) @ v) / den
        wx = -np.vdot(p, (dA[0] + alpha * dA[1] + alpha * alpha * dA[2] + omega * dG) @ v) / den
        r.update(xi=float(xi), x_mm=float((1 - theta) * self.x[i - 1] + theta * self.x[i]),
                 alpha_r=float(alpha), alpha_i=0., beta_r=self.beta, beta_i=0.,
                 omega_r=float(omega.real), omega_i=float(omega.imag), growth=0.,
                 omega_alpha_i=float(wa.imag), omega_xi_i=float(wx.imag),
                 original_real_frequency_polynomial_residual=float(np.linalg.norm((K + omega.real * G) @ v)))
        return r, v, p

    def gradient(self, r):
        return np.array([r['omega_xi_i'], r['omega_alpha_i']]) * S

    def correct(self, pred, plane, q, omega, maxit=20):
        z, tr = pred.copy(), []
        for it in range(maxit):
            r, v, p = self.evaluate(z, q, omega)
            F = np.array([r['omega_i'], np.dot(z - pred, plane)])
            tr.append(dict(point=r, F=F.tolist()))
            if abs(F[0]) < 1e-11 and abs(F[1]) < 1e-10 and r['original_real_frequency_polynomial_residual'] < 1e-8:
                return z, r, v, p, tr
            dz = np.linalg.solve(np.vstack([self.gradient(r), plane]), -F)
            if np.linalg.norm(dz) > .2:
                dz *= .2 / np.linalg.norm(dz)
            z += dz
            q, omega = v, r['omega_r']
        raise RuntimeError('Temporal neutral corrector limit')


def trace(seed, tag, family=None):
    if family is None:
        family = Family(seed['case'], seed['beta'])
    z = np.array([seed['point']['xi'], seed['point']['alpha_r']]) / S
    q = np.load(seed['vector_file'])
    r, q, p = family.evaluate(z, q, seed['point']['omega_r'])
    assert abs(r['omega_i']) < 1e-8
    grad = family.gradient(r)
    t = np.array([grad[1], -grad[0]])
    t /= np.linalg.norm(t)
    if t[0] < 0:
        t = -t
    z0, q0, t0 = z.copy(), q.copy(), t.copy()
    wave = ROOT / 'work' / tag
    wave.mkdir(exist_ok=True)
    out = ROOT / 'outputs' / (tag + '_guide.json')
    assert not out.exists()
    state = dict(case=seed['case'], beta=seed['beta'], seed=seed, scale=S.tolist(),
                 status='running', records=[], failures=[], model=getattr(family, 'model', 'Exact frequency eigenpencil; native-operator interpolation guide'))

    def store(z, r, q, p, t):
        j = len(state['records'])
        qf, pf = wave / f'q{j:04d}.npy', wave / f'p{j:04d}.npy'
        np.save(qf, q)
        np.save(pf, p)
        state['records'].append(dict(index=j, z=z.tolist(), point=r, tangent=t.tolist(),
                                    vector_file=str(qf), left_vector_file=str(pf)))
        if j % 20 == 0:
            out.write_text(json.dumps(state, indent=2) + '\n')

    store(z, r, q, p, t)
    step, arc = .015, 0.
    for loop in range(5000):
        if abs(z[0] * S[0] - round(z[0] * S[0])) < 1e-9:
            r, q, p = family.evaluate(z, q, r['omega_r'], side=np.sign(t[0]))
            grad = family.gradient(r)
            tt = np.array([grad[1], -grad[0]])
            tt /= np.linalg.norm(tt)
            if tt[0] * t[0] < 0:
                tt = -tt
            t = tt
        oldz, oldq, oldsection = z.copy(), q.copy(), float(np.dot(z - z0, t0))
        accepted = False
        for retry in range(20):
            try:
                pred = z + step * t
                xi, px = z[0] * S[0], pred[0] * S[0]
                boundary = int(np.floor(xi + 1e-8)) + 1 if t[0] > 0 else int(np.ceil(xi - 1e-8)) - 1
                assert 1 <= boundary <= family.count, 'Original domain boundary'
                bridge = (px - boundary) * t[0] >= 0
                if bridge:
                    pred = z + (boundary / S[0] - z[0]) / t[0] * t
                    zz, rr, v, pp, tr = family.correct(pred, np.array([1., 0.]), q, r['omega_r'])
                else:
                    zz, rr, v, pp, tr = family.correct(pred, t, q, r['omega_r'])
                grad = family.gradient(rr)
                tt = np.array([grad[1], -grad[0]])
                tt /= np.linalg.norm(tt)
                if bridge:
                    if tt[0] * t[0] < 0:
                        tt = -tt
                else:
                    if np.dot(tt, t) < 0:
                        tt = -tt
                    assert np.dot(tt, t) > .4, 'Temporal sharp turn'
                assert mac(oldq, v) > .95 and np.linalg.norm(zz - oldz) < 2 * step
                assert np.linalg.norm(zz - pred) < .5 * step
                z, r, q, p, t = zz, rr, v, pp, tt
                arc += np.linalg.norm(z - oldz)
                store(z, r, q, p, t)
                accepted = True
                break
            except Exception as error:
                state['failures'].append(dict(index=len(state['records']) - 1, step=step, reason=repr(error)))
                step /= 2
                if step < 1e-6:
                    break
        if not accepted:
            state['status'] = 'unresolved_minimum_step'
            break
        section, distance = float(np.dot(z - z0, t0)), float(np.linalg.norm(z - z0))
        if len(state['records']) > 20 and distance < 1e-7 and mac(q0, q) > 1 - 1e-7 and np.dot(t, t0) > .7:
            state.update(status='closed_verified', closure=dict(passed=True, distance=distance, tangent_dot=float(np.dot(t, t0)), raw_MAC=mac(q0, q)))
            break
        if len(state['records']) > 30 and arc > .3 and oldsection < 0 <= section and distance < 2 * step and np.dot(t, t0) > .7:
            normal = np.array([-t0[1], t0[0]])
            pred = oldz + (-oldsection) / (section - oldsection) * (z - oldz) + .001 * normal
            pred -= np.dot(pred - z0, t0) * t0
            zz, rr, v, pp, tr = family.correct(pred, t0, q, r['omega_r'])
            if np.linalg.norm(zz - z0) < 1e-7 and mac(q0, v) > 1 - 1e-7:
                store(zz, rr, v, pp, t0)
                state.update(status='closed_verified', closure=dict(passed=True, distance=float(np.linalg.norm(zz-z0)),
                             raw_MAC=mac(q0,v), trace=tr, initial_distance=float(np.linalg.norm(pred-z0))))
                break
        if len(state['records']) % 100 == 0:
            print('TEMPORAL', tag, len(state['records']), r['xi'], r['omega_r'], flush=True)
        if len(tr) <= 4 and mac(oldq, q) > .99:
            step = min(.06, step * 1.15)
    else:
        state['status'] = 'unresolved_point_limit'
    state['finished'] = True
    out.write_text(json.dumps(state, indent=2) + '\n')
    print('TEMPORAL_FINAL', tag, state['status'], len(state['records']), flush=True)
    return state, out


def one(job):
    case, beta = job
    from trace import covered
    seeds = json.loads((ROOT / 'outputs/seeds.json').read_text())['results']
    guides = []
    results = []
    for j, seed in enumerate(seeds):
        if seed['case'] != case or seed['beta'] != beta or seed['status'] != 'neutral_corrected':
            continue
        if any(g['status'] == 'closed_verified' and covered(seed, g) for g, _ in guides):
            continue
        tag = f'{case}_beta{beta:+g}_temporal_c{len(guides)}'
        g, path = trace(seed, tag)
        guides.append((g, path))
        results.append(dict(guide=str(path), status=g['status']))
        (ROOT / 'outputs' / f'{case}_beta{beta:+g}_temporal_components.json').write_text(json.dumps(results, indent=2) + '\n')
    return results


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('case')
    parser.add_argument('beta', type=float)
    args = parser.parse_args()
    one((args.case, args.beta))
