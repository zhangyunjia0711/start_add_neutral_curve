"""Recover multiple fixed-beta seeds from verified neutral-surface points."""
import json
import numpy as np
from common import ROOT, sha
from family import Family, SCALE

BASE = ROOT.parent
BETAS = [-30., -20., -10., -7., 0., 10., 20.]


def sources():
    paths = []
    for case in ['bluntplate', 'rough']:
        for f in [-112, -110, 125, 150, 175, 200]:
            paths.append(BASE / 'neutral_xbeta_extended_frequencies_20261027/outputs' / f'{case}_f{f}kHz_S_final.json')
        for f in [-75, -25, 25, 75]:
            paths.append(BASE / 'neutral_xbeta_extra_frequencies_20261027/outputs' / f'{case}_f{f}kHz_S_final.json')
        paths.append(BASE / 'neutral_xbeta_omega0_20261026/outputs' / f'{case}_final.json')
    manifest = json.loads((BASE / 'neutral_m3_frequency_limits_20261027/outputs/final_numerical_audit.json').read_text())
    paths += [__import__('pathlib').Path(r['path']) for r in manifest['curves']]
    result = []
    for path in paths:
        d = json.loads(path.read_text())
        assert d['status'] in ['closed_verified', 'open_at_domain_boundaries'], (path, d['status'])
        result.append((path, d))
    return result


def main():
    data = sources()
    out = ROOT / 'outputs/seeds.json'
    assert not out.exists(), 'Keep completed seed attempts'
    state = dict(betas=BETAS, sources=[dict(path=str(p), sha256=sha(p)) for p, _ in data], results=[])
    for case in ['bluntplate', 'rough']:
        for beta in BETAS:
            if beta == -7:
                continue
            family = Family(case, beta)
            candidates = []
            for path, d in data:
                if d['case'] != case:
                    continue
                points = [p for p in d['points'] if p['kind'] == 0 and abs(p['xi'] - round(p['xi'])) < 1e-9]
                if not points:
                    continue
                ps = sorted(points, key=lambda p: abs(p['beta_r'] - beta))
                if min(p['beta_r'] for p in points) <= beta <= max(p['beta_r'] for p in points):
                    candidates.append((abs(ps[0]['beta_r'] - beta), path, ps[0]))
                    # Preserve opposite x crossings when they are well separated.
                    second = next((p for p in ps[1:] if abs(p['xi'] - ps[0]['xi']) > 5 and abs(p['beta_r'] - beta) < .25), None)
                    if second:
                        candidates.append((abs(second['beta_r'] - beta), path, second))
            if not candidates:
                for path, d in data:
                    if d['case'] == case:
                        p = min((p for p in d['points'] if p['kind'] == 0), key=lambda p: abs(p['beta_r'] - beta))
                        candidates.append((abs(p['beta_r'] - beta), path, p))
                candidates = sorted(candidates, key=lambda v: v[0])[:3]
            seen = []
            for _, path, source in sorted(candidates, key=lambda v: v[0]):
                if any(abs(source['xi'] - r['xi']) < 1e-9 and abs(source['alpha_r'] - r['alpha_r']) < 1e-4 for r in seen):
                    continue
                seen.append(source)
                j = len(state['results'])
                record = dict(case=case, beta=beta, source_path=str(path), source=source, status='pending')
                try:
                    q = np.load(source['vector_file'])
                    z, point, v, p, t, trace = family.seed(source, q)
                    assert abs(point['xi'] - source['xi']) < 1e-7
                    qfile, pfile = ROOT / 'outputs' / f'seed{j:03d}_q.npy', ROOT / 'outputs' / f'seed{j:03d}_p.npy'
                    np.save(qfile, v)
                    np.save(pfile, p)
                    e, J, _, _, _, _ = family.evaluate(z, v, p)
                    numeric = []
                    for k in range(3):
                        h = 1e-5
                        delta = np.eye(3)[k] * h
                        if k == 0:
                            ep, *_ = family.evaluate(z + delta, v, p)
                            numeric.append((ep - e) / h)
                        else:
                            ep, *_ = family.evaluate(z + delta, v, p)
                            em, *_ = family.evaluate(z - delta, v, p)
                            numeric.append((ep - em) / (2 * h))
                    numeric = np.array([np.real(numeric), np.imag(numeric)])
                    error = np.linalg.norm(J - numeric) / max(np.linalg.norm(J), 1e-30)
                    assert error < .01, ('Jacobian check', error)
                    record.update(status='neutral_corrected', point=point, z=z.tolist(), tangent=t.tolist(),
                                  vector_file=str(qfile), left_vector_file=str(pfile), trace=trace, jacobian_relative_error=float(error))
                    print('SEED', j, case, beta, point['xi'], point['omega_r'], flush=True)
                except Exception as error:
                    record.update(status='unresolved', reason=repr(error))
                    print('SEED_FAIL', j, case, beta, repr(error), flush=True)
                state['results'].append(record)
                out.write_text(json.dumps(state, indent=2) + '\n')
    print('SEEDS_COMPLETE', len(state['results']), sum(r['status'] == 'neutral_corrected' for r in state['results']), flush=True)


if __name__ == '__main__':
    main()
