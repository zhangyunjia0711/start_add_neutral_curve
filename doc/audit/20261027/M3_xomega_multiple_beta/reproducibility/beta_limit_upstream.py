"""Bridge coarse upstream profiles before attempting each native beta-fold solve."""
import json
import numpy as np
from common import ROOT, CASES
from physical import Pair
from beta_limit import solve


def main():
    d = json.loads((ROOT / 'outputs/rough_positive_beta_limit.json').read_text())
    row = d['critical']
    q = np.load(d['vector_file'])
    x0 = int(round(row['xi']))
    records, failures, pairs = [], [], {}
    out = ROOT / 'outputs/rough_positive_beta_limit_upstream.json'
    for x in np.arange(x0 - .25, 1.99, -.25):
        try:
            i = min(int(np.floor(x)), CASES['rough']['stations'] - 1)
            t = float(x - i)
            if i not in pairs:
                pairs[i] = Pair('rough', i, 0., 'rough_beta_limit_upstream')
            p = pairs[i].pencil(t)
            row, q, _ = solve(p, row, q)
            vf = ROOT / 'outputs' / f'rough_beta_limit_upstream_x{x:g}.npy'
            np.save(vf, q)
            row.update(xi=float(x), x_mm=(1 - t) * CASES['rough']['x_mm'][i - 1] + t * CASES['rough']['x_mm'][i],
                       vector_file=str(vf), generation=str(pairs[i].build(t) / 'generation.json'))
            records.append(row)
            while len(pairs[i].cache) > 4:
                del pairs[i].cache[next(iter(pairs[i].cache))]
            if len(records) % 10 == 0:
                out.write_text(json.dumps(dict(records=records, failures=failures), indent=2) + '\n')
                print('UPSTREAM_BETA', x, row['beta_r'], row['omega_r'], flush=True)
        except Exception as error:
            failures.append(dict(xi=float(x), reason=repr(error)))
            break
    out.write_text(json.dumps(dict(records=records, failures=failures), indent=2) + '\n')
    print('UPSTREAM_BETA_FINAL', len(records), max(r['beta_r'] for r in records) if records else None, failures, flush=True)


if __name__ == '__main__':
    main()
