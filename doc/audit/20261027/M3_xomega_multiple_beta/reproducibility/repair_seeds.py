"""Additional original-station seeds, retaining all unsuccessful diagnostics."""
import json
import numpy as np
from common import ROOT
from family import Family, SCALE


def main():
    path = ROOT / 'outputs/seeds.json'
    state = json.loads(path.read_text())
    original = ROOT / 'outputs/seeds_initial.json'
    assert not original.exists()
    original.write_text(path.read_text())
    source = ROOT.parent / 'neutral_m3_frequency_limits_20261027/outputs/rough_f-200.0kHz_S_final.json'
    data = json.loads(source.read_text())
    family = Family('rough', -30.)
    results = []
    for source_point in sorted([p for p in data['points'] if p['kind'] == 0], key=lambda p: p['beta_r']):
        j = len(state['results'])
        record = dict(case='rough', beta=-30., source_path=str(source), source=source_point, status='pending')
        try:
            z, r, q, p, t, tr = family.seed(source_point, np.load(source_point['vector_file']))
            assert r['raw_MAC'] > .9
            qfile, pfile = ROOT / 'outputs' / f'seed{j:03d}_q.npy', ROOT / 'outputs' / f'seed{j:03d}_p.npy'
            np.save(qfile, q)
            np.save(pfile, p)
            record.update(status='neutral_corrected', point=r, z=z.tolist(), tangent=t.tolist(),
                          vector_file=str(qfile), left_vector_file=str(pfile), trace=tr)
            print('REPAIRED_SEED', j, r, flush=True)
        except Exception as error:
            record.update(status='unresolved', reason=repr(error))
        state['results'].append(record)
        results.append(record)
        if record['status'] == 'neutral_corrected':
            break
    path.write_text(json.dumps(state, indent=2) + '\n')
    (ROOT / 'outputs/seed_repair.json').write_text(json.dumps(results, indent=2) + '\n')


if __name__ == '__main__':
    main()
