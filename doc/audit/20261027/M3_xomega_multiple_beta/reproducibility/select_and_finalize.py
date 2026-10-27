"""Audit every corrected seed against complete oriented contours, then finalize."""
import concurrent.futures
import json
from pathlib import Path
import numpy as np
from common import ROOT
from trace import covered
from finalize import finalize


def main():
    seeds = json.loads((ROOT / 'outputs/seeds.json').read_text())['results']
    selection = []
    for cp in sorted((ROOT / 'outputs').glob('*_components.json')):
        component = json.loads(cp.read_text())
        valid = []
        for source in component['guides']:
            g = json.loads(Path(source).read_text())
            if g['status'] != 'closed_verified':
                continue
            if np.dot(g['records'][0]['tangent'], g['records'][-1]['tangent']) <= .7:
                continue
            if any(covered(g['seed'], h) and covered(h['seed'], g) for _, h in valid):
                continue
            valid.append((source, g))
        matched, unmatched = [], []
        for j, seed in enumerate(seeds):
            if seed['case'] == component['case'] and seed['beta'] == component['beta'] and seed['status'] == 'neutral_corrected':
                matches = [p for p, g in valid if covered(seed, g)]
                (matched if matches else unmatched).append(dict(seed_index=j, guides=matches))
        selection.append(dict(case=component['case'], beta=component['beta'], guides=[p for p, _ in valid],
                              matched_seeds=matched, unmatched_seeds=unmatched))
        print('SELECTION', component['case'], component['beta'], len(valid), 'unmatched', len(unmatched), flush=True)
    (ROOT / 'outputs/guide_selection.json').write_text(json.dumps(selection, indent=2) + '\n')
    if any(r['unmatched_seeds'] for r in selection):
        print('UNRESOLVED_SEEDS_RETAINED; finalize independently completed contours first', flush=True)
    paths = [p for r in selection for p in r['guides']]
    results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:
        futures = {pool.submit(finalize, p): p for p in paths}
        for future in concurrent.futures.as_completed(futures):
            path = futures[future]
            try:
                d = future.result()
                record = dict(guide=path, status=d['status'], final=path.replace('_guide.json', '_final.json'),
                              points=len(d['points']), tips=len(d['physical_tips']))
            except Exception as error:
                record = dict(guide=path, status='unresolved', reason=repr(error))
                print('FINAL_FAIL', record, flush=True)
            results.append(record)
            (ROOT / 'outputs/finalization_jobs.json').write_text(json.dumps(results, indent=2) + '\n')
    assert all(r['status'] == 'closed_verified' for r in results), 'Physical finalization failures retained'


if __name__ == '__main__':
    main()
