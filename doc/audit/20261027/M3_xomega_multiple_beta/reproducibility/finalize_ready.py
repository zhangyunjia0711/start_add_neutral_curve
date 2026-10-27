"""Finalize independently completed contours, with persistent per-attempt logs."""
import concurrent.futures
import contextlib
import json
from pathlib import Path
from common import ROOT
from finalize import finalize


def one(path):
    tag = Path(path).stem
    log = ROOT / 'outputs' / (tag + '_finalization.log')
    k = 1
    while log.exists():
        log = ROOT / 'outputs' / (tag + f'_finalization_attempt{k}.log')
        k += 1
    with log.open('w') as stream, contextlib.redirect_stdout(stream):
        try:
            d = finalize(path)
            result = dict(guide=path, status=d['status'], final=path.replace('_guide.json', '_final.json'),
                          case=d['case'], beta=d['beta'], points=len(d['points']), tips=len(d['physical_tips']), log=str(log))
        except Exception as error:
            result = dict(guide=path, status='unresolved', reason=repr(error), log=str(log))
            print(result)
    print('FINAL_READY', tag, result['status'], result.get('reason', '')[:160], flush=True)
    return result


def main():
    paths = []
    for case in ['bluntplate', 'rough']:
        for beta in [-30., -20., -10., 0., 10., 20.]:
            tag = f'{case}_beta{beta:+g}'
            previous = ROOT / 'outputs' / (tag + '_components.json')
            temporal = ROOT / 'outputs' / (tag + '_temporal_components.json')
            if temporal.exists():
                paths.extend(r['guide'] for r in json.loads(temporal.read_text()) if r['status'] == 'closed_verified')
            elif previous.exists():
                paths.extend(json.loads(previous.read_text())['guides'])
    # The downstream negative-frequency pocket is independently completed.
    extra = ROOT / 'outputs/rough_beta-20_c2_retry3_guide.json'
    if extra.exists():
        paths.append(str(extra))
    paths = sorted(set(p for p in paths if json.loads(Path(p).read_text())['status'] == 'closed_verified'))
    results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:
        for result in pool.map(one, paths):
            results.append(result)
            (ROOT / 'outputs/ready_finalization.json').write_text(json.dumps(results, indent=2) + '\n')


if __name__ == '__main__':
    main()
