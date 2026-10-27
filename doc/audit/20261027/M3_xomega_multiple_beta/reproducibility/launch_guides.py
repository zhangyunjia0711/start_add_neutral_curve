"""Bounded parallel numerical jobs, each with an independent output log."""
import concurrent.futures
import json
import os
import subprocess
import sys
from common import ROOT


def run(job):
    case, beta = job
    tag = f'{case}_beta{beta:+g}'
    log = ROOT / 'outputs' / (tag + '_guide.log')
    assert not log.exists(), 'Keep an earlier execution log'
    with log.open('w') as stream:
        result = subprocess.run([sys.executable, str(ROOT / 'work/run_guides.py'), case, str(beta)],
                                stdout=stream, stderr=subprocess.STDOUT,
                                env=dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1'))
    print('GUIDE_JOB', tag, result.returncode, flush=True)
    return dict(case=case, beta=beta, exit=result.returncode, log=str(log))


def main():
    jobs = [(case, beta) for beta in [-30., -20., -10., 0., 10., 20.]
            for case in ['bluntplate', 'rough'] if not (case == 'rough' and beta in [-30., 20.])]
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        for result in pool.map(run, jobs):
            results.append(result)
            (ROOT / 'outputs/guide_jobs.json').write_text(json.dumps(results, indent=2) + '\n')


if __name__ == '__main__':
    main()
