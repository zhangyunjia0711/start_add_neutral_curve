"""Resolve every distinct component represented by the corrected seed catalog."""
import argparse
import json
from common import ROOT
from trace import trace, covered


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('case')
    parser.add_argument('beta', type=float)
    parser.add_argument('--retry', type=int, default=0)
    args = parser.parse_args()
    seeds = json.loads((ROOT / 'outputs/seeds.json').read_text())['results']
    seeds = [(j, r) for j, r in enumerate(seeds) if r['case'] == args.case and r['beta'] == args.beta and r['status'] == 'neutral_corrected']
    tag = f'{args.case}_beta{args.beta:+g}'
    report = ROOT / 'outputs' / (tag + '_components.json')
    if args.retry and report.exists():
        import shutil
        backup = report.with_name(report.stem + f'_attempt{args.retry}.json')
        assert not backup.exists()
        shutil.copy2(report, backup)
    guides = []
    results = []
    for j, seed in seeds:
        match = next((i for i, (g, _) in enumerate(guides) if g['status'] in ['closed_verified', 'open_at_domain_boundaries'] and covered(seed, g)), None)
        if match is not None:
            results.append(dict(seed_index=j, status='already_covered', component=match))
            continue
        k = len(guides)
        suffix = f'_c{k}' + (f'_retry{args.retry}' if args.retry else '')
        path = ROOT / 'outputs' / (tag + suffix + '_guide.json')
        if path.exists():
            g = json.loads(path.read_text())
            assert g['seed'] == seed and g['finished']
        else:
            g, path = trace(seed, tag + suffix)
        guides.append((g, path))
        results.append(dict(seed_index=j, status=g['status'], component=k, guide=str(path)))
        report.write_text(json.dumps(dict(case=args.case, beta=args.beta, seeds=results,
            guides=[str(p) for _, p in guides]), indent=2) + '\n')
    report.write_text(json.dumps(dict(case=args.case, beta=args.beta, seeds=results,
        guides=[str(p) for _, p in guides]), indent=2) + '\n')
    print('COMPONENTS', tag, [(g['status'], len(g['records'])) for g, _ in guides], flush=True)


if __name__ == '__main__':
    main()
