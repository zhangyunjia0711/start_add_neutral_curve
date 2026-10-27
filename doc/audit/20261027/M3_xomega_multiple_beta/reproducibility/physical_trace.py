"""Actual-profile continuation through the optimized nose's problematic cells."""
import json
from collections import OrderedDict
import numpy as np
from common import ROOT
from physical import Pair
from temporal import Family as Temporal, trace
from solver import export
from scipy.sparse.linalg import norm


class PhysicalFamily(Temporal):
    model = 'Guide evaluated on reassembled physical local profiles; finite-difference xi derivatives; no CFD mesh refinement'

    def __init__(self, case, beta, tag):
        super().__init__(case, beta)
        self.pairs, self.tag = {}, tag
        self.fixed_cache = OrderedDict()

    def physical(self, i, theta):
        if i not in self.pairs:
            self.pairs[i] = Pair(self.case, i, 0., self.tag)
        key = (i, float(theta).hex())
        if key in self.fixed_cache:
            self.fixed_cache.move_to_end(key)
            return self.fixed_cache[key]
        directory = self.pairs[i].build(float(theta))
        A, _ = export(directory, 5, self.beta, 0., label=self.tag + '_fixed_beta')
        B, _ = export(directory, 5, self.beta, 1., label=self.tag + '_fixed_beta')
        assert norm(A[1] - B[1]) < 1e-10 and norm(A[2] - B[2]) < 1e-10
        result = A, B[0] - A[0]
        self.fixed_cache[key] = result
        while len(self.fixed_cache) > 16:
            self.fixed_cache.popitem(last=False)
        return result

    def pencil(self, xi, side=1):
        assert 1 <= xi <= self.count
        if abs(xi - round(xi)) < 1e-9:
            xi = float(round(xi))
        i = min(int(np.floor(xi)), self.count - 1)
        if side < 0 and abs(xi - round(xi)) < 1e-9:
            i = max(1, int(round(xi)) - 1)
        theta = float(xi - i)
        h = 2e-5
        low, high = max(0., theta - h), min(1., theta + h)
        A, G = self.physical(i, theta)
        al, gl = self.physical(i, low)
        ah, gh = self.physical(i, high)
        return A, G, [(b - a) / (high - low) for a, b in zip(al, ah)], (gh - gl) / (high - low), i, theta


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--seed', type=int, default=28)
    parser.add_argument('--tag', default='rough_beta-20_physical_c0')
    args = parser.parse_args()
    seed = json.loads((ROOT / 'outputs/seeds.json').read_text())['results'][args.seed]
    trace(seed, args.tag, PhysicalFamily(seed['case'], seed['beta'], args.tag))


if __name__ == '__main__':
    main()
