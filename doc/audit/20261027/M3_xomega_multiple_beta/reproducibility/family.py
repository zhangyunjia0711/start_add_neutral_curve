"""Original polynomial neutral equations in real (station, alpha, omega)."""
from functools import lru_cache
import numpy as np
from scipy.sparse import bmat, csc_matrix
from scipy.sparse.linalg import splu
from common import CASES
from solver import Pencil, qep, mac

SCALE = np.array([50., 20., 10.])


class Family:
    def __init__(self, case, beta):
        self.case, self.beta = case, float(beta)
        self.count = CASES[case]['stations']
        self.x = np.array(CASES[case]['x_mm'])
        self.native = lru_cache(maxsize=8)(self._native)

    def _native(self, i):
        p = Pencil(self.case, i, 0.)
        return p.spatial(self.beta), p.G

    def pencil(self, xi, side=1):
        if abs(xi - round(xi)) < 1e-9:
            xi = float(round(xi))
        assert 1 <= xi <= self.count, ('Domain boundary', xi)
        i = min(int(np.floor(xi)), self.count - 1)
        if side < 0 and abs(xi - round(xi)) < 1e-9:
            i = max(1, int(round(xi)) - 1)
        theta = xi - i
        a, ga = self.native(i)
        b, gb = self.native(i + 1)
        A = [(1 - theta) * x + theta * y for x, y in zip(a, b)]
        G = (1 - theta) * ga + theta * gb
        return A, G, [y - x for x, y in zip(a, b)], gb - ga, i, theta

    def evaluate(self, z, qref, pref, side=1):
        xi, alpha, omega = z * SCALE
        A, G, dA, dG, i, theta = self.pencil(float(xi), side)
        Q = A[0] + alpha * A[1] + alpha * alpha * A[2] + omega * G
        D = [dA[0] + alpha * dA[1] + alpha * alpha * dA[2] + omega * dG,
             A[1] + 2 * alpha * A[2], G]
        n = len(qref)
        M = bmat([[Q, csc_matrix(pref[:, None])],
                  [csc_matrix(qref.conj()[None, :]), csc_matrix((1, 1), dtype=complex)]], format='csc')
        lu = splu(M)
        rhs = np.r_[np.zeros(n, dtype=complex), 1. + 0j]
        solution, adjoint = lu.solve(rhs), lu.solve(rhs, trans='H')
        q, p, error = solution[:-1], adjoint[:-1], solution[-1]
        gradient = np.array([-np.vdot(p, matrix @ q) * scale for matrix, scale in zip(D, SCALE)])
        J = np.array([gradient.real, gradient.imag])
        _, singular, V = np.linalg.svd(J, full_matrices=True)
        assert singular[-1] > 1e-10 * singular[0], ('Neutral Jacobian rank loss', singular, z)
        tangent = V[-1]
        q, p = q / np.linalg.norm(q), p / np.linalg.norm(p)
        result = dict(xi=float(xi), x_mm=float((1 - theta) * self.x[i - 1] + theta * self.x[i]),
                      alpha_r=float(alpha), alpha_i=0., beta_r=self.beta, beta_i=0.,
                      omega_r=float(omega), omega_i=0., growth=0.,
                      right_residual=float(np.linalg.norm(Q @ q)), left_residual=float(np.linalg.norm(Q.conj().T @ p)),
                      raw_MAC=mac(qref, q), border_error_r=float(error.real), border_error_i=float(error.imag),
                      native_bracket=[i, i + 1], theta=float(theta))
        return error, J, tangent, q, p, result

    def correct(self, pred, plane, q, p, maxit=25):
        z = pred.copy()
        constrained = np.linalg.norm(plane[1:]) > 1e-12
        xi = pred[0] * SCALE[0]
        cell = min(int(np.floor(xi)), self.count - 1)
        if abs(xi - round(xi)) < 1e-9 and plane[0] < 0:
            cell = max(1, int(round(xi)) - 1)
        lo, hi = cell / SCALE[0], (cell + 1) / SCALE[0]
        trace = []
        for it in range(maxit):
            side = -1 if constrained and abs(z[0] * SCALE[0] - (cell + 1)) < 1e-8 else 1
            e, J, t, v, w, r = self.evaluate(z, q, p, side)
            F = np.array([e.real, e.imag, np.dot(z - pred, plane)])
            trace.append(dict(point=r, F=F.tolist()))
            if abs(e) < 1e-12 and abs(F[2]) < 1e-9 and r['right_residual'] < 1e-8:
                return z, r, v, w, t, trace
            A = np.vstack([J, plane])
            weights = 1 / np.maximum(np.linalg.norm(A, axis=1), 1e-30)
            dz = np.linalg.solve(A * weights[:, None], -F * weights)
            if np.linalg.norm(dz) > .2:
                dz *= .2 / np.linalg.norm(dz)
            # The operator is continuous across native cells. A Newton iterate
            # may cross the predictor's cell at a real x fold; freezing its cell
            # can cause repeated reversal rather than following the contour.
            z += dz
        raise RuntimeError('Neutral corrector iteration limit')

    def seed(self, source, q):
        i = int(round(source['xi']))
        original = Pencil(self.case, i, source['omega_r'])
        _, q, _, p = qep(original.spatial(source['beta_r']), source['alpha_r'], q, True)
        z = np.array([i, source['alpha_r'], source['omega_r']]) / SCALE
        return self.correct(z, np.array([1., 0., 0.]), q, p)
