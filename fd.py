import numpy as np


def laplacian(p):
    lp = np.zeros_like(p)
    lp[1:-1, 1:-1] = p[2:, 1:-1] + p[:-2, 1:-1] + p[1:-1, 2:] + p[1:-1, :-2] \
                   - 4 * p[1:-1, 1:-1]
    return lp

def jacobian(p, q):
    ja = np.zeros_like(p)
    j1 = (q[1:-1, 2:] - q[1:-1, :-2]) * (p[2:, 1:-1] - p[:-2, 1:-1]) - \
         (q[2:, 1:-1] - q[:-2, 1:-1]) * (p[1:-1, 2:] - p[1:-1, :-2])
    j2 = (q[2:, 2:]  - q[2:, :-2])  * p[2:, 1:-1] - \
         (q[:-2, 2:] - q[:-2, :-2]) * p[:-2, 1:-1] - \
         (q[2:, 2:]  - q[:-2, 2:])  * p[1:-1, 2:] + \
         (q[2:,:-2]  - q[:-2, :-2]) * p[1:-1, :-2]
    j3 = q[1:-1, 2:]  * (p[2:, 2:]  - p[:-2, 2:]) - \
         q[1:-1, :-2] * (p[2:, :-2] - p[:-2, :-2]) - \
         q[2:, 1:-1]  * (p[2:, 2:]  - p[2:, :-2]) + \
         q[:-2, 1:-1] * (p[:-2, 2:] - p[:-2, :-2])
    ja[1:-1, 1:-1] += (j1 + j2 + j3) / 12
    return ja

def l2norm(u, h):
    return np.sqrt(h ** u.ndim * (u ** 2).sum())


