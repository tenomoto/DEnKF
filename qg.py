from math import tau
import numpy as np
import mg


def laplacian(p):
    lp = np.zeros_like(p)
    lp[1:-1, 1:-1] = mg.four_point_sum(p) - 4 * p[1:-1, 1:-1]
    return lp

def jacobian(p, q):
    ja = np.zeros_like(p)
    j1 = (q[1:-1, 2:] - q[1:-1, :-2]) * (p[2:, 1:-1] - p[:-2, 1:-1]) - \
         (q[2:, 1:-1] - q[:-2, 1:-1]) * (p[1:-1, 2:] - p[1:-1, :-2]))
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

def solve_helmholz(q, psi_old, d, f):
    psi, res, = mg.v_cycle(psi_old, q, 
    pass

def step(q, psi_old, y, f=1600, eps=1.0e-5, a=2.0e-12, itermax=(4, 20, 20, 100)):
    d = y[1] - y[0]
    d2 = d ** 2
    psi = solve_helmholz(psi_old, q, d, f, itermax)
    psix = psi[2:, :] - psi[:2, :]
    lap3psi = laplacian(laplacian(q + f * psi) / d2)) / d2 
    jac = jacobian(psi, q) / d2
    dqdt = -psix - f * jac - a * lap3psi + tau * np.sin(tau * y[None,])
    return dqdt, psi
    
