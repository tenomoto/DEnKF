from math import tau
import numpy as np
from fd import laplacian, jacobian
from mg import v_cycle
from ode import rk4
import sys


def step(q, psi, y, beta, f, eps, a, tau0, itermax, tol):
    d = y[1] - y[0]
    d2 = d * 2
    dd = d ** 2
    psi[:,:], _ = v_cycle(psi, q, d, f, itermax, tol)
    psix = np.zeros_like(psi)
    dqdt = np.zeros_like(psi)
    psix[1:-1, 1:-1] = (psi[2:, 1:-1] - psi[:-2, 1:-1]) / d2
    lap3psi = laplacian(laplacian(q + f * psi) / dd) / dd 
    jac = jacobian(psi, q) / dd
    dqdt[1:-1, 1:-1] = (-beta * psix - eps * jac - a * lap3psi + tau0 * np.sin(tau * y[None,]))[1:-1, 1:-1]
    return dqdt

if __name__ == "__main__":
    n = 129 
    dt = 1.5
    tsave = 0 
    nstep = 1000
    nsave = 100
    itermax = 1, 1, 100
    datadir="free"
    x = np.linspace(0.0, 1.0, n)
    y = np.linspace(0.0, 1.0, n)
    q = np.zeros([n, n])
    seed = 514
    rng = np.random.default_rng(514)
    q[1:-1, 1:-1] = 2.0 * rng.random([n-2, n-2]) - 1.0
    psi = np.zeros([n, n])
    beta, f, eps, a, tau0 = 0.0, 0.0, 1.0, 2.0e-12, 0.0
#    beta, f, eps, a, tau0 = 1.0, 1600, 1.0e-5, 2.0e-12, -tau
    tol = 1.0e-4
    params = beta, f, eps, a, tau0, itermax, tol
    for i in range(nstep+1): 
        if i >= tsave and i % nsave == 0:
            np.save(f"q{i:06d}.npy", q)
            np.save(f"p{i:06d}.npy", psi)
        print(f"step {i} p: min={psi.min():5.2e} max={psi.max():5.2e} q: min={q.min():5.2e} max={q.max():5.2e}")
        dq = np.zeros([n, n])
        dq[1:-1, 1:-1] = rk4(step, q, dt, psi, y, *params)[1:-1, 1:-1]
        np.save(f"dq.npy", dq)
        q[1:-1, 1:-1] += rk4(step, q, dt, psi, y, *params)[1:-1, 1:-1]
