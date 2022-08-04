from math import tau
import numpy as np
from fd import laplacian, jacobian
from mg import v_cycle
from ode import rk4


def step(q, psi, y, c=1.0, f=1600, eps=1.0e-5, a=2.0e-12, tau0=tau, itermax=(1, 1, 5000), tol=1.0e-4):
    d = y[1] - y[0]
    d2 = d * 2
    dd = d ** 2
    psi[:,:], _ = v_cycle(psi, q, d, f, itermax, tol)
    psix = np.zeros_like(psi)
    dqdt = np.zeros_like(psi)
    psix[1:-1, 1:-1] = (psi[2:, 1:-1] - psi[:-2, 1:-1]) / d2
    lap3psi = laplacian(laplacian(q + f * psi) / dd) / dd 
    jac = jacobian(psi, q) / dd
    dqdt[1:-1, 1:-1] = (-c * psix - eps * jac - a * lap3psi + tau0 * np.sin(tau * y[None,]))[1:-1, 1:-1]
    return dqdt

if __name__ == "__main__":
    n = 129 
    dt = 1.5
    nstep = 1000
    nsave = 100
    x = np.linspace(0.0, 1.0, n)
    y = np.linspace(0.0, 1.0, n)
    q = np.zeros([n, n])
#    seed = 514
#    rng = np.random.default_rng(514)
#    q = 2.0 * rng.random([n, n]) - 1.0
    psi = np.zeros([n, n])
    for i in range(nstep+1): 
        print(f"step {i} p: min={psi.min():5.2e} max={psi.max():5.2e} q: min={q.min():5.2e} max={q.max():5.2e}")
#        params = 0.0, 0.0, 1.0, 2.0e-12, 0.0
#        q[1:-1, 1:-1] += rk4(step, q, dt, psi, y, *params)[1:-1, 1:-1]
        q[1:-1, 1:-1] += rk4(step, q, dt, psi, y)[1:-1, 1:-1]
        if i % nsave == 0:
            np.save(f"q{i:05d}.npy", q)
            np.save(f"psi{i:05d}.npy", psi)

    
