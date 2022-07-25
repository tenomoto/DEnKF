from math import tau
import numpy as np
from fd import laplacian, jacobian
from mg import v_cycle
from ode import rk4


def step(q, psi, y, f=1600, eps=1.0e-5, a=2.0e-12, itermax=(1, 1, 5000), tol=1.0e-4):
    d = y[1] - y[0]
    d2 = d ** 2
    psi, _ = v_cycle(psi, q, d, f, itermax, tol)
    psix = np.zeros_like(psi)
    psix[1:-1, 1:-1] = (psi[2:, 1:-1] - psi[:-2, 1:-1]) / d
    lap3psi = laplacian(laplacian(q + f * psi) / d2) / d2 
    jac = jacobian(psi, q) / d2
    dqdt = -psix - eps * jac - a * lap3psi + tau * np.sin(tau * y[None,])
    return dqdt

if __name__ == "__main__":
    n = 129 
    dt = 1.5
    nstep = 100
    nsave = 10
    x = np.linspace(0.0, 1.0, n)
    y = np.linspace(0.0, 1.0, n)
    q = np.zeros([n, n])
    psi = np.zeros([n, n])
    dqdt = step(q, psi, y)
    params = psi, y
    for i in range(nstep+1): 
        print(f"step {i} qmin={q.min()} qmax={q.max()}")
        q[1:-1, 1:-1] += rk4(step, q, dt, *params)[1:-1, 1:-1]
        if i % nsave == 0:
            np.save(f"q{i:05d}.npy", q)

    
