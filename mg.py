import numpy as np
import matplotlib.pyplot as plt
import sys


def four_point_sum(p):
    return p[2:, 1:-1] + p[:-2, 1:-1] + p[1:-1, 2:] + p[1:-1, :-2]

def l2norm(u, h):
    return np.sqrt(h ** u.ndim * (u ** 2).sum())

def jacobi_step(pin, q, d, f, niter=100, tol=1e-5):
    d2 = d ** 2
    p = pin.copy()
    r = np.zeros_like(pin)
    res = 0
    for i in range(niter):
        p[1:-1, 1:-1] = (four_point_sum(p) - d2 * q[1:-1, 1:-1]) / (4 + d2 * f)
        r[1:-1, 1:-1] = q[1:-1, 1:-1] - (
                four_point_sum(p) - 4*p[1:-1, 1:-1]) / d2 + f * p[1:-1, 1:-1]  
        res = l2norm(r, d)
        if res < tol: break
    return p, res, i

def restrict(pin):
    p = pin[::2, ::2].copy()
    return p

def prolong(pin):
    imax, jmax = pin.shape
    imax, jmax = 2 * imax - 1, 2 * jmax - 1
    p = np.zeros([imax, jmax])
# copy values at matching elements
    p[2:-2:2, 2:-2:2] = pin[1:-1, 1:-1]
# one-dimensional interpolation
    p[2:imax-2:2, 1:-1:2] = 0.5 * (p[2:imax-2:2, 2::2] + p[2:imax-1:2, :-2:2])
    p[1:-1:2, 2:jmax-2:2] = 0.5 * (p[2::2, 2:jmax-2:2] + p[:-2:2, 2:jmax-2:2])
    p[1:imax-1:2, 1:-1:2] = 0.5 * (p[1:imax-1:2, 2::2] + p[1:imax-1:2, :-2:2])
    return p

def v_cycle(p0, q, d, f, itermax=(10, 10, 10), debug=False):
    p, res, niter = jacobi_step(p0, q, d, f, itermax[0])
    if debug: print(f"pre: res={res}, niter={niter}")
    qlist = []
    h = d
    for i in range(nlev):
        p = restrict(p)
        if i == 0:
            qlist.append(restrict(q))
        else:
            qlist.append(restrict(qlist[i-1]))
        h *= 2
        p, res, niter = jacobi_step(p, qlist[i], h, f, itermax[1])
        if debug: print(f"restrict {i}: res={res}, niter={niter}")
    for i in range(nlev):
        p = prolong(p)
        h /= 2
        p, res, niter = jacobi_step(p, qlist[nlev - 2 - i], h, 0, itermax[2])
        if debug: print(f"prolong {i}: res={res}, niter={niter}")
    p, res, niter = jacobi_step(p, q, d, f, itermax[3])
    if debug: print(f"post: res={res}, niter={niter}")
    return p, res

if __name__ == "__main__":
#    n = 5
#    p = np.zeros([n, n])
#    x = np.linspace(-1, 1, n)
#    y = x
#    for i in range(1,n-1):
#        for j in range(1,n-1):
#            p[i, j] = np.exp(-(x[i]**2 + y[j]**2))
#    plt.matshow(p)
#    plt.show()
#    p1 = prolong(p)
#    plt.matshow(p1)
#    plt.show()
    n = 129 
    nlev = int(np.log2(n - 1)) - 1
    itermax = 4, 20, 20, 100
    d = 1.0 / (n - 1)
    x = np.linspace(0, 1, n)
    y = np.linspace(0, 1, n)
    X, Y = np.meshgrid(x, y)

    ptrue = (X**2 - X**4) * (Y**4 - Y**2)
    q = -2 * ((1 - 6 * X*X) * Y*Y * (1- Y*Y) + \
            (1 - 6 * Y*Y) * X*X * (1- X*X))
    p0 = np.zeros_like(q)
    p, res = v_cycle(p0, q, d, 0.0, itermax)
    err = l2norm(p - ptrue, d)
    print(f"res={res} err={err}")

    plt.rcParams["font.size"] = 18
    fig, ax = plt.subplots(figsize=[10, 8])
#    c = ax.contourf(x, y, q)
#    c = ax.contourf(x, y, p)
#    ax.set_title(f"itermax={itermax} res={res:.2e} l2={err:.2e}")
#    c = ax.contourf(x, y, ptrue)
#    ax.set_title(r"$(x^2-x^4)(y^4-y^2)$")
    c = ax.contourf(x, y, p-ptrue, cmap="coolwarm", vmin=-3e-3, vmax=3e-3)
    ax.set_title(r"mgrid$-$true")
    fig.colorbar(c)
#    fig.savefig("ptrue.png", bbox_inches="tight", dpi=300)
#    fig.savefig("p.png", bbox_inches="tight", dpi=300)
    fig.savefig("p-ptrue.png", bbox_inches="tight", dpi=300)
    plt.show()
