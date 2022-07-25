import numpy as np
import matplotlib.pyplot as plt
from fd import four_point_sum, laplacian, l2norm
import sys


def four_point_sum(p):
    return p[2:, 1:-1] + p[:-2, 1:-1] + p[1:-1, 2:] + p[1:-1, :-2]

def laplacian(p):
    lp = np.zeros_like(p)
    lp[1:-1, 1:-1] = four_point_sum(p) - 4 * p[1:-1, 1:-1]
    return lp

def l2norm(u, h):
    return np.sqrt(h ** u.ndim * (u ** 2).sum())

def jacobi_step(pin, q, d, f, niter=100, tol=1e-5):
    d2 = d ** 2
    p = pin.copy()
    r = np.zeros_like(pin)
    res = 0
    k = 0
    for i in range(niter):
        k += 1
        r = (laplacian(p) / d2 - f * p - q)[1:-1, 1:-1]
#        p[1:-1, 1:-1] = (four_point_sum(p) - d2 * q[1:-1, 1:-1]) / (4 + d2 * f)
        p[1:-1, 1:-1] += d2 * r / (4 + d2 * f)
        res = l2norm(r, d)
        if res < tol: break
    return p, res, k

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

def v_cycle(p0, q, d, f, itermax=(10, 10, 10), tol=1.0e-5, debug=False):
    n = p0.shape[0]
    nlev = int(np.log2(n - 1))
    p, res, niter = jacobi_step(p0, q, d, f, itermax[0], tol)
    if debug: print(f"pre: res={res}, niter={niter}")
    qlist = [q]
    h = d
    for i in range(1, nlev):
        p = restrict(p)
        qlist.append(restrict(qlist[i-1]))
        h *= 2
#        if debug: print(f"down {i} {p.shape} {qlist[i].shape}")
        p, res, niter = jacobi_step(p, qlist[i], h, f, itermax[1], tol)
        if debug: print(f"restrict {i}: res={res}, niter={niter}")
    for i in range(nlev-2, -1, -1):
        p = prolong(p)
        h /= 2
#        if debug: print(f"up {i} {p.shape} {qlist[i].shape}")
        p, res, niter = jacobi_step(p, qlist[i], h, f, itermax[2], tol)
        if debug: print(f"prolong {i}: res={res}, niter={niter}")
    return p, res

def prolong_test():
    n = 5
    p = np.zeros([n, n])
    x = np.linspace(-1, 1, n)
    y = x
    for i in range(1,n-1):
        for j in range(1,n-1):
            p[i, j] = np.exp(-(x[i]**2 + y[j]**2))
    plt.matshow(p)
    plt.show()
    p1 = prolong(p)
    plt.matshow(p1)
    plt.show()

def mg_test(plot=True, itermax=(1, 1, 20000), tol=1e-5, debug=False):
    n = 129 
    ncycle = 1
    d = 1.0 / (n - 1)
    x = np.linspace(0, 1, n)
    y = np.linspace(0, 1, n)
    X, Y = np.meshgrid(x, y, indexing="ij")

    ptrue = (X**2 - X**4) * (Y**4 - Y**2)
    q = -2 * ((1 - 6 * X*X) * Y*Y * (1 - Y*Y) + \
              (1 - 6 * Y*Y) * X*X * (1 - X*X))

    p0 = np.zeros_like(q)
    p, res = v_cycle(p0, q, d, 0.0, itermax, tol, debug)
    err = l2norm(p - ptrue, d)
    print(f"res={res} err={err}")
    for i in range(1, ncycle):
        p, res = v_cycle(p, q, d, 0.0, itermax, tol, debug)
        err = l2norm(p - ptrue, d)
        print(f"cycle={i} res={res} err={err}")

    if plot:
        plt.rcParams["font.size"] = 12
        fig, axs = plt.subplots(1, 3, figsize=[14, 4])
        z = [ptrue, p, ptrue-p]
        title = [r"$(x^2-x^4)(y^4-y^2)$",
                 f"res={res:.2e} l2={err:.2e}",
                 r"mgrid$-$true"]
        for j in range(len(z)):
            ax = axs[j]
            if j < 2:
                c = ax.contourf(x, y, z[j], levels=np.linspace(-0.07, 0.0, 8))
            else:
                c = ax.contourf(x, y, z[j], cmap="coolwarm", vmin=-5e-6, vmax=5e-6)
            ax.set_title(title[j])
            ax.set_aspect("equal")
            fig.colorbar(c, ax=ax)
        fig.suptitle(f"Multigrid Jacobi itermax={itermax}")
        fig.tight_layout()
        fig.savefig("multigrid_jacobi.png", bbox_inches="tight", dpi=300)
        plt.show()

if __name__ == "__main__":
#    prolong_test()
    mg_test(plot=False, itermax=(1,1,10000), tol=1.0e-4, debug=True)
