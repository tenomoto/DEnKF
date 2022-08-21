import numpy as np
import os
import sys

imax, jmax = 129, 129
n = imax * jmax
nens = 25
nstep = 4
inflation_factor = 1.06
rsig = 4.0
exp = "cycle"
expdir = f"{os.path.dirname(__file__)}/{exp}"
fdir = f"{expdir}/forecast"
odir = f"{expdir}/obs"
adir = f"{expdir}/analysis-python"
lfile = f"{expdir}/loc.dat"
oofile = f"{odir}/obsoff.txt"


def analysis(t):
    loc = np.fromfile(lfile, dtype=np.float32).reshape(n, n)
    xfq = np.fromfile(f"{fdir}/q{t:06d}.dat", dtype=np.float32).reshape(nens, n).T
    xfp = np.fromfile(f"{fdir}/p{t:06d}.dat", dtype=np.float32).reshape(nens, n).T
    xfqbar = xfq.mean(axis=1)
    xfpbar = xfp.mean(axis=1)
    pf12q = (xfq - xfqbar[:, None]) / np.sqrt(nens - 1.0)
    pf12p = (xfp - xfpbar[:, None]) / np.sqrt(nens - 1.0)
    pfqp = pf12q @ pf12p.T
    pfpp = pf12p @ pf12p.T
    pfqp = inflation_factor * loc * pfqp
    pfpp = inflation_factor * loc * pfpp

    obsoff = np.loadtxt(oofile).astype(np.int32)
    y = np.fromfile(f"{odir}/o{t:06d}.dat", dtype=np.float32)
    nobs = y.size
    obsloc = np.arange(obsoff[t//nstep], n, n // nobs) - 1
    hmat = np.zeros([nobs, n])
    rmat = np.zeros([nobs, nobs])
    rmat[np.diag_indices_from(rmat)] = rsig

    for m in range(nobs):
        hmat[m, obsloc[m]] = 1.0
    hhphr_inv = hmat.T @ np.linalg.inv(hmat @ pfpp @ hmat.T + rmat)
    kmatq = pfqp @ hhphr_inv
    kmatp = pfpp @ hhphr_inv
    innov = y - hmat @ xfpbar
    xaqbar = xfqbar + kmatq @ innov
    xapbar = xfpbar + kmatp @ innov
    hmatpf12p = hmat @ pf12p
    pa12q = pf12q - 0.5 * kmatq @ hmatpf12p
    pa12p = pf12p - 0.5 * kmatp @ hmatpf12p
    xaq = xaqbar[:, None] + pa12q * np.sqrt(nens - 1.0)
    xap = xapbar[:, None] + pa12p * np.sqrt(nens - 1.0)
    xaq.T.astype(np.float32).tofile(f"{adir}/q{t:06d}.dat")
    xap.T.astype(np.float32).tofile(f"{adir}/p{t:06d}.dat")
    xapbar.astype(np.float32).tofile(f"{adir}/pbar{t:06d}.dat")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage :: python {sys.argv[0]} step")
        sys.exit()
    t = int(sys.argv[1])
    analysis(t)
