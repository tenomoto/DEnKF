import math
import numpy as np
import denkf
import sys


dtobs = 5
nobs = 4
r = 0.01
nens = 55 
nstep = 500

def step(x):
    return np.roll(x, 1)

def gen_state(n, kmax, rng):
    s = np.zeros(n)
    phi0 = math.tau * np.arange(1, n+1) / n
    for k in range(kmax+1):
        ak = rng.random()
        phik = math.tau * rng.random()
        s += ak * np.sin(k * phi0 + phik)
    s /= s.std()
    return s

def gen_ens(nens, n, kmax, rng):
    xmat = np.zeros([n, nens])
    for i in range(nens):
        xmat[:, i] = gen_state(n, kmax, rng)
    xmean = xmat.mean(axis=1)
    return xmat - xmean[:, None]

def gen_obs(xtrue, rng, obs_loc, ystd):
    nobs = obs_loc.size
    y = xtrue[obs_loc] + ystd * rng.standard_normal(nobs)
    return y

def simulate_obs(x, obs_loc):
    return x[obs_loc]

def gen_ymat(pf12, obs_loc):
    n, nens = pf12.shape
    nobs = obs_loc.size
    ymat = np.zeros([nobs, nens])
    for j in range(nens):
        ymat[:, j] = simulate_obs(pf12[:, j], obs_loc)
    return ymat

if __name__ == "__main__":
    n = 1000
    kmax = 25
    nsave = nstep // dtobs
    seed = 514
    dtsave = dtobs * 10

    obs_int = n // nobs
    obs_loc = np.arange(obs_int // 2 - 1, n, obs_int)

    rng = np.random.default_rng(seed)

    pf12 = gen_ens(nens, n, kmax, rng) / np.sqrt(nens - 1)
    xf = gen_state(n, kmax, rng)
    xtrue = gen_state(n, kmax, rng) + xf

    sprd= np.zeros([nsave])
    l2 = np.zeros([nsave])
    i = 0
    for t in range(nstep):
        if t % dtobs == 0:
            y = gen_obs(xtrue, rng, obs_loc, np.sqrt(r))
            o_minus_b = y - simulate_obs(xf, obs_loc)
            ymat = gen_ymat(pf12, obs_loc)
            xa, pa12 = denkf.analysis(xf, o_minus_b, pf12, ymat, r)
            xens = xa[:, None] + pa12 * np.sqrt(n - 1)
            l2[i] = (xa - xtrue).std()
            sprd[i] = np.sqrt(np.diag(pa12 @ pa12.T)).mean()
            if t % dtsave == 0:
                np.save(f"xf{t:04d}.npy", xf)
                np.save(f"xa{t:04d}.npy", xa)
                np.save(f"y{t:04d}.npy", y)
                np.save(f"xtrue{t:04d}.npy", xtrue)
#            print(t, l2[i], sprd[i])
            i += 1
        xtrue = step(xtrue)
        for j in range(nens):
            xens[:, j] = step(xens[:, j])
        xf = xens.mean(axis=1)
        pf12 = (xens - xf[:, None]) / np.sqrt(n - 1)
    np.save("l2.npy", l2)
    np.save("sprd.npy", sprd)

