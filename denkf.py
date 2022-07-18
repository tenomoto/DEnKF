import numpy as np
import scipy.linalg as la
import sys


def calc_kalman_gain_mform(pf12, ymat, rmat):
    return pf12 @ ymat.T @ la.inv(ymat @ ymat.T + rmat)

def calc_kalman_gain_nform(pf12, ymat, rinv):
    ic = ymat.T @ rinv @ ymat
    ic[np.diag_indices_from(ic)] += 1.0
    return pf12 @ la.inv(ic) @ ymat.T @ rinv

def update_state(xf, kalman_gain, o_minus_b):
    return xf + kalman_gain @ o_minus_b

def update_anomalies(pf12, kalman_gain, ymat):
    return pf12 - 0.5 * kalman_gain @ ymat

def gen_obs_covariance(rin, nobs, invert=False):
    if type(rin) == float:
        if invert: rin = 1.0 / rin
        rout = np.diag(np.full(nobs, rin))
    elif type(rin) == np.ndarray:
        if rin.dim == 1:
            if invert: rin = 1.0 / rin
            rout = np.diag(rin)
        else:
            if invert: rin = la.inv(rin)
            rout = rin
    return rout

def analysis(xf, o_minus_b, pf12, ymat, r):
    nstate, nens = pf12.shape
    nobs = o_minus_b.size
    if nobs > nens:
        rmat = gen_obs_covariance(r, nobs)
        kalman_gain = calc_kalman_gain_mform(pf12, ymat, rmat)
    else:
        rinv = gen_obs_covariance(r, nobs, invert=True)
        kalman_gain = calc_kalman_gain_nform(pf12, ymat, rinv)
    xa = update_state(xf, kalman_gain, o_minus_b)
    pa12 = update_anomalies(pf12, kalman_gain, ymat)
    return xa, pa12
