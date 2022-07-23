import numpy as np


def jacobian(p, q, d):
    j1 = ((q[2:,:] - q[:-2, :]) * (p[:, 2] - p[:, :-2])
            -(q[:, 2:] - q[:, :-2]) * (p[2:, :] - p[:-2, :]))
    return (j1 + j2 + j3) / (12 * d**2)

