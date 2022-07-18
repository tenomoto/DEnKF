import numpy as np
import matplotlib.pyplot as plt


plt.rcParams["font.size"] = 18
for t in range(0, 500, 50):
    xtrue = np.load(f"xtrue{t:04d}.npy")
    xf = np.load(f"xf{t:04d}.npy")
    xa = np.load(f"xa{t:04d}.npy")
    y = np.load(f"y{t:04d}.npy")
    n = xf.size
    nobs = y.size
    obs_int = n // nobs
    obs_loc = np.arange(obs_int // 2 - 1, n, obs_int)

    fig, ax = plt.subplots(figsize=[10, 5])
    ax.plot(xtrue, label="true")
    ax.plot(xf, label="guess")
    ax.plot(xa, label="analysis")
    ax.scatter(obs_loc, y)
    ax.set_title(f"t={t}")
    ax.legend()
    fig.savefig(f"x{t:04d}.png", bbox_inches="tight", dpi=300)
#plt.show()
