import numpy as np
import matplotlib.pyplot as plt


t = 400
q = np.load(f"q{t:05d}.npy").T
psi = np.load(f"p{t:05d}.npy").T
n = q.shape[0]
x = np.linspace(0, 1, n)
y = np.linspace(0, 1, n)
plt.rcParams["font.size"] = 18
fig, axs = plt.subplots(1, 2, figsize=[14, 6])
z = [q, psi]
zmax = [0.6, 1.2e-3]
title = ["pv", r"$\psi$"]
for i in range(len(axs)):
    ax = axs[i]
    c = ax.pcolormesh(x, y, z[i],
            cmap="RdYlBu_r", vmin=-zmax[i], vmax=zmax[i])
    ax.set_title(title[i])
    ax.set_aspect("equal")
    fig.colorbar(c, ax=ax, shrink=0.8)
fig.suptitle(r"$t=$"+f"{t}")
plt.show()

