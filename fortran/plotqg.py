import numpy as np
import matplotlib.pyplot as plt


t = 1000
n = 129
q = np.fromfile(f"q{t:06d}.dat", dtype=np.float32).reshape(n, n)
psi = np.fromfile(f"p{t:06d}.dat", dtype=np.float32).reshape(n, n)
x = np.linspace(0, 1, n)
y = np.linspace(0, 1, n)
plt.rcParams["font.size"] = 18
fig, axs = plt.subplots(1, 2, figsize=[14, 6])
z = [q, psi]
#zmax = [1.5e6, 6.0e2]
zmax = [1.0e5, 1.0e1]
title = ["pv", r"$\psi$"]
for i in range(len(axs)):
    ax = axs[i]
    c = ax.pcolormesh(x, y, z[i])
#    c = ax.pcolormesh(x, y, z[i],
#            cmap="RdYlBu_r", vmin=-zmax[i], vmax=zmax[i])
    ax.set_title(title[i])
    ax.set_aspect("equal")
    fig.colorbar(c, ax=ax, shrink=0.8)
fig.suptitle(r"$t=$"+f"{t}")
#fig.savefig(f"pq{t:06d}.png", bbox_inches="tight", dpi=300)
plt.show()

