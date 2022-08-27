import numpy as np
import matplotlib.pyplot as plt


t = 10000
dt = 1.5
n = 129
q = np.fromfile(f"q{t:06d}.dat", dtype=np.float32).reshape(n, n)
psi = np.fromfile(f"p{t:06d}.dat", dtype=np.float32).reshape(n, n)
x = np.linspace(0, 1, n)
y = np.linspace(0, 1, n)
plt.rcParams["font.size"] = 24
fig, axs = plt.subplots(1, 2, figsize=[14, 6])
z = [q, psi]
zmax = [1.0e5, 3.0e1]
#zmax = [4.0e4, 1.0e1]
title = [r"$q$", r"$\psi$"]
for i in range(len(axs)):
    ax = axs[i]
#    c = ax.pcolormesh(x, y, z[i], cmap="RdYlBu_r")
    c = ax.pcolormesh(x, y, z[i],
            cmap="RdYlBu_r", vmin=-zmax[i], vmax=zmax[i])
    ax.set_title(title[i])
    ax.set_aspect("equal")
    fig.colorbar(c, ax=ax, shrink=0.8)
fig.suptitle(r"$t=$"+f"{t*dt:<5.0f}")
fig.tight_layout()
#fig.savefig(f"pq{t:06d}.png", bbox_inches="tight", dpi=300)
fig.savefig(f"pq{t:06d}.pdf", bbox_inches="tight")
plt.show()

