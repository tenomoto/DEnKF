import numpy as np
import matplotlib.pyplot as plt
import sys


if len(sys.argv) < 2:
    print(f"Usage :: python {sys.argv[0]} exp")
    sys.exit()
exp = sys.argv[1]

steps = [300, 600, 900, 1200]
dstep = 4
nens = 25
nobs = 300
imax, jmax = 129, 129
n = imax * jmax
dt = 1.25

dobs = n // nobs
obsoff = np.loadtxt(f"{exp}/obs/obsoff.txt").astype(np.int32)
obsloc = np.zeros(nobs, dtype=int)
x = np.linspace(0, 1, imax)
y = np.linspace(0, 1, jmax)

vmin = [-30.0, -1.5, 0.0, -0.3]
vmax = [ 30.0,  1.5, 0.5,  0.3]
cmap = [ "RdYlBu_r", "RdYlBu_r", "Oranges", "RdYlBu_r"]
label = ["true obs", "analysis-true", "analysis std", "increment"]

fig, axs = plt.subplots(len(cmap), len(steps), figsize=[8, 8], constrained_layout=True)
pc = []
for step, j in zip(steps, np.arange(len(steps))):
    true = np.fromfile(f"{exp}/true/p{step:06d}.dat", dtype=np.float32).reshape(jmax, imax)
    psia = np.fromfile(f"{exp}/analysis/p{step:06d}.dat", dtype=np.float32).reshape(nens, jmax, imax)
    pabar = psia.mean(axis=0)
    pastd = psia.std(axis=0)
    psif = np.fromfile(f"{exp}/forecast/p{step:06d}.dat", dtype=np.float32).reshape(nens, jmax, imax)
    pfbar = psif.mean(axis=0)

    obs = np.fromfile(f"{exp}/obs/o{step:06d}.dat", dtype=np.float32)
    for i in range(nobs):
        obsloc[i] = obsoff[step//dstep] + i * dobs
    iloc = obsloc % imax
    jloc = obsloc // jmax

    z = [true, pabar-true, pastd, pabar-pfbar]
    for i in range(len(z)):
        ax = axs[j, i]
        pc.append(ax.pcolormesh(x, y, z[i], vmin=vmin[i], vmax=vmax[i], cmap=cmap[i]))
        ax.set_aspect("equal")
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        if j == 0:
            ax.set_title(label[i])
    axs[j, 0].scatter(x[iloc], y[jloc], c=obs, vmin=vmin[0], vmax=vmax[0],
            s=8, linewidths=0.5, edgecolors="gray", cmap=cmap[0])
    axs[j, 0].set_ylabel(r"$t=$"+f"{step*dt:<4.0f}")
for i in range(len(cmap)):
    fig.colorbar(pc[i], ax=axs[:, i], location="bottom",
            shrink=0.9, pad=0.01)
#fig.suptitle(f"{exp}")
fig.savefig(f"diff_{exp}.png", bbox_inches="tight", dpi=300)
fig.savefig(f"diff_{exp}.pdf", bbox_inches="tight")
plt.show()
