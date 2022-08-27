import numpy as np
import matplotlib.pyplot as plt
import sys


if len(sys.argv) < 3:
    print(f"Usage :: python {sys.argv[0]} step exp")
    sys.exit()
step = int(sys.argv[1])
exp = sys.argv[2]
dstep = 4
nens = 25
imax, jmax = 129, 129
n = imax * jmax
dt = 1.25

true = np.fromfile(f"{exp}/true/p{step:06d}.dat", dtype=np.float32).reshape(jmax, imax)
psia = np.fromfile(f"{exp}/analysis/p{step:06d}.dat", dtype=np.float32).reshape(nens, jmax, imax)
pabar = psia.mean(axis=0)
pastd = psia.std(axis=0)
psif = np.fromfile(f"{exp}/forecast/p{step:06d}.dat", dtype=np.float32).reshape(nens, jmax, imax)
pfbar = psif.mean(axis=0)

obsoff = np.loadtxt(f"{exp}/obs/obsoff.txt").astype(np.int32)
obs = np.fromfile(f"{exp}/obs/o{step:06d}.dat", dtype=np.float32)
dobs = n // obs.size
obsloc = np.zeros(obs.size, dtype=int)
for i in range(obs.size):
    obsloc[i] = obsoff[step//dstep] + i * dobs
iloc = obsloc % imax
jloc = obsloc // jmax

x = np.linspace(0, 1, imax)
y = np.linspace(0, 1, jmax)

fig, axs = plt.subplots(2, 3, figsize=[12, 8])
z = [true, pabar, pfbar, pabar-true, pastd, pabar-pfbar, ]
label = ["true", "analysis", "forecast", "analysis-true", "analysis std", "increment"]
if step >= 100:
    vmin = [-30.0, -30.0, -30.0, -1.5, 0.0, -0.3]
    vmax = [ 30.0,  30.0,  30.0,  1.5, 0.5,  0.3]
elif step >= 20:
    vmin = [-30.0, -30.0, -30.0, -3.0, 0.0, -3.0]
    vmax = [ 30.0,  30.0,  30.0,  3.0, 2.0,  3.0]
else:
    vmin = [-30.0, -30.0, -30.0, -5.0, 0.0, -10.0]
    vmax = [ 30.0,  30.0,  30.0,  5.0, 5.0,  10.0]
cmap = [ "RdYlBu_r", "RdYlBu_r", "RdYlBu_r", "RdYlBu_r", "Oranges", "RdYlBu_r"]
for i in range(len(z)):
    ax = axs[i//3, i % 3]
    pc = ax.pcolormesh(x, y, z[i], vmin=vmin[i], vmax=vmax[i], cmap=cmap[i])
    ax.set_aspect("equal")
    ax.set_title(label[i])
    fig.colorbar(pc, ax=ax, shrink=0.7)
axs[0, 0].scatter(x[iloc], y[jloc], c=obs, vmin=vmin[0], vmax=vmax[0],
        linewidths=0.5, edgecolors="white", cmap=cmap[0])
fig.suptitle(f"{exp} cycle {step // dstep} step {step} t={step * dt:<4.0f}",
        fontsize=18)
fig.tight_layout()
fig.savefig(f"diff{step}_{exp}.png", bbox_inches="tight", dpi=300)
plt.show()
