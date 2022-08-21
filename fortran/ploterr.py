import numpy as np
import matplotlib.pyplot as plt
import sys

nstep = 4
nens = 25
ncycle = 300
imax, jmax = 129, 129
n = imax * jmax
rstd = 2.0
exp = "cycle"

l2 = np.zeros(ncycle+1)
sd = np.zeros(ncycle+1)
for i in range(ncycle+1):
    step = i * nstep
    true = np.fromfile(f"{exp}/true/p{step:06d}.dat", dtype=np.float32).reshape(jmax, imax)
    psia = np.fromfile(f"{exp}/analysis/p{step:06d}.dat", dtype=np.float32).reshape(nens, jmax, imax)
    psiabar = psia.mean(axis=0)
    l2[i] = np.linalg.norm((psiabar - true)) / (imax - 1.0)
    sd[i] = psia.std(axis=0).mean()

plt.rcParams["font.size"] = 18
fig, ax = plt.subplots(figsize=[10, 7])
ax.plot(l2, label="analysis error")
ax.plot(sd, label="analysis stddev")
ax.axhline(rstd, c="gray", ls="--", label="obs stddev")
ax.legend()
ax.set_xlabel("cycle")
ax.set_ylabel(r"$\ell_2$")
fig.savefig("l2.png", bbox_inches="tight", dpi=300)
plt.show()
