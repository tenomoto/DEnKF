import numpy as np
import matplotlib.pyplot as plt
import sys

nstep = 4
nens = 25
ncycle = [300, 300, 300, 100]
imax, jmax = 129, 129
n = imax * jmax
rstd = 2.0

l2 = np.zeros(max(ncycle)+1)
sd = np.zeros(max(ncycle)+1)

plt.rcParams["font.size"] = 18
fig, ax = plt.subplots(figsize=[10, 7])

exp = ["qp", "q", "p", "noloc"]
lab = [r"$q\psi$", r"$q$", r"$\psi$", r"$q\psi$ noloc"]
wid = [10, 3, 6, 6]
col = ["tab:blue", "tab:orange", "tab:red", "tab:green"]
for j in range(len(exp)):
    for i in range(ncycle[j]+1):
        step = i * nstep
        true = np.fromfile(f"cycle_{exp[j]}/true/p{step:06d}.dat",
                dtype=np.float32).reshape(jmax, imax)
        psia = np.fromfile(f"cycle_{exp[j]}/analysis/p{step:06d}.dat",
                dtype=np.float32).reshape(nens, jmax, imax)
        psiabar = psia.mean(axis=0)
        l2[i] = np.linalg.norm((psiabar - true)) / (imax - 1.0)
        sd[i] = psia.std(axis=0).mean()
    print(f"{exp[j]}: {l2[:10]}")
    ax.semilogy(l2[:ncycle[j]+1], c=col[j], lw=wid[j], label=f"{lab[j]}")
    ax.semilogy(sd[:ncycle[j]+1], c=col[j], lw=wid[j], ls=":")
ax.axhline(rstd, c="gray", ls="--", label="obs stddev")
ax.legend(fontsize=14)
ax.set_xlabel("cycle")
ax.set_ylabel(r"$\ell_2$")
#ax.set_ylim([-0.1, 6.0])
fig.savefig("l2_compare.png", bbox_inches="tight", dpi=300)
fig.savefig("l2_compare.pdf", bbox_inches="tight")
plt.show()
