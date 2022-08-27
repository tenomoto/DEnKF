import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("..")
from fd import l2norm


n = 129
d = 1 / (n - 1)
x = np.linspace(0, 1, n)
y = np.linspace(0, 1, n)
X, Y = np.meshgrid(x, y, indexing="ij")
ptrue = (X**2 - X**4) * (Y**4 - Y**2)

buf = np.fromfile("mg.dat")
nlev = 7
off = 0
m = 129
for i in range(nlev):
    off += m * m
    m = m // 2 + 1
p = buf[off:off+n*n].reshape(n, n)

plt.rcParams["font.size"] = 12
fig, axs = plt.subplots(1, 3, figsize=[14, 4])
z = [ptrue, p, ptrue-p]
err = l2norm(p - ptrue, d)
title = [r"true $(x^2-x^4)(y^4-y^2)$",
         f"mgrid l2={err:.2e}",
         r"mgrid$-$true"]
for j in range(len(z)):
    ax = axs[j]
    if j < 2:
#        c = ax.contourf(x, y, z[j], levels=np.linspace(-0.07, 0.0, 8))
        c = ax.matshow(z[j])
    else:
#        c = ax.contourf(x, y, z[j], cmap="coolwarm", vmin=-5e-6, vmax=5e-6)
        c = ax.matshow(z[j], cmap="coolwarm", vmin=-5e-6, vmax=5e-6)
    ax.set_title(title[j])
    ax.set_aspect("equal")
    if j > 0:
        fig.colorbar(c, ax=ax)
#fig.suptitle(f"Multigrid Jacobi itermax={itermax}")
fig.tight_layout()
fig.savefig("multigrid_jacobi.pdf", bbox_inches="tight")
#plt.show()

