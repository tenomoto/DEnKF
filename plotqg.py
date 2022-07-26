import numpy as np
import matplotlib.pyplot as plt


t = 1000
q = np.load(f"q{t:05d}.npy").T
psi = np.load(f"psi{t:05d}.npy").T
n = q.shape[0]
x = np.linspace(0, 1, n)
y = np.linspace(0, 1, n)
plt.rcParams["font.size"] = 18
fig, axs = plt.subplots(1, 2, figsize=[14, 6])
z = [q, psi]
labels = ["pv", r"$\psi$"]
for i in range(len(axs)):
    ax = axs[i]
    c = ax.contourf(x, y, z[i])
    ax.set_aspect("equal")
    fig.colorbar(c, ax=ax)
plt.show()

