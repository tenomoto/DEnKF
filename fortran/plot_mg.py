import numpy as np
import matplotlib.pyplot as plt

nlev = 7
buf = np.fromfile("mg.dat")

plt.rcParams["font.size"] = 18
fig = plt.figure(figsize=[35, 10], tight_layout=True)
n = 129
off = 0
for i in range(nlev):
    ax = fig.add_subplot(2, nlev, i+1)
    ax.matshow(buf[off:off+n*n].reshape(n, n))
    ax.set_aspect("equal")
    off += n * n
    n = n // 2 + 1
n = 129
for i in range(nlev):
    ax = fig.add_subplot(2, nlev, i+nlev+1)
    ax.matshow(buf[off:off+n*n].reshape(n, n))
    ax.set_aspect("equal")
    off += n * n
    n = n // 2 + 1
fig.savefig("mg_test.png", bbox_inches="tight", dpi=300)
plt.show()

