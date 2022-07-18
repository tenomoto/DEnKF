import numpy as np
import matplotlib.pyplot as plt
from la import nstep, nens, nobs, dtobs, r

timestep = np.arange(0, nstep, dtobs)

l2 = np.load("l2.npy")
sprd = np.load("sprd.npy")

plt.rcParams["font.size"] = 18
fig, ax = plt.subplots(figsize=[10, 5])
ax.plot(timestep, l2, label="RMSE")
ax.plot(timestep, sprd, label="spread")
ax.axhline(np.sqrt(r), color="gray", linestyle=":")
ax.legend()
ax.set_xlabel("time step")
ax.set_ylabel("RMSE/spread")
ax.set_title(f"DEnKF M{nens} {nobs} obs every {dtobs} steps")
fig.savefig(f"la_denkf_M{nens}_{nobs}obs_dt{dtobs}.png", bbox_inches="tight", dpi=300)
#plt.show()

