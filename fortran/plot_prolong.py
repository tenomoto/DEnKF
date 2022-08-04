import numpy as np
import matplotlib.pyplot as plt

q = np.fromfile("prolong.dat").reshape(9, 9)

plt.rcParams["font.size"] = 18
fig, ax = plt.subplots(figsize=[10, 5])
ax.matshow(q)
plt.show()

