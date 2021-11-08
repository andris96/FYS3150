import re
import glob

import numpy as np
import matplotlib.pyplot as plt

#foldername = "outputs"
filenames_list = glob.glob("*.txt")

fig, ax = plt.subplots(1, 1)
for filename in filenames_list:
    if "fractions_f_" in filename:
        f = filename[-8:-4] # Extracting f value
        data = np.loadtxt(filename)
        omega_V, fraction = data[:, 0], data[:, 1]
        ax.plot(omega_V, fraction, label=f"f={f}")
        ax.legend()

plt.xlabel(r"$\omega_V$ [MHz]")
plt.ylabel("fraction")
fig.savefig("fractions_plots.pdf")
plt.close()