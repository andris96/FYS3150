import re
import glob

import numpy as np
import matplotlib.pyplot as plt

#foldername = "outputs"
filenames_list = glob.glob("*.txt")

fig, ax = plt.subplots(1, 1)
for filename in filenames_list:
    if "fractions_fine_scan" in filename:
        mode = "without" if "without" in filename else "with" 
        data = np.loadtxt(filename)
        omega_V, fraction = data[:, 0], data[:, 1]
        ax.plot(omega_V, fraction, label=f"W{mode[1:]} particle interaction")
        ax.legend()

plt.xlabel(r"$\omega_V$ [MHz]")
plt.ylabel("fraction")
fig.savefig("fractions_fine_scan_plots.pdf")
plt.close()