import numpy as np
import matplotlib.pyplot as plt



z = np.loadtxt("motion_z_RK4.txt")
t = np.loadtxt("time_interval.txt")


plt.plot(t,z, label = "Motion in Z direction")
plt.xlabel("Time [$\mu s$]")
plt.ylabel("Position [$\mu m$]")
plt.legend()
plt.savefig("zplot.pdf")
