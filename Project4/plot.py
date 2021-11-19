#can change name of this file if needed


import matplotlib.pyplot as plt
import numpy as np

T1Oe = np.loadtxt("T1O_e_values.txt")
T1Om = np.loadtxt("T1O_m_values.txt")
T1Re = np.loadtxt("T1R_e_values.txt")
T1Rm = np.loadtxt("T1R_m_values.txt")
T2Oe = np.loadtxt("T2O_e_values.txt")
T2Om = np.loadtxt("T2O_m_values.txt")
T2Re = np.loadtxt("T2R_e_values.txt")
T2Rm = np.loadtxt("T2R_m_values.txt")

cycles = np.loadtxt("cycles.txt")

# Just temporary plotting, should make them look nice eventually

plt.plot(cycles, T1Oe, label = "Ordered")
plt.plot(cycles, T1Re, label = "Random")
plt.legend()
plt.show()

plt.plot(cycles, T1Om, label = "ordered")
plt.plot(cycles, T1Rm, label = "random")
plt.legend()
plt.show()

plt.plot(cycles, T1Oe, label = "ordered")
plt.plot(cycles, T1Re, label = "random")
plt.legend()
plt.show()

plt.plot(cycles, T2Om, label = "ordered")
plt.plot(cycles, T2Rm, label = "random")
plt.legend()
plt.show()


