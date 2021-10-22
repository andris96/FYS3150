import numpy as np
import matplotlib.pyplot as plt

# Load data
z = np.loadtxt("motion_z_RK4.txt")
t = np.loadtxt("time_interval.txt")

r1_with_interactions = np.loadtxt("motion_r_1_with_interactions.txt")
v1_with_interactions = np.loadtxt("motion_v_1_with_interactions.txt")
r2_with_interactions = np.loadtxt("motion_r_2_with_interactions.txt")
v2_with_interactions = np.loadtxt("motion_v_2_with_interactions.txt")

r1_without_interactions = np.loadtxt("motion_r_1_without_interactions.txt")
v1_without_interactions = np.loadtxt("motion_v_1_without_interactions.txt")
r2_without_interactions = np.loadtxt("motion_r_2_without_interactions.txt")
v2_without_interactions = np.loadtxt("motion_v_2_without_interactions.txt")

# Plot the motion of a single particle in the z direction as a funstion of time
plt.plot(t, z, label="Motion in Z direction")
plt.xlabel("Time [$\mu s$]")
plt.ylabel("Position, z [$\mu m$]")
plt.legend()
plt.savefig("zplot.pdf")

# Plot the motion of the two particles in the xy-plane
plt.plot(r1_with_interactions[:,0], r1_with_interactions[:,1], label="Particle 1")
plt.plot(r2_with_interactions[:,0], r2_with_interactions[:,1], label="Particle 2")
plt.plot(r1_without_interactions[:,0], r1_without_interactions[:,1], label="Particle 1")
plt.plot(r2_without_interactions[:,0], r2_without_interactions[:,1], label="Particle 2")
plt.xlabel("x [$\mu m$]")
plt.ylabel("y [$\mu m$]")
plt.legend()
plt.savefig("xyplot.pdf")

# Phase space plots..
## NOT COMPLETED...

# plt.plot(r1_with_interactions[:,0], v1_with_interactions[:,0] label="Particle 1")
# plt.plot(v2_with_interactions[:,0], v2_with_interactions[:,0] label="Particle 2")
# plt.xlabel("x [$\mu m$]")
# plt.ylabel("$v_x$ [$\mu m$]")
# plt.legend()
# plt.savefig("xv_xplot.pdf")

# plt.plot(r1_with_interactions[:,0], r1_with_interactions[:,1] label="Particle 1")
# plt.plot(r2_with_interactions[:,0], r2_with_interactions[:,1] label="Particle 2")
# plt.xlabel("x [$\mu m$]")
# plt.ylabel("y [$\mu m$]")
# plt.legend()
# plt.savefig("xyplot.pdf")

# plt.plot(r1_with_interactions[:,0], r1_with_interactions[:,1] label="Particle 1")
# plt.plot(r2_with_interactions[:,0], r2_with_interactions[:,1] label="Particle 2")
# plt.xlabel("x [$\mu m$]")
# plt.ylabel("y [$\mu m$]")
# plt.legend()
# plt.savefig("xyplot.pdf")
