import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Load data
z = np.loadtxt("motion_z_RK4.txt")
z_a = np.loadtxt("motion_z_analytical.txt")
t = np.loadtxt("time_interval.txt")

r1_with_interactions = np.loadtxt("motion_r_1_with_interactions.txt")
v1_with_interactions = np.loadtxt("motion_v_1_with_interactions.txt")
r2_with_interactions = np.loadtxt("motion_r_2_with_interactions.txt")
v2_with_interactions = np.loadtxt("motion_v_2_with_interactions.txt")

r1_without_interactions = np.loadtxt("motion_r_1_without_interactions.txt")
v1_without_interactions = np.loadtxt("motion_v_1_without_interactions.txt")
r2_without_interactions = np.loadtxt("motion_r_2_without_interactions.txt")
v2_without_interactions = np.loadtxt("motion_v_2_without_interactions.txt")

# Structering data for easier access
data_2p_system = {
    "With" : {
        "1" : { "r" : r1_with_interactions, "v" : v1_with_interactions
        },
        "2" : { "r" : r2_with_interactions, "v" : v2_with_interactions
        }
    },
    "Without" : {
        "1" : { "r" : r1_without_interactions, "v" : v1_without_interactions
        },
        "2" : { "r" : r2_without_interactions, "v" : v2_without_interactions
        }
    }
}

# Plot the motion of a single particle in the z direction as a funstion of time
#plt.plot(t, z, label="Motion in Z direction")
plt.plot(t, z_a, label="Motion in Z direction, analytical")
plt.xlabel("Time [$\mu s$]")
plt.ylabel("Position, z [$\mu m$]")
plt.legend()
plt.savefig("z_plot.pdf")
plt.clf()

# Plot the motion of the two particles in the xy-plane
fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True, figsize=plt.figaspect(0.5))
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)

axes = [ax1, ax2]
for idx, mode in enumerate(["With", "Without"]):
    for particle in ["1", "2"]:    
        x = data_2p_system[mode][particle]["r"][:, 0]
        y = data_2p_system[mode][particle]["r"][:, 1]
        axes[idx].plot(x, y, label=f"Particle {particle}")
        axes[idx].set_title(f"{mode} particle interactions")
        axes[idx].legend()

fig.tight_layout(pad=2)
plt.xlabel("x [$\mu m$]")
plt.ylabel("y [$\mu m$]")
fig.savefig("xy_plots.pdf")
plt.close()


# Phase space plots for the two particles: position vs. velocity for each direction
# NOT COMPLETED!! Need minor fixes to make it "pretty"
fig, axes = plt.subplots(3, 2, sharex=True, sharey=True, figsize=(8,12))
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)

for row, dir_ in enumerate(["x", "y", "z"]):
    for col, mode in enumerate(["With", "Without"]):
        for particle in ["1", "2"]:
            axes[row, col].plot(data_2p_system[mode][particle]["r"][:,row], 
                                data_2p_system[mode][particle]["v"][:,row], 
                                label=f"Particle {particle}")
            axes[row, col].set_title(f"{mode} particle interactions")
            axes[row, col].legend()

fig.tight_layout(pad=2)
plt.xlabel("x [$\mu m$]")
plt.ylabel("y [$\mu m$]")
fig.savefig("phase_plots.pdf")
plt.close()

# 3D plot of trajectories
# NOT COMPLETE! Both particles not appearing
fig = plt.figure(figsize=plt.figaspect(0.5))
axes = [ax1, ax2]
for idx, mode in enumerate(["With", "Without"]):
    for particle in ["1", "2"]:
        x = data_2p_system[mode][particle]["r"][:, 0]
        y = data_2p_system[mode][particle]["r"][:, 1]
        z = data_2p_system[mode][particle]["r"][:, 2]

        axes[idx] = fig.add_subplot(1, 2, idx+1, projection='3d')

        axes[idx].plot(x, y, z, marker="x", label=f"Particle {particle}")
        axes[idx].set_title(f"{mode} particle interactions")
        axes[idx].legend()

plt.savefig("3D_plot.pdf")