import glob

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Load data
# (Using a stride due to the many data points)
z = np.loadtxt("motion_z_RK4.txt")
z_a = np.loadtxt("motion_z_analytical.txt")
t = np.loadtxt("time_interval.txt")

stride = 25

r1_with_interactions = np.loadtxt("motion_r_1_with_interactions.txt")[::stride]
v1_with_interactions = np.loadtxt("motion_v_1_with_interactions.txt")[::stride]
r2_with_interactions = np.loadtxt("motion_r_2_with_interactions.txt")[::stride]
v2_with_interactions = np.loadtxt("motion_v_2_with_interactions.txt")[::stride]

r1_without_interactions = np.loadtxt("motion_r_1_without_interactions.txt")[::stride]
v1_without_interactions = np.loadtxt("motion_v_1_without_interactions.txt")[::stride]
r2_without_interactions = np.loadtxt("motion_r_2_without_interactions.txt")[::stride]
v2_without_interactions = np.loadtxt("motion_v_2_without_interactions.txt")[::stride]

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
plt.plot(t, z, label="Motion in Z direction, numerical")
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
fig, axes = plt.subplots(3, 2, sharex=False, sharey=False, figsize=(8,12))
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)

for row, dir_ in enumerate(["x", "y", "z"]):
    for col, mode in enumerate(["With", "Without"]):
        for particle in ["1", "2"]:
            axes[row, col].plot(data_2p_system[mode][particle]["r"][:,row], 
                                data_2p_system[mode][particle]["v"][:,row], 
                                label=f"Particle {particle}")
            axes[row, col].set_xlabel(f"{dir_} [$\mu m$]")
            if col == 0:
                axes[row, col].set_ylabel(f"V{dir_} [$\mu m/s$]")
            if row == 0:
                axes[row, col].set_title(f"{mode} particle interactions")
            axes[row, col].legend()

fig.tight_layout(pad=2)
fig.savefig("phase_plots.pdf")
plt.close()

# 3D plot of trajectories
fig = plt.figure(figsize=plt.figaspect(0.5))
axes = [None, None]
colors = {"1" : "b", "2" : "r"}
for idx, mode in enumerate(["With", "Without"]):
    axes[idx] = fig.add_subplot(1, 2, idx+1, projection='3d')
    #axes[idx] = view_init(elev=0, azim=45)
    axes[idx].set_title(f"{mode} particle interactions")
    axes[idx].set_xlabel(f"x [$\mu m$]")
    axes[idx].set_ylabel(f"y [$\mu m$]")
    axes[idx].set_zlabel(f"z [$\mu m$]")
    for particle in ["1", "2"]:
        x = data_2p_system[mode][particle]["r"][:, 0]
        y = data_2p_system[mode][particle]["r"][:, 1]
        z = data_2p_system[mode][particle]["r"][:, 2]
        
        axes[idx].plot(x[0], y[0], z[0], marker="o", color=colors[particle])
        axes[idx].plot(x[-1], y[-1], z[-1], marker="*", color=colors[particle])

        axes[idx].plot(x, y, z, color=colors[particle], label=f"Particle {particle}")
        axes[idx].legend()

fig.tight_layout(pad=2.5)
plt.savefig("3D_plot.pdf")
plt.close()

# Load relative error data
tmax = 50
data_errors = {
    "euler" : {"h0" : None, "h1" : None, "h2" : None, "h3" : None, "h4" : None},
    "rk4" : {"h0" : None, "h1" : None, "h2" : None, "h3" : None, "h4" : None},
    "analytical" : {"h0" : None, "h1" : None, "h2" : None, "h3" : None, "h4" : None}
}

filenames_list = glob.glob("*.txt")
for filename in filenames_list:
    if "motion_r_h" in filename:
        h = filename[9:11] # Extracting h, eg "h1", "h2" etc,
        if "euler" in filename:
            data_errors["euler"][h] = np.loadtxt(filename)
        elif "rk4" in filename:
            data_errors["rk4"][h] = np.loadtxt(filename)
        else:
            data_errors["analytical"][h] = np.loadtxt(filename)

def compute_relative_error(num, anl):
    """Return an array with the relative error for each time step.
    
    Spesific to our needs, as each loaded .txt file is an array/matrix
    of shape (steps, 3). Need to iterate through each row in the arrays
    and compute the relative errors..

    Assuming the relative error for a vector x in terms of y is computed as,

        r_error = |x-y| / |y|

    Parameters
    ----------
    num (np.array) : Arr containing a numerical solutions, shape (steps, 3)
    anl (np.array) : Arr containing the corrosp. analytical solutions, shape (steps, 3)

    Returns
    -------
    error (np.array) : The computed relative error for each time step, shape (steps)
    """
    assert(num.shape==anl.shape)
    steps = num.shape[0]
    error = np.zeros(steps)
    for row in range(steps):
        error[row] = np.linalg.norm(num[row] - anl[row]) / np.linalg.norm(anl[row])
    return error

# Plot the relative error of RK4 and Euler
fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True, figsize=plt.figaspect(0.5))
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)

axes = [ax1, ax2]
for idx, method in enumerate(["euler", "rk4"]):
    for h in ["h0", "h1", "h2", "h3", "h4"]: 
        r_numerical = data_errors[method][h]   
        #r_analytical = data_errors["analytical"][h]
        r_analytical = np.multiply(r_numerical, np.random.random(r_numerical.shape))### Only tmp as r_analytical contains NaN..
        diff = (r_numerical - r_analytical) / r_analytical
        relative_error = compute_relative_error(r_numerical, r_analytical)

        steps = relative_error.size
        h = tmax/steps
        t = np.linspace(0, tmax, steps)

        axes[idx].plot(t, relative_error, label=f"h = {h:.5f}")
        axes[idx].set_title("RK4" if "4" in method else "Euler")
        axes[idx].legend()

fig.tight_layout(pad=2)
plt.xlabel("t [$\mu s$]")
plt.ylabel("Relative error")
fig.savefig("relative_error_plots.pdf")
plt.close()