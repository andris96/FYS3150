import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def heatmap2D(z, xmin, xmax, ymin, ymax, filename, xlabel, ylabel, cbar_label, cmap):
    """Produce heatmap plot of $z and save as $filename .pdf
    
    Partly adapted from:
    https://matplotlib.org/stable/gallery/images_contours_and_fields/image_annotated_heatmap.html
    
    params:
    -------
    z (np.array) : 2D array of shape (N, N)
    xmin, xmax, ymin, ymax (float) : Setting boundaries for X- and Y-axis
    filename (string) : Name of the file the figure should be saved to (without extension)
    xlabel (string) : Label of x-axis
    ylabel (string) : Label of y-axis
    cbar_label (string) : Label of the color bar or z-axis
    cmap (string) : Matplotlib color map
    """
    ax = plt.gca()

    im = ax.imshow(z, cmap=cmap, vmin=z.min(), vmax=z.max(), extent=(xmin, xmax, ymin, ymax))
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel(cbar_label, rotation=-90, va="bottom")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.figure.savefig(f"{filename}.pdf")
    plt.close()

# Initializing params
dt = 2.5E-5
T = 0.002
M = 201
tsteps = int(T/dt)
h = 0.005
x_points = np.arange(0, 1+h, h)
y_points = np.arange(0, 1+h, h)
x, y = np.meshgrid(x_points, y_points, sparse=True)

#####################################################
# Problem 7

# Plots of the deviation from of the total probablity
# from 1.0 as a functions of time. 

# Case 1 : no slit barrier
# Case 2 : double-slit barrier
#####################################################
p_sum_cases = [np.fromfile("p_sum.bin")]
#p_sum_cases = [np.fromfile("p_sum_case1.bin"), np.fromfile("p_sum_case2.bin")]
t = np.linspace(0, T ,tsteps)
for case, p_sum in enumerate(p_sum_cases):
    plt.plot(t, p_sum)
    plt.xlabel("Time")
    plt.ylabel("Probability")
    plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
    plt.title("") # Title or no title?
    plt.savefig(f"p_sum_{case+1}.pdf")


#####################################################
# Problem 8

# Heatmap plots that illustrate the time evolution of
# the 2D probability function p^n, as well as for 
# Re(U) and Im(U), for times t=0, t=0.001 and t=0.002
#####################################################
p = np.fromfile("p.bin").reshape(tsteps, M, M)
U_real = np.fromfile("U_real.bin").reshape(tsteps, M, M)
U_img = np.fromfile("U_imag.bin").reshape(tsteps, M, M)

data = {"p" : p, "U_real" : U_real, "U_img" : U_img}

for t_idx in [0, tsteps//2, -1]:
    for quantity in data:
        heatmap2D(z=data[quantity][t_idx, :, :],
                  xmin=x.min(), xmax=x.max(), 
                  ymin=y.min(), ymax=y.max(),
                  filename=f"heatmap_{quantity}_{t_idx}", 
                  xlabel="x", ylabel="y",
                  cbar_label=f"{quantity}", cmap="viridis")

#####################################################
# Problem 9

# Plot of the (normalized) probability distribution
# at x=0.8 for time t=0.002
#####################################################
p_single = None #np.fromfile("p_single.bin").reshape(tsteps, M, M)
p_triple = None #np.fromfile("p_trippel.bin").reshape(tsteps, M, M)
ps = {"single slit" : p_single, "double slit" : p, "triple slit" : p_triple}

for slit in ps:
    if ps[slit] is not None: # remove this if-condition when p_single.bin and p_triple.bin are valid files
        x_idx = int(0.8 / ((1 - 0) / M))
        p_given_x = ps[slit][-1, :, x_idx] # p(y|x=0.8;t=0.002), t=0.002 eqv. to last time idx
        p_given_x_normalized = p_given_x / np.linalg.norm(p_given_x, ord=1)

        plt.plot(y, p_given_x_normalized)
        plt.xlabel("y")
        plt.ylabel("P", rotation=90)
        plt.title(f"p(y|x=0.8;t=0.002), {slit}")
        plt.savefig(f"1D_prob_dist_{slit.split()[0]}_{slit.split()[1]}.pdf")

#####################################################
# Problem X

# Make animations of the double slit simulation
#####################################################
U_data = U_real 

# Some settings
fontsize = 12
t_min = t[0]
x_min, x_max = x_points[0], x_points[-1]
y_min, y_max = y_points[0], y_points[-1]

# Create figure
fig = plt.figure()
ax = plt.gca()

#plt.imshow(U_data[-1])
#plt.savefig("U.pdf")

# Create a colour scale normalization according to the max U_data value in the first frame
norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(U_data[0]))

# Plot the first frame
img = ax.imshow(U_data[0], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm)

# Axis labels
plt.xlabel("x", fontsize=fontsize)
plt.ylabel("y", fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)

# Add a text element showing the time
time_txt = plt.text(0.95, 0.95, "t = {:.3e}".format(t_min), color="white", 
                    horizontalalignment="right", verticalalignment="top", fontsize=fontsize)

def animation(i):
    """Takes care of updating the U data and other things for each frame"""
    # Normalize the colour scale to the current frame?
    #U_max = np.max(U_data[i])
    #U_max_root = np.sqrt(U_max)
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(U_data[i]))
    img.set_norm(norm)

    # Update U data
    img.set_data(U_data[i])

    # Update the time label
    current_time = t_min + i * dt
    time_txt.set_text("t = {:.3e}".format(current_time))

    return img

# Use matplotlib.animation.FuncAnimation to put it all together
anim = FuncAnimation(fig, animation, interval=1, frames=np.arange(0, len(U_data), 2), repeat=False, blit=0)
anim.save('./animation.mp4', writer="ffmpeg", bitrate=-1, fps=30)
