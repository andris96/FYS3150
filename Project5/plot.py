import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

dt = 2.5E-5
T = 0.002
tsteps = int(T/dt)

data = np.loadtxt("p.txt")

t = np.linspace(0,T,tsteps)

plt.plot(t,data)
plt.savefig("p.pdf")



h = 0.005
x_points = np.arange(0, 1+h, h)
y_points = np.arange(0, 1+h, h)
x, y = np.meshgrid(x_points, y_points, sparse=True)


U_data = np.loadtxt("U_real.txt")

U_data = U_data.reshape(80,201,201)



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

norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(U_data[0]))

# Plot the first frame
img = ax.imshow(U_data[0], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm)

# Axis labels
plt.xlabel("y", fontsize=fontsize)
plt.ylabel("x", fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)

# Add a text element showing the time
time_txt = plt.text(0.95, 0.95, "t = {:.3e}".format(t_min), color="white", 
                    horizontalalignment="right", verticalalignment="top", fontsize=fontsize)

                    # Function that takes care of updating the z data and other things for each frame


def animation(i):
    # Normalize the colour scale to the current frame?
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(U_data[i]))
    img.set_norm(norm)

    # Update z data
    img.set_data(U_data[i])

    # Update the time label
    current_time = t_min + i * dt
    time_txt.set_text("t = {:.3e}".format(current_time))

    return img

# Use matplotlib.animation.FuncAnimation to put it all together
anim = FuncAnimation(fig, animation, interval=1, frames=np.arange(0, len(U_data), 1), repeat=False, blit=0)

anim.save('./animation.mp4', writer="ffmpeg", bitrate=-1, fps=30)
