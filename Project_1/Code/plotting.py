import numpy as np
import matplotlib.pyplot as plt

#Using numpy.loadtxt to copy the values from the file. Skipping the first row.
#this is a very simplistic program as i didn't have time to make higher quality

data = np.loadtxt("data.txt",skiprows=1)

x = np.zeros(len(data))
u = np.zeros(len(data))


#Saving the x-values in one array and u in another.
for i in range(len(data)):
    x[i] = data[i][0]
    u[i] = data[i][1]

#plotting the graph and saving the plot
plt.figure()
plt.plot(x,u,"-o")
plt.title("Analytical u(x)")
plt.xlabel("x")
plt.ylabel("u(x)")
plt.savefig("data_plott.pdf")
plt.show()

#same procedure for the other plots
data_v1 = np.loadtxt("v_x1.txt",skiprows=1)
data_v2 = np.loadtxt("v_x2.txt",skiprows=1)
data_v3 = np.loadtxt("v_x3.txt",skiprows=1)

v1 = np.zeros(len(data_v1))
v2 = np.zeros(len(data_v2))
v3= np.zeros(len(data_v3))
x1 = np.zeros(len(data_v1))
x2 = np.zeros(len(data_v2))
x3 = np.zeros(len(data_v3))

for i in range(len(data_v1)-2):
    x1[i] = data_v1[i][0]
    v1[i] = data_v1[i][1]

for i in range(len(data_v2)-2):
    x2[i] = data_v2[i][0]
    v2[i] = data_v2[i][1]

for i in range(len(data_v3)-2):
    x3[i] = data_v3[i][0]
    v3[i] = data_v3[i][1]


plt.plot(x1,v1,"o")
plt.title("Numerical approximation v(x) with n = 10")
plt.xlabel("x")
plt.ylabel("v(x)")
plt.savefig("v1_plot.pdf")
plt.show()

plt.plot(x2,v2,"o")
plt.title("Numerical approximation v(x) with n = 100")
plt.xlabel("x")
plt.ylabel("v(x)")
plt.savefig("v2_plot.pdf")
plt.show()

plt.plot(x3,v3,"o")
plt.title("Numerical approximation v(x) with n = 1000")
plt.xlabel("x")
plt.ylabel("v(x)")
plt.savefig("v3_plot.pdf")
plt.show()
