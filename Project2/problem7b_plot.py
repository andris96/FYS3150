import numpy as np
import matplotlib.pyplot as plt

#this code is basically copy-pasted from 7a, because I'm lazy

n = 100
vectors = np.loadtxt("eigenvec100.txt",skiprows=1)

x = np.linspace(0,1,n)
vec1 = np.zeros(len(x))
vec2 = np.zeros(len(x))
vec3 = np.zeros(len(x))
vec1_analytical = np.zeros(len(x))
vec2_analytical = np.zeros(len(x))
vec3_analytical = np.zeros(len(x))


#Taking the N sized vectors in vec1, vec2, and vec3
#The boundaries are already set as 0.0 above
for i in range(1,n-1):
    vec1[i] = vectors[i][0]
    vec2[i] = vectors[i][1]
    vec3[i] = vectors[i][2]



#plotting the analytical solutions
for i in range(1,n-1):
    vec1_analytical[i] = np.sin(i*np.pi/n)
    vec2_analytical[i] = np.sin(2*i*np.pi/n)
    vec3_analytical[i] = np.sin(3*i*np.pi/n)

#normalizing the vectors
vec1_analytical = vec1_analytical/np.linalg.norm(vec1_analytical)
vec2_analytical = vec2_analytical/np.linalg.norm(vec2_analytical)
vec3_analytical = vec3_analytical/np.linalg.norm(vec3_analytical)

#plotting
plt.plot(x, vec1,"b", label = "First eigenvector")
plt.plot(x, vec2,"r", label = "Second eigenvector")
plt.plot(x, vec3,"g", label = "Third eigenvector")
plt.plot(x, vec1_analytical,"b-.", label = "First analytical eigenvector")
plt.plot(x, vec2_analytical,"r-.", label = "Second analytical eigenvector")
plt.plot(x, vec3_analytical,"g-.", label = "Third analytical eigenvector")
plt.ylabel("$v_i$")
plt.xlabel("$\hat{x}_i$")
plt.title("Eigenvectors as function of x_i with n = 100 steps")
plt.legend()
plt.savefig("eigenvec100.pdf")