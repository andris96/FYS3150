import numpy as np
import matplotlib.pyplot as plt

#plotting the results from problem 6, I could make it more dynamic by making
#the c++ program give an output file and read it here, but I decide not to spend much time
#on that as project 1 already involved doing that
x = np.array([10,15,20,25,30,35,40,45,50,55,60,65,70])
y = np.array([148,351,649,1028,1477,2056,2698,3391,4174,5117,6092,7126,8341])

#making a plot of f(x) = x^2 to compare
x1 = np.linspace(10,70,100)
def pow2(x):
    return x**2

#plotting
plt.plot(x,y, label = "iterations")
plt.plot(x1,pow2(x1), label = "x^2")
plt.title("Number of iterations as a function of matrix size")
plt.xlabel("N")
plt.ylabel("Iterations")
plt.legend()
plt.savefig("iterations_plot.pdf")

