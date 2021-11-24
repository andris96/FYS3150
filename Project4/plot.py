import matplotlib.pyplot as plt
import numpy as np

############################################
# Plotting the evolution of <e> and <m> for 
# an increasing number of Monte Carlo cycles
# Problem 5 : burn-in time
############################################

cycles = np.loadtxt("cycles.txt")
for T in ["1", "2"]:
    for value in ["e", "m"]:

        random = np.loadtxt("T" + T + "R" + "_" + value + "_values.txt")
        ordered = np.loadtxt("T" + T + "O" + "_" + value + "_values.txt")

        title = "T = 1.0 J/kB" if T=="1" else "T = 2.4 J/kB"
        ylabel = r"<$\epsilon$>" if value=="e" else "<|m|>"
        filename = "plot_burn_in_T" + T + "_" + value +  ".pdf"

        plt.figure()
        plt.plot(cycles, ordered, label = "Ordered")
        plt.plot(cycles, random, label = "Random")
        plt.xlabel("Cycles")
        plt.ylabel(ylabel)
        plt.title(title)
        plt.legend()
        plt.savefig(filename)


############################################
# Plotting a histogram containing epsilon
# samples with T = 1 and T = 2.4
# Problem 6 : probability distribution
############################################

samplesT1 = np.loadtxt("samplesT1.txt")
samplesT2 = np.loadtxt("samplesT2.txt")
width = 0.001

plt.figure()
plt.hist(samplesT1, bins = np.arange(min(samplesT1), max(samplesT1) + width, width))
plt.xlabel("$\epsilon$")
plt.ylabel("distribution")
plt.savefig("samplesT1.pdf")

plt.figure()
plt.hist(samplesT2, bins = np.arange(min(samplesT2), max(samplesT1) + width, width))
plt.xlabel("$\epsilon$")
plt.ylabel("distribution")
plt.savefig("samplesT2.pdf")


############################################
# Plotting expectation values as function of
# temperature, for different lattize sizes
# Problem 8 : phase transitions
############################################
quantities = ["e", "m", "Cv", "X"]
figures = {qnt : plt.figure() for qnt in quantities}
axes = {qnt : figures[qnt].add_subplot(1,1,1) for qnt in quantities}

temperatures = np.linspace(2.1, 2.4, 50)
for L in ["40", "60", "80", "100"]:
    data = np.loadtxt("expectation_valuesL" + L + ".txt")
    for col, qnt in enumerate(quantities):
        axes[qnt].plot(temperatures, data[:, col], label=f"L={L}")

for qnt in quantities:
    ylabel = r"<$\epsilon$>" if qnt=="e" else "<" + qnt + ">"
    axes[qnt].set_xlabel("T")
    axes[qnt].set_ylabel(ylabel)
    axes[qnt].legend()
    figures[qnt].savefig("plot_" + qnt + ".pdf")