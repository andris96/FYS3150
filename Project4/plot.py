import matplotlib.pyplot as plt
import numpy as np
import glob

############################################
# Plotting the evolution of <e> and <m> for 
# an increasing number of Monte Carlo cycles
# Problem 5 : burn-in time

# Note: The output files from burn_in.cpp
# differ for with and without 
# parallelization (PL). With PL gives one 
#.txt file for each thread, which has to be
# with the other files and restructured here 
# in python...
############################################
PL = True # need to choose this

if PL:
    evalues_from_all_threads_list = []
    filenames_list = glob.glob("*.txt")
    for filename in filenames_list:
        if "thread" in filename:
            tmp_arr = np.loadtxt(filename)
            if len(tmp_arr.shape) > 1:
                for row in tmp_arr:
                    evalues_from_all_threads_list.append(row)
            else:
                evalues_from_all_threads_list.append(tmp_arr)
    evalues = np.array(evalues_from_all_threads_list)

    evalues = evalues[np.argsort(evalues[:, 0])] # sort by values in first col
    cycles = evalues[:, 0]

    # row-format in evalues: {max_cycles T1Oe T1Om T1Re T1Rm T2Oe T2Om T2Re T2Rm}
    evalues_dict = {
        "1" : {
            "e" : {
                "O" : evalues[:, 1],
                "R" : evalues[:, 3],
            },
            "m" : {
                "O" : evalues[:, 2],
                "R" : evalues[:, 4],
            }
        },
        "2" : {
            "e" : {
                "O" : evalues[:, 5],
                "R" : evalues[:, 7],
            },
            "m" : {
                "O" : evalues[:, 6],
                "R" : evalues[:, 8],
            }
        }
    }

    for T in ["1", "2"]:
        for value in ["e", "m"]:

            ordered = evalues_dict[T][value]["O"]
            random = evalues_dict[T][value]["R"]
            
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

else:
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
widthT1 = 0.005
widthT2 = 0.001

plt.figure()
plt.hist(samplesT1, bins = np.arange(min(samplesT1), max(samplesT1) + widthT1, widthT1))
plt.xlabel("$\epsilon$")
plt.ylabel("distribution")
plt.savefig("samplesT1.pdf")

plt.figure()
plt.hist(samplesT2, bins = np.arange(min(samplesT2), max(samplesT2) + widthT2, widthT2))
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

# must set thse values correct according to last run in main.cpp!
n_T = 10 
temperatures = np.linspace(2.3, 2.4, n_T)
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