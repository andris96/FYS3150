import matplotlib.pyplot as plt
import numpy as np

# Load data
data_burn_in = {
    "1" : {
        "R" : {
            "e" : np.loadtxt("T1R_e_values.txt"),
            "m" : np.loadtxt("T1R_m_values.txt")
        },
        "O" : {
            "e" : np.loadtxt("T1O_e_values.txt"),
            "m" : np.loadtxt("T1O_m_values.txt")
        },
    }, 
    "2" : {
        "R" : {
            "e" : np.loadtxt("T2R_e_values.txt"),
            "m" : np.loadtxt("T2R_m_values.txt")
        },
        "O" : {
            "e" : np.loadtxt("T2O_e_values.txt"),
            "m" : np.loadtxt("T2O_m_values.txt")
        },
    } 
}

cycles = np.loadtxt("cycles.txt")

# Plotting the evolution of <e> and <m> for an increasing number of 
# Monte Carlo cycles, ie. burn-in (Problem)
for T in ["1", "2"]:
        for value in ["e", "m"]:
            title = "T = 1.0 J/kB" if T=="1" else "T = 2.4 J/kB"
            ylabel = r"<$\epsilon$>" if value=="e" else "<|m|>"
            file_name = "plot_burn_in_T" + T + "_" + value +  ".pdf"

            plt.figure()
            plt.plot(cycles, data_burn_in[T]["O"][value], label = "Ordered")
            plt.plot(cycles, data_burn_in[T]["R"][value], label = "Random")
            plt.xlabel("Cycles")
            plt.ylabel(ylabel)
            plt.title(title)
            plt.legend()
            plt.savefig(file_name)




