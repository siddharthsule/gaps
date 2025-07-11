import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import iqr

import matplotlib as mpl
mpl.rc_file("mplstyleerc")

# Data to plot
nev = np.array([1, 2, 5, 10, 20, 50, 100, 200, 500, 1000,
                2000, 5000, 10000, 20000, 50000,
               100000, 200000, 500000, 1000000])

# Params
nreps = 100

dir_to_results = "../results-times/"

# import data
cpu_full = np.genfromtxt(dir_to_results + "cpu-time.dat", delimiter=',')
gpu_full = np.genfromtxt(dir_to_results + "gaps-time.dat", delimiter=',')

# Calculate the average and standard deviation for all repetitions
cpu_median = np.zeros((len(nev), 4))
gpu_median = np.zeros((len(nev), 4))

cpu_iqr = np.zeros((len(nev), 4))
gpu_iqr = np.zeros((len(nev), 4))

for i in range(len(nev)):

    """
    # To Check if there are any outliers
    if nev[i] == 5000:
        print(gpu_full[i*nreps:i*nreps+nreps])

        fig, ax = plt.subplots()
        ax.hist(gpu_full[i*nreps:i*nreps+nreps, 0], bins=50)
        fig.savefig("histo.pdf")
    """

    cpu_median[i] = np.median(cpu_full[i*nreps:i*nreps+nreps], axis=0)
    gpu_median[i] = np.median(gpu_full[i*nreps:i*nreps+nreps], axis=0)

    cpu_iqr[i] = iqr(cpu_full[i*nreps:i*nreps+nreps], axis=0)
    gpu_iqr[i] = iqr(gpu_full[i*nreps:i*nreps+nreps], axis=0)

cpu = cpu_median
gpu = gpu_median

# Convert the arrays to pandas DataFrames for easier printing
columns_lin = ['Matrix  Element', 'Parton  Shower', 'Observables', 'Total']
cpu_df = pd.DataFrame(cpu, index=nev, columns=columns_lin)
gpu_df = pd.DataFrame(gpu, index=nev, columns=columns_lin)

# Calculate the ratios and convert to integer
cpu_gpu_ratio = (cpu / gpu)

# Convert the ratio arrays to DataFrame
cpu_gpu_ratio_df = pd.DataFrame(cpu_gpu_ratio, index=nev, columns=columns_lin)

# Print the DataFrame
print("CPU / GPU Ratio for different Number of Events:")
print(cpu_gpu_ratio_df)
print("\n")

# Initialize p as a 3D array
p = np.zeros((2, 2, 2))

# Define the labels for printing
labels = ["GPU Matrix Element Gradient", "GPU Parton Shower Gradient",
          "GPU Observables Gradient", "GPU Total Gradient"]

# Loop over the columns - to prevent lots of repeated statements
"""
For i = 0, i // 2 is 0 and i % 2 is 0.
For i = 1, i // 2 is 0 and i % 2 is 1.
For i = 2, i // 2 is 1 and i % 2 is 0.
For i = 3, i // 2 is 1 and i % 2 is 1.

Kept it here because I thought it was a neat way of doing loops!
"""
for i in range(4):

    # Linea Fit
    p1, c1 = np.polyfit(np.log(nev[14:]), np.log(gpu[14:, i]), 1, cov=True)

    # Linear Fit
    p[i//2, i % 2, :] = p1

    # Print the results
    print(labels[i], ":", round(p1[0], 3), "±", round(np.sqrt(c1[0, 0]), 3))

# print(p)

# Create a new figure with a 4x2 grid of subplots
fig, axs = plt.subplots(2, 2, figsize=(10, 6.4))

# Define the column names
columns = [['Matrix  Element', 'Parton  Shower'], ['Observables', 'Total']]

# Add this line to adjust the space between subplots
fig.subplots_adjust(wspace=1, hspace=1)

# Add linspace for the linear fit
x = np.linspace(40000, 1300000, 1000)

# Loop over the columns and plot the data
for i in range(2):
    for j in range(2):
        ax = axs[i, j]
        cpu_errorbar = ax.errorbar(
            nev, cpu[:, 2*i + j], yerr=cpu_iqr[:, 2*i + j], fmt='o', color='C0')
        gpu_errorbar = ax.errorbar(
            nev, gpu[:, 2*i + j], yerr=gpu_iqr[:, 2*i + j], fmt='o', color='C2')
        cpu_plot = ax.plot(nev, cpu[:, 2*i + j], color='C0', alpha=0.3)
        gpu_plot = ax.plot(nev, gpu[:, 2*i + j], color='C2', alpha=0.3)
        fit_plot = ax.plot(x, np.exp(p[i, j, 1]) * x**p[i, j, 0], color='C1',
                           linestyle='--')

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Number of events')
        ax.set_ylabel('Execution time (s)')
        ax.set_title(columns[i][j])
        ax.grid(True)

        # Add a vertical line
        ax.axvline(x=5120, color='C2', linestyle='--')

        # Create a proxy artist for the axvline to use in the legend
        v100_gpu_cores_line = mpl.lines.Line2D(
            [], [], color='C2', label="V100 GPU Cores")

        # Create a list of handles and labels manually, including the proxy artist
        handles = [cpu_errorbar, gpu_errorbar,
                   v100_gpu_cores_line, fit_plot[0]]
        labels = ['CPU', 'GPU',  "V100 GPU Cores",
                  "Linear Fit, Gradient = " + str(round(p[i, j, 0], 2))]

        # Add the legend with the updated handles and labels
        ax.legend(handles, labels)

        ax.legend(handles, labels)

fig.tight_layout()
plt.savefig('executionTime.pdf')
