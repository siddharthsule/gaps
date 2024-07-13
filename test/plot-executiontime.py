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
cpp_full = np.genfromtxt(dir_to_results + "cpp-time.dat", delimiter=',')
cud_full = np.genfromtxt(dir_to_results + "gaps-time.dat", delimiter=',')

# Calculate the average and standard deviation for all repetitions
cpp_median = np.zeros((len(nev), 4))
cud_median = np.zeros((len(nev), 4))

cpp_iqr = np.zeros((len(nev), 4))
cud_iqr = np.zeros((len(nev), 4))

for i in range(len(nev)):

    """
    # To Check if there are any outliers
    if nev[i] == 5000:
        print(cud_full[i*nreps:i*nreps+nreps])

        fig, ax = plt.subplots()
        ax.hist(cud_full[i*nreps:i*nreps+nreps, 0], bins=50)
        fig.savefig("histo.pdf")
    """

    cpp_median[i] = np.median(cpp_full[i*nreps:i*nreps+nreps], axis=0)
    cud_median[i] = np.median(cud_full[i*nreps:i*nreps+nreps], axis=0)

    cpp_iqr[i] = iqr(cpp_full[i*nreps:i*nreps+nreps], axis=0)
    cud_iqr[i] = iqr(cud_full[i*nreps:i*nreps+nreps], axis=0)

cpp = cpp_median
cud = cud_median

# Convert the arrays to pandas DataFrames for easier printing
columns_lin = ['Matrix  Element', 'Parton  Shower', 'Observables', 'Total']
cpp_df = pd.DataFrame(cpp, index=nev, columns=columns_lin)
cud_df = pd.DataFrame(cud, index=nev, columns=columns_lin)

# Calculate the ratios and convert to integer
cpp_cud_ratio = (cpp / cud)

# Convert the ratio arrays to DataFrame
cpp_cud_ratio_df = pd.DataFrame(cpp_cud_ratio, index=nev, columns=columns_lin)

# Print the DataFrame
print("CPU / GPU Ratio for different Number of Events:")
print(cpp_cud_ratio_df)
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
    p1, c1 = np.polyfit(np.log(nev[14:]), np.log(cud[14:, i]), 1, cov=True)

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
        cpp_errorbar = ax.errorbar(
            nev, cpp[:, 2*i + j], yerr=cpp_iqr[:, 2*i + j], fmt='o', color='C0')
        cud_errorbar = ax.errorbar(
            nev, cud[:, 2*i + j], yerr=cud_iqr[:, 2*i + j], fmt='o', color='C2')
        cpp_plot = ax.plot(nev, cpp[:, 2*i + j], color='C0', alpha=0.3)
        cud_plot = ax.plot(nev, cud[:, 2*i + j], color='C2', alpha=0.3)
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
        handles = [cpp_errorbar, cud_errorbar,
                   v100_gpu_cores_line, fit_plot[0]]
        labels = ['CPU', 'GPU',  "V100 GPU Cores",
                  "Linear Fit, Gradient = " + str(round(p[i, j, 0], 2))]

        # Add the legend with the updated handles and labels
        ax.legend(handles, labels)

        ax.legend(handles, labels)

fig.tight_layout()
plt.savefig('executionTime.pdf')
