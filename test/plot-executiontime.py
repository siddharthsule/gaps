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
print("C++ / CUDA Ratio:")
print(cpp_cud_ratio_df)
print("\n")

# Create a new figure with a 4x2 grid of subplots
fig, axs = plt.subplots(2, 2, figsize=(10, 6.4))

# Define the column names
columns = [['Matrix  Element', 'Parton  Shower'], ['Observables', 'Total']]

# Add this line to adjust the space between subplots
fig.subplots_adjust(wspace=1, hspace=1)

# Loop over the columns and plot the data
for i in range(2):
    for j in range(2):
        ax = axs[i, j]
        ax.errorbar(nev, cpp[:, 2*i + j], yerr=cpp_iqr[:, 2*i + j], fmt='o', label='C++', color='C0')
        ax.errorbar(nev, cud[:, 2*i + j], yerr=cud_iqr[:, 2*i + j], fmt='o', label='CUDA', color='C2')
        ax.plot(nev, cpp[:, 2*i + j], color='C0', alpha=0.3)
        ax.plot(nev, cud[:, 2*i + j], color='C2', alpha=0.3)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Number of events')
        ax.set_ylabel('Execution time (s)')
        ax.set_title(columns[i][j])
        ax.legend()
        ax.grid(True)

fig.tight_layout()
plt.savefig('executionTime.pdf')