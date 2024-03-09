import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import matplotlib as mpl
mpl.rc_file("mplstyleerc")

# Data to plot
nev = np.array([1, 2, 5, 10, 20, 50, 100, 200, 500, 1000,
                2000, 5000, 10000, 20000, 50000,
               100000, 200000, 500000, 1000000])

# Params
nloops = 10

dir_to_results = "../results-times/"

# import data
cpp_full = np.genfromtxt(dir_to_results + "cpp-time.dat", delimiter=',')
cud_full = np.genfromtxt(dir_to_results + "gaps-time.dat", delimiter=',')

cpp_avg = np.zeros((len(nev), 4))
cud_avg = np.zeros((len(nev), 4))

cpp_std = np.zeros((len(nev), 4))
cud_std = np.zeros((len(nev), 4))

for i in range(len(nev)):
    cpp_avg[i] = np.mean(cpp_full[i*nloops:i*nloops+nloops], axis=0)
    cud_avg[i] = np.mean(cud_full[i*nloops:i*nloops+nloops], axis=0)

    cpp_std[i] = np.std(cpp_full[i*nloops:i*nloops+nloops], axis=0)
    cud_std[i] = np.std(cud_full[i*nloops:i*nloops+nloops], axis=0)

cpp = cpp_avg
cud = cud_avg

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
        ax.errorbar(nev, cpp[:, 2*i + j], yerr=cpp_std[:, 2*i + j], fmt='o', label='C++', color='C0')
        ax.errorbar(nev, cud[:, 2*i + j], yerr=cud_std[:, 2*i + j], fmt='o', label='CUDA', color='C2')
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