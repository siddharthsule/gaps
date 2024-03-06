import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import matplotlib as mpl
mpl.rc_file("mplstyleerc")

# Data to plot
nev = np.array([1, 2, 5, 10, 20, 50, 100, 200, 500, 1000,
                2000, 5000, 10000, 20000, 50000,
               100000, 200000, 500000, 1000000])
columns = ['Matrix  Element', 'Parton  Shower',
           'Observables', 'Total']  # Column names

# Params
nloops = 10

dir_to_results = "../results-original/"

# import data
cpp_full = np.genfromtxt(dir_to_results + "cpp-time.dat", delimiter=',')
cud_full = np.genfromtxt(dir_to_results + "cuda-time.dat", delimiter=',')
#csm_full = np.genfromtxt("../results-single/cuda-time.dat", delimiter=',')

cpp_avg = np.zeros((len(nev), 4))
cud_avg = np.zeros((len(nev), 4))
#csm_avg = np.zeros((len(nev), 4))

cpp_std = np.zeros((len(nev), 4))
cud_std = np.zeros((len(nev), 4))
#csm_std = np.zeros((len(nev), 4))

for i in range(len(nev)):
    cpp_avg[i] = np.mean(cpp_full[i*nloops:i*nloops+nloops], axis=0)
    cud_avg[i] = np.mean(cud_full[i*nloops:i*nloops+nloops], axis=0)
    #csm_avg[i] = np.mean(csm_full[i*nloops:i*nloops+nloops], axis=0)

    cpp_std[i] = np.std(cpp_full[i*nloops:i*nloops+nloops], axis=0)
    cud_std[i] = np.std(cud_full[i*nloops:i*nloops+nloops], axis=0)
    #csm_std[i] = np.std(csm_full[i*nloops:i*nloops+nloops], axis=0)

cpp = cpp_avg
cud = cud_avg
#csm = csm_avg

# Convert the arrays to pandas DataFrames for easier printing
cpp_df = pd.DataFrame(cpp, index=nev, columns=columns)
cud_df = pd.DataFrame(cud, index=nev, columns=columns)

# Calculate the ratios and convert to integer
cpp_cud_ratio = (cpp / cud)

# Convert the ratio arrays to DataFrame
cpp_cud_ratio_df = pd.DataFrame(cpp_cud_ratio, index=nev, columns=columns)

# Print the DataFrame
print("C++ / CUDA Ratio:")
print(cpp_cud_ratio_df)
print("\n")

# Create a new figure with a 4x2 grid of subplots
fig, axs = plt.subplots(4, 2, figsize=(10, 12.8))

# Add this line to adjust the space between subplots
fig.subplots_adjust(wspace=0.5, hspace=0.5)

for i in range(4):  # Loop over columns
    
    for j in range(2):  # Loop over rows

        if i == 0:

            lop = np.genfromtxt("../doc/extra/loPointKernel.csv", delimiter=',')
            lop = lop[:, 1] / 1e9

            # From Single Profiling Session
            axs[i, j].scatter(nev, lop, color='C3', label='LO  Point  Kernel')

        axs[i, j].errorbar(nev, cpp[:, i], yerr=cpp_std[:, i], fmt='o', color="C0", label='C++')
        axs[i, j].plot(nev, cpp[:, i], color='C0', alpha=0.3)

        axs[i, j].errorbar(nev, cud[:, i], yerr=cud_std[:, i], fmt='o', color="C2", label='CUDA')
        axs[i, j].plot(nev, cud[:, i], color='C2', alpha=0.3)

        #axs[i, j].scatter(nev, csm[:, i], color='C1', label='CUDA  Single  Emission')
        #axs[i, j].plot(nev, csm[:, i], color='C1', alpha=0.5)

        axs[i, j].set_xlabel('Number  of  events')
        axs[i, j].set_ylabel('Execution  time [s]')

        axs[i, 0].set_title(f'{columns[i]}, linear scale')
        axs[i, 1].set_title(f'{columns[i]}, log-log scale')
        axs[i, j].legend()
        axs[i, j].grid(True)

        if j == 1:

            axs[i, j].set_xscale('log')
            axs[i, j].set_yscale('log')

fig.tight_layout()
plt.savefig('executiontime_combined.pdf')

# ------------------------------------------------------------------------------

# Plotting Number of Completed Events per cycle and Number of newly Completed
# Events per cycle

# Plot
fig, ax = plt.subplots(1, 2, figsize=(9, 3.75))

nev2 = [1000, 10000, 100000, 1000000]

for n in nev2:

    filename = "../results-original/cycles/cuda-cycles-" + str(n) + "-" + str(1)

    temp = np.genfromtxt(filename + ".dat", delimiter='\n')
    temp /= n  # Divide by number of events

    comp = np.zeros((200))
    diff = np.zeros((200))
    max = 0

    comp[:len(temp)] = temp
    comp[len(temp):] = temp[-1]

    for i in range(1, len(temp) - 1):

        diff[i] = temp[i] - temp[i-1]

    if len(temp) > max:
        max = len(temp)

    comp = comp[:max]
    comp = np.append(comp, 0.0)

    diff = diff[:max]
    diff = np.append(diff, 0.0)

    cycles = np.arange(1, len(comp) + 1)

    ax[0].scatter(cycles, comp, label= str(n) + '  Events ')
    ax[0].plot(cycles, comp, alpha=0.5)
    ax[0].set_xlabel('Cycle')
    ax[0].set_ylabel('Number  of  Completed  Events / Total')
    ax[0].set_title('Number  of  Completed  Events  per  Cycle')
    ax[0].legend()
    ax[0].grid(True)

    ax[1].scatter(cycles, diff, label=str(n) + '  Events ')
    ax[1].plot(cycles, diff, alpha=0.5)
    ax[1].set_xlabel('Cycle')
    ax[1].set_ylabel('Number  of  Newly  Completed  Events / Total')
    ax[1].set_title('Number  of  Newly  Completed   Events  per  Cycle')
    ax[1].legend()
    ax[1].grid(True)

fig.tight_layout()
plt.savefig("completedEvents.pdf")
