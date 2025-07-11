import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as stats
import matplotlib as mpl
mpl.rc_file("../../test/mplstyleerc")

# Data to plot
eng = np.array([100., 200., 500., 1000., 2000., 5000., 7000., 10000., 13000.])

"""
# ADD TO RUNGAPS
# Execution time vs Energy (Ensure kine-limit=1e-9, and maxPartons=100!)
if args.runtype == 'energy':
    # Remove previous results and make folder
    if os.path.exists('results'):
        shutil.rmtree('results')
    os.makedirs('results', exist_ok=True)

    # Clear previous log files
    if os.path.exists('gpu-time.dat'):
        os.remove('gpu-time.dat')

    # Run the comparison 100 times, for different number of events
    energyarray = [100., 200., 500., 1000.,
                   2000., 5000., 7000., 10000., 13000.]
    for e in energyarray:
        # Run and store the output in a log file
        for i in range(1, 11):
            print(f"Running GAPS with E = {e} GeV")
            subprocess.run(['./gpu/bin/gpu', str(args.nevents), str(e)])

    # Move the log files to the results directory
    shutil.move('gpu-time.dat', 'results/')
    shutil.move('gpu.yoda', 'results/')
"""

# Params
nreps = 10

# import data
full = np.genfromtxt("gpu-time-energies.dat", delimiter=',')

print(full)

# Calculate the average and standard deviation for all repetitions
median = np.zeros((len(eng), 4))
iqr = np.zeros((len(eng), 4))

for i in range(len(eng)):

    median[i] = np.median(full[i*nreps:i*nreps+nreps], axis=0)
    iqr[i] = stats.iqr(full[i*nreps:i*nreps+nreps], axis=0)

print(median)

# Convert the arrays to pandas DataFrames for easier printing
columns_lin = ['Matrix  Element', 'Parton  Shower', 'Observables', 'Total']
df = pd.DataFrame(median, index=eng, columns=columns_lin)

# Print the DataFrame
print("Execution time at different Energies:")
print(df)
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
    p1, c1 = np.polyfit(np.log(eng[3:]), np.log(median[3:, i]), 1, cov=True)

    # Linear Fit
    p[i//2, i % 2, :] = p1

    # Print the results
    print(labels[i], ":", round(p1[0], 3), "±", round(np.sqrt(c1[0, 0]), 3))


# Create a new figure with a 4x2 grid of subplots
fig, axs = plt.subplots(2, 2, figsize=(10, 6.4))

# Define the column names
columns = [['Matrix  Element', 'Parton  Shower'], ['Observables', 'Total']]

# Add this line to adjust the space between subplots
fig.subplots_adjust(wspace=1, hspace=1)

# Add linspace for the linear fit
x = np.linspace(900, 14000, 10000)

# Loop over the columns and plot the data
for i in range(2):
    for j in range(2):
        ax = axs[i, j]
        err = ax.errorbar(eng, median[:, 2*i + j],
                          yerr=iqr[:, 2*i + j], fmt='o', color='C0')
        lin = ax.plot(eng, median[:, 2*i + j], alpha=0.3, color='C0')
        fit = ax.plot(x, np.exp(p[i, j, 1]) * x**p[i, j, 0], color='C1',
                      linestyle='--')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Energy (GeV)')
        ax.set_ylabel('Execution time (s)')
        ax.set_title(columns[i][j])
        ax.grid(True)

        handles = [err, fit[0]]
        labels = ['Data', "Linear Fit, Gradient = " +
                  str(round(p[i, j, 0], 2))]

        ax.legend(handles, labels)


fig.tight_layout()
plt.savefig('executionTimevsEnergy.pdf')
