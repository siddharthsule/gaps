import os
import sys
import shutil
import argparse
import subprocess
import numpy as np
import pandas as pd
from scipy.stats import iqr
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc_file("mplstyleerc")

# To do CPU Cluster or GPU or Plot Results?
parser = argparse.ArgumentParser(description='Run Full GAPS Analysis')
parser.add_argument('option', choices=[
                    'cpu', 'gpu', 'plot'], help='Select the option to run')
parser.add_argument('--ncpu', type=int, default=64,
                    help='Number of CPU cores to use (default: 64)')
parser.add_argument('--nreps', type=int, default=10,
                    help='Number of repetitions for each run (default: 10)')
args = parser.parse_args()
option = args.option
ncpu = args.ncpu
nreps = args.nreps

# Different Values of Events
n_list = [10000, 100000, 1000000]

# GPU Tuning - best value for Threads per Block
thr = 128

# Emission filenames
em_gpu = 'gpu-emissions.dat'
em_cpu = 'cpu-cluster-emissions.dat'

# ------------------------------------------------------------------------------
# If CPU or GPU

if option == 'cpu' or option == 'gpu':

    # Go to base directory
    os.chdir("..")

    # Clear old Emission Data
    emissions_file = em_gpu if option == 'gpu' else em_cpu
    if os.path.exists(emissions_file):
        os.remove(emissions_file)

    # Run the simulation for each number of events
    for n in n_list:

        # Common for both
        command = f"./rungaps -p LHC -nlo -codecarbon -n {n}"

        # GPU Specific
        if option == 'gpu':
            command += f" -t {thr}"
        elif option == 'cpu':
            command += f" -r cpu-cluster -ncpu {ncpu}"

        # Run the Simulation
        for i in range(nreps):
            subprocess.run(command, shell=True)

    # Replace all the CodeCarbon Data into data we need to plot
    try:
        data = pd.read_csv(emissions_file, sep=",", header=0)
        data = data[data["duration"] > 0.5]
        data = data[["duration", "cpu_power", "gpu_power", "energy_consumed"]]
        data["power"] = data["cpu_power"] + data["gpu_power"]
        data = data[["duration", "energy_consumed"]]
        data.to_csv(emissions_file, sep=",", index=False)
    except FileNotFoundError:
        print(f"Error: {emissions_file} not found. Ensure codecarbon.")

    # Return to test directory
    os.chdir("test")

    # Exit after CPU/GPU processing is complete
    print(f"Completed {option} analysis. Results saved to {emissions_file}")
    exit(0)


# ------------------------------------------------------------------------------
# Plotting Results

if option == 'plot':
    # First, ensure CPU and GPU files exist
    if not os.path.exists("../" + em_cpu):
        print(f"Error: {em_cpu} not found.")
        exit(1)
    if not os.path.exists("../" + em_gpu):
        print(f"Error: {em_gpu} not found.")
        exit(1)

    gpu = pd.read_csv("../" + em_gpu, header=0, delimiter=',')
    cpu = pd.read_csv("../" + em_cpu, header=0, delimiter=',')

    print("GPU Emissions Data:")
    print(gpu)

    print("CPU Emissions Data:")
    print(cpu)

    def median_and_iqr(data):
        """
        Calculate the median and interquartile range (IQR) of the data.
        """

        # the dataset is 10 sets of 10 repetitions. We want the median and IQR for each set.
        median = np.zeros((len(data) // nreps, data.shape[1]))
        iqr_values = np.zeros((len(data) // nreps, data.shape[1]))

        for i in range(len(data) // nreps):
            start = i * nreps
            end = start + nreps

            # Get subset of data
            subset = data[start:end]

            # Remove Orders of Magnitude Differences
            subset = subset[subset["duration"] <
                            10 * subset["duration"].median()]
            subset = subset[subset["duration"] >
                            0.1 * subset["duration"].median()]

            # Remove extremely differing values that might affect the statistics
            duration_filter = (subset["duration"] < subset["duration"].quantile(0.95)) & \
                (subset["duration"] > subset["duration"].quantile(0.05))
            subset_filtered = subset[duration_filter]

            print(
                f"Processing subset {i+1}/{len(data) // nreps}: {subset_filtered}")

            # Ensure we have enough data points after filtering
            if len(subset_filtered) >= 3:  # Need at least 3 points for meaningful statistics
                median[i] = np.median(subset_filtered, axis=0)
                iqr_values[i] = iqr(subset_filtered, axis=0)
            else:
                # Fallback to original data if filtering removed too much
                median[i] = np.median(subset, axis=0)
                iqr_values[i] = iqr(subset, axis=0)

        return median, iqr_values

    gpu_med, gpu_iqr = median_and_iqr(gpu)
    cpu_med, cpu_iqr = median_and_iqr(cpu)

    print("CPU Median and IQR:")
    print(cpu_med)
    print(cpu_iqr)

    print("GPU Median and IQR:")
    print(gpu_med)
    print(gpu_iqr)

    # Placing the plots in the plane
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3))

    # Plotting CPU Data
    ax1.errorbar(n_list, cpu_med[:, 0], yerr=cpu_iqr[:, 0],
                 fmt='o-', label='CPU', color='C0')
    ax2.errorbar(n_list, cpu_med[:, 1], yerr=cpu_iqr[:, 1],
                 fmt='o-', label='CPU', color='C0')

    # Plotting GPU Data
    ax1.errorbar(n_list, gpu_med[:, 0], yerr=gpu_iqr[:, 0],
                 fmt='o-', label='GPU', color='C2')
    ax2.errorbar(n_list, gpu_med[:, 1], yerr=gpu_iqr[:, 1],
                 fmt='o-', label='GPU', color='C2')

    # Set the labels and scales for the plots
    ax1.set_xscale('log')
    # ax1.set_yscale('log')
    ax1.set_xlabel('Number of events')
    ax1.set_ylabel('Execution Time (s)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    ax2.set_xscale('log')
    # ax2.set_yscale('log')
    ax2.ticklabel_format(style='scientific',
                         axis='y', scilimits=(0, 0))
    ax2.set_xlabel('Number of events')
    ax2.set_ylabel('Tot. Energy (kWh)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig("time-and-energy.pdf")

    print("Plot saved as time-and-energy.pdf")
