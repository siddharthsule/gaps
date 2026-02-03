import os
import sys
import shutil
import argparse
import subprocess
import numpy as np
import pandas as pd
from scipy.stats import iqr
from scipy.optimize import curve_fit
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
n_list = [10000, 100000, 200000, 300000, 400000,
          500000, 600000, 700000, 800000, 900000, 1000000]

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

# Function to model f(t) = t^a, to fit to the data
def f(t, a, p, c):
    """
    Power law with offset: f(t) = a * t^p + c
    """
    return a * t ** p + c


def reduced_chi2(y_obs, y_exp, y_err, num_params):
    """
    Calculate the reduced chi-squared statistic.
    """
    chi2 = np.sum(((y_obs - y_exp) / y_err) ** 2)
    dof = len(y_obs) - num_params
    return chi2 / dof


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

    # print("GPU Emissions Data:")
    # print(gpu)

    # print("CPU Emissions Data:")
    # print(cpu)

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

            # print(
            #     f"Processing subset {i+1}/{len(data) // nreps}: {subset_filtered}")

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

    # print("CPU Median and IQR:")
    # print(cpu_med)
    # print(cpu_iqr)

    # print("GPU Median and IQR:")
    # print(gpu_med)
    # print(gpu_iqr)

    # Fit the power law function to the data
    # Convert n_list to numpy array for fitting
    n_array = np.array(n_list)

    # Fit CPU time data
    popt_cpu_time, _ = curve_fit(f, n_array, cpu_med[:, 0], p0=[1.0, 1.0, 0.0])
    a_cpu_time, p_cpu_time, c_cpu_time = popt_cpu_time

    # Fit CPU energy data
    popt_cpu_energy, _ = curve_fit(
        f, n_array, cpu_med[:, 1], p0=[1.0, 1.0, 0.0])
    a_cpu_energy, p_cpu_energy, c_cpu_energy = popt_cpu_energy

    # Fit GPU time data
    popt_gpu_time, _ = curve_fit(f, n_array, gpu_med[:, 0], p0=[1.0, 1.0, 0.0])
    a_gpu_time, p_gpu_time, c_gpu_time = popt_gpu_time

    # Fit GPU energy data
    popt_gpu_energy, _ = curve_fit(
        f, n_array, gpu_med[:, 1], p0=[1.0, 1.0, 0.0])
    a_gpu_energy, p_gpu_energy, c_gpu_energy = popt_gpu_energy

    # Try Linear fit from 100k to 1M for GPU time
    popt_linear_gpu_time, _ = curve_fit(
        lambda x, m, c: m * x + c, n_array, gpu_med[:, 0], p0=[1e-6, 0.0])
    m_linear_gpu_time, c_linear_gpu_time = popt_linear_gpu_time
    popt_linear_gpu_energy, _ = curve_fit(
        lambda x, m, c: m * x + c, n_array, gpu_med[:, 1], p0=[1e-9, 0.0])
    m_linear_gpu_energy, c_linear_gpu_energy = popt_linear_gpu_energy

    print("Fitting Power Law to Execution Time and Energy Data")
    print("f(n) = a * n^p + c, mainly interested in p exponent")
    print("")
    print("CPU:")
    print(
        f"  Execution :  p = {p_cpu_time:.2f}, reduced chi2 = {reduced_chi2(cpu_med[:, 0], f(n_array, *popt_cpu_time), cpu_iqr[:, 0], len(popt_cpu_time)):.2f}")
    print(
        f"  Energy    :  p = {p_cpu_energy:.2f}, reduced chi2 = {reduced_chi2(cpu_med[:, 1], f(n_array, *popt_cpu_energy), cpu_iqr[:, 1], len(popt_cpu_energy)):.2f}")
    print("")
    print("GPU:")
    print(
        f"  Execution :  p = {p_gpu_time:.2f}, reduced chi2 = {reduced_chi2(gpu_med[:, 0], f(n_array, *popt_gpu_time), gpu_iqr[:, 0], len(popt_gpu_time)):.2f}")
    print(
        f"  Energy    :  p = {p_gpu_energy:.2f}, reduced chi2 = {reduced_chi2(gpu_med[:, 1], f(n_array, *popt_gpu_energy), gpu_iqr[:, 1], len(popt_gpu_energy)):.2f}")
    print("")
    print("GPU (linear fit)")
    print(
        f"  Execution :  m = {m_linear_gpu_time:.6f}, reduced chi2 = {reduced_chi2(gpu_med[:, 0], m_linear_gpu_time * n_array + c_linear_gpu_time, gpu_iqr[:, 0], len(popt_linear_gpu_time)):.2f}")
    print(
        f"  Energy    :  m = {m_linear_gpu_energy:.9f}, reduced chi2 = {reduced_chi2(gpu_med[:, 1], m_linear_gpu_energy * n_array + c_linear_gpu_energy, gpu_iqr[:, 1], len(popt_linear_gpu_energy)):.2f}")
    print("")

    # Generate smooth curves for plotting
    n_smooth = np.linspace(min(n_array), max(n_array), 100)

    # Placing the plots in the plane
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3))

    # Plot CPU data
    ax1.errorbar(
        n_array, cpu_med[:, 0], yerr=cpu_iqr[:, 0], fmt='o', label='CPU', color='C0')
    ax2.errorbar(
        n_array, cpu_med[:, 1], yerr=cpu_iqr[:, 1], fmt='o', label='CPU', color='C0')

    # Plot GPU data
    ax1.errorbar(
        n_array, gpu_med[:, 0], yerr=gpu_iqr[:, 0], fmt='o', label='GPU', color='C2')
    ax2.errorbar(
        n_array, gpu_med[:, 1], yerr=gpu_iqr[:, 1], fmt='o', label='GPU', color='C2')

    # Plotting the fits
    ax1.plot(n_smooth, f(n_smooth, a_gpu_time, p_gpu_time,
             c_gpu_time), '--', color='C2', alpha=0.4)
    ax1.plot(n_smooth, f(n_smooth, a_cpu_time, p_cpu_time,
             c_cpu_time), '--', color='C0', alpha=0.4)
    ax2.plot(n_smooth, f(n_smooth, a_gpu_energy, p_gpu_energy,
             c_gpu_energy), '--', color='C2', alpha=0.4)
    ax2.plot(n_smooth, f(n_smooth, a_cpu_energy, p_cpu_energy,
             c_cpu_energy), '--', color='C0', alpha=0.4)

    # Under the plots write the power law exponents (based on fit)
    x = 700000
    ax1.text(x, f(x, a_cpu_time, p_cpu_time, c_cpu_time)+8,
             f"p = {p_cpu_time:.2f}", rotation=30*p_cpu_time, fontsize=8)
    ax1.text(x, f(x, a_gpu_time, p_gpu_time, c_gpu_time)-8,
             f"p = {p_gpu_time:.2f}", rotation=20*p_gpu_time, fontsize=8)
    ax2.text(x, f(x, a_cpu_energy, p_cpu_energy, c_cpu_energy)+0.0008,
             f"p = {p_cpu_energy:.2f}", rotation=30*p_cpu_energy, fontsize=8)
    ax2.text(x, f(x, a_gpu_energy, p_gpu_energy, c_gpu_energy)-0.0005,
             f"p = {p_gpu_energy:.2f}", rotation=15*p_gpu_energy, fontsize=8)

    # Set the labels and scales for the plots
    ax1.set_xlabel('Number of events')
    ax1.set_ylabel('Execution Time (s)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    ax2.ticklabel_format(style='scientific',
                         axis='y', scilimits=(0, 0))
    ax2.set_xlabel('Number of events')
    ax2.set_ylabel('Tot. Energy (kWh)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig("time-and-energy.pdf")

    print("Plot saved as time-and-energy.pdf")
