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

# Different Values of Events (evenly distributed on log10 scale)
# 15 points from 10^4 to 10^6
n_list = np.logspace(4, 6, 15, dtype=int).tolist()

# GPU Tuning - best value for Threads per Block
thr = 128

# Emission filenames
em_cpu = 'cpu-cluster-emissions.dat'
em_gpu = 'gpu-emissions.dat'

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
def f(t, a, p):
    """
    Pure power law: f(t) = a * t^p
    """
    return a * t ** p


def reduced_chi2(y_obs, y_exp, y_err, num_params):
    """
    Calculate the reduced chi-squared statistic.
    """
    chi2 = np.sum(((y_obs - y_exp) / y_err) ** 2)
    dof = len(y_obs) - num_params
    return chi2 / dof


def median_and_iqr(data):
    """
    Calculate the median and interquartile range (IQR) of the data.
    """
    median = np.zeros((len(data) // nreps, data.shape[1]))
    iqr_values = np.zeros((len(data) // nreps, data.shape[1]))

    for i in range(len(data) // nreps):
        start = i * nreps
        end = start + nreps
        subset = data[start:end]

        # Remove outliers using IQR method (Tukey's fences)
        duration = subset["duration"]
        q1 = duration.quantile(0.25)
        q3 = duration.quantile(0.75)
        iqr_val = q3 - q1
        lower_bound = q1 - 1.5 * iqr_val
        upper_bound = q3 + 1.5 * iqr_val

        subset = subset[(duration >= lower_bound) & (duration <= upper_bound)]

        print(
            f"Subset {i}: {len(subset)} points after outlier removal (original {nreps})")
        print(subset)

        # Use filtered data if enough points remain, otherwise use original
        if len(subset) >= 3:
            median[i] = np.median(subset, axis=0)
            iqr_values[i] = iqr(subset, axis=0)
        else:
            median[i] = np.median(data[start:end], axis=0)
            iqr_values[i] = iqr(data[start:end], axis=0)

        print(
            f"Subset {i}: median={median[i][0]:.6f}, error_bars=[{median[i][0] - iqr_values[i][0]/2:.6f}, {median[i][0] + iqr_values[i][0]/2:.6f}]")

    return median, iqr_values


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

    # Convert energy from kWh to Wh
    gpu['energy_consumed'] *= 1000
    cpu['energy_consumed'] *= 1000

    # print("GPU Emissions Data:")
    # print(gpu)

    # print("CPU Emissions Data:")
    # print(cpu)

    gpu_med, gpu_iqr = median_and_iqr(gpu)
    cpu_med, cpu_iqr = median_and_iqr(cpu)

    # Convert n_list to numpy array for plotting
    n_array = np.array(n_list)

    # Filter data for 25k events or more
    mask_25k = n_array >= 50000
    n_fit = n_array[mask_25k]
    cpu_med_fit = cpu_med[mask_25k]
    gpu_med_fit = gpu_med[mask_25k]

    # Fit linear in log space (m*x + c) for log-log plots - using only 25k+ data
    # This means log(y) = m*log(x) + c, so y = exp(c) * x^m
    log_n_fit = np.log(n_fit)

    # GPU time log fit
    log_gpu_time_fit = np.log(gpu_med_fit[:, 0])
    m_gpu_time, c_gpu_time = np.polyfit(log_n_fit, log_gpu_time_fit, 1)

    # GPU energy log fit
    log_gpu_energy_fit = np.log(gpu_med_fit[:, 1])
    m_gpu_energy, c_gpu_energy = np.polyfit(log_n_fit, log_gpu_energy_fit, 1)

    # Generate smooth curves for plotting (only in the fit range)
    n_smooth = np.linspace(min(n_fit), max(n_fit), 100)

    # Create 1x2 subplot layout
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3))

    # Log-Log scale with linear fits in log space
    # Plot CPU data
    ax1.errorbar(
        n_array, cpu_med[:, 0], yerr=cpu_iqr[:, 0]/2, fmt='o', label='CPU', color='C0')

    ax2.errorbar(
        n_array, cpu_med[:, 1], yerr=cpu_iqr[:, 1]/2, fmt='o', label='CPU', color='C0')

    # Plot GPU data
    ax1.errorbar(
        n_array, gpu_med[:, 0], yerr=gpu_iqr[:, 0]/2, fmt='o', label='GPU', color='C2')

    ax2.errorbar(
        n_array, gpu_med[:, 1], yerr=gpu_iqr[:, 1]/2, fmt='o', label='GPU', color='C2')

    # Plot fit lines
    ax1.plot(n_smooth, np.exp(m_gpu_time * np.log(n_smooth) + c_gpu_time),
             '--', color='C2', alpha=0.5, linewidth=1.5)

    ax2.plot(n_smooth, np.exp(m_gpu_energy * np.log(n_smooth) + c_gpu_energy),
             '--', color='C2', alpha=0.5, linewidth=1.5)

    # Add tilted text annotations for gradients
    x = 0.7
    ax1.text(x, 0.65, f'$p_{{\\mathrm{{GPU}}}} = {m_gpu_time:.2f}$',
             transform=ax1.transAxes, rotation=30, va='top', ha='left', color='C2', fontsize=9)

    ax2.text(x, 0.65, f'$p_{{\\mathrm{{GPU}}}} = {m_gpu_energy:.2f}$',
             transform=ax2.transAxes, rotation=30, va='top', ha='left', color='C2', fontsize=9)

    # Set labels and formatting
    ax1.set_xlabel('Number of events')
    ax1.set_ylabel('Execution Time (s)')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.legend()
    ax1.grid(True, alpha=0.2, which='both')

    ax2.set_xlabel('Number of events')
    ax2.set_ylabel('Tot. Energy (Wh)')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.legend()
    ax2.grid(True, alpha=0.2, which='both')

    fig.tight_layout()
    fig.savefig("time-and-energy.pdf")

    print("Plot saved as time-and-energy.pdf")
