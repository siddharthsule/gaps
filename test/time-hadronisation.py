import os
import argparse
import subprocess
import numpy as np
import pandas as pd
from scipy.stats import iqr
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc_file("mplstyleerc")

parser = argparse.ArgumentParser(
    description='Run Hadronisation Timing Analysis')
parser.add_argument('option', choices=[
                    'cpu', 'gpu', 'plot'], help='Select the option to run')
parser.add_argument('--nreps', type=int, default=10,
                    help='Number of repetitions for each run (default: 10)')
parser.add_argument('--ncpu', type=int, default=1,
                    help='Number of CPU threads (default: 1)')
args = parser.parse_args()
option = args.option
nreps = args.nreps
ncpu = args.ncpu

# 15 points from 10^4 to 10^6
n_list = np.logspace(4, 6, 15, dtype=int).tolist()

# GPU Tuning - best value for Threads per Block
thr = 128

time_cpu = 'cpu-time.dat'
time_gpu = 'gpu-time.dat'

# Column order in the timing data files
time_cols = ['me_time', 'shower_time',
             'hadronisation_time', 'analysis_time', 'total_time']

# ------------------------------------------------------------------------------
# Run simulations

if option == 'cpu' or option == 'gpu':
    os.chdir("..")

    time_file = time_gpu if option == 'gpu' else time_cpu
    if os.path.exists(time_file):
        os.remove(time_file)

    for n in n_list:
        command = f"./rungaps -p LEP -nlo -hadronise -n {n} -do_partitioning no"
        if option == 'gpu':
            command += f" -t {thr}"
        elif option == 'cpu':
            command += f" -r cpu-cluster -ncpu {ncpu}"

        for i in range(nreps):
            subprocess.run(command, shell=True)

    os.chdir("test")
    print(f"Completed {option} analysis. Results saved to {time_file}")
    exit(0)


# ------------------------------------------------------------------------------
# Plotting

def median_and_iqr(data, n_per_group):
    n_groups = len(data) // n_per_group
    median = np.zeros((n_groups, data.shape[1]))
    iqr_values = np.zeros((n_groups, data.shape[1]))

    for i in range(n_groups):
        start = i * n_per_group
        end = start + n_per_group
        subset = data[start:end]

        total = subset["total_time"]
        q1 = total.quantile(0.25)
        q3 = total.quantile(0.75)
        iqr_val = q3 - q1
        lower_bound = q1 - 1.5 * iqr_val
        upper_bound = q3 + 1.5 * iqr_val
        subset_filtered = subset[(total >= lower_bound)
                                 & (total <= upper_bound)]

        print(
            f"Subset {i}: {len(subset_filtered)} points after outlier removal (original {n_per_group})")

        if len(subset_filtered) >= 3:
            median[i] = np.median(subset_filtered, axis=0)
            iqr_values[i] = iqr(subset_filtered, axis=0)
        else:
            median[i] = np.median(data[start:end], axis=0)
            iqr_values[i] = iqr(data[start:end], axis=0)

        print(
            f"Subset {i}: median total={median[i][-1]:.4f}s, IQR={iqr_values[i][-1]:.4f}s")

    return median, iqr_values


if option == 'plot':
    if not os.path.exists("../" + time_cpu):
        print(f"Error: {time_cpu} not found.")
        exit(1)
    if not os.path.exists("../" + time_gpu):
        print(f"Error: {time_gpu} not found.")
        exit(1)

    gpu = pd.read_csv("../" + time_gpu, header=None,
                      delimiter=',', names=time_cols)
    cpu = pd.read_csv("../" + time_cpu, header=None,
                      delimiter=',', names=time_cols)

    gpu_med, gpu_iqr_vals = median_and_iqr(gpu, nreps)
    cpu_med, cpu_iqr_vals = median_and_iqr(cpu, nreps * ncpu)

    n_array = np.array(n_list)
    n_groups = len(n_list)
    x = np.arange(n_groups)

    comp_labels = ['ME', 'Shower', 'Hadronisation', 'Analysis']
    comp_colors = ['C0', 'C1', 'C2', 'C3']

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), sharey=True)

    for k, (label, color) in enumerate(zip(comp_labels, comp_colors)):
        bottom_cpu = cpu_med[:, :k].sum(
            axis=1) if k > 0 else np.zeros(n_groups)
        bottom_gpu = gpu_med[:, :k].sum(
            axis=1) if k > 0 else np.zeros(n_groups)

        ax1.bar(x, cpu_med[:, k], bottom=bottom_cpu, color=color, label=label)
        ax2.bar(x, gpu_med[:, k], bottom=bottom_gpu, color=color, label=label)

    # # Error bars on the total time
    # ax1.errorbar(x, cpu_med[:, 4], yerr=cpu_iqr_vals[:, 4] / 2,
    #              fmt='none', color='black', capsize=2, linewidth=0.8)
    # ax2.errorbar(x, gpu_med[:, 4], yerr=gpu_iqr_vals[:, 4] / 2,
    #              fmt='none', color='black', capsize=2, linewidth=0.8)

    tick_labels = [f'$10^{{{np.log10(n):.1f}}}$' for n in n_array]

    for ax, title in [(ax1, 'CPU'), (ax2, 'GPU')]:
        ax.set_title(title)
        ax.set_xlabel('Number of events')
        ax.set_xticks(x)
        ax.set_xticklabels(tick_labels, rotation=45, ha='right', fontsize=7)
        ax.grid(True, alpha=0.2, axis='y')

    ax1.set_ylabel('Time (s)')
    ax1.legend(loc='upper left', fontsize=8)

    # ax1.set_ylim(0, None)
    # ax2.set_ylim(0, None)

    ax1.set_yscale('log')
    ax2.set_yscale('log')

    fig.tight_layout()
    fig.savefig("time-hadronisation.pdf")
    print("Plot saved as time-hadronisation.pdf")

    # --------------------------------------------------------------------------
    # Log-log plots: CPU vs GPU time per component (ME, Shower, Had., Analysis)

    # Fit power law in log space (log(y) = m*log(x) + c => y = exp(c) * x^m)
    # using only 50k+ events, as in time-and-energy.py
    mask_fit = n_array >= 50000
    log_n_fit = np.log(n_array[mask_fit])
    n_smooth = np.linspace(n_array[mask_fit].min(), n_array[mask_fit].max(), 100)

    # Components to show (skip ME); ('label', column index in time_cols)
    loglog_comps = [('Shower', 1), ('Hadronisation', 2), ('Observables', 3)]

    fig2, axes = plt.subplots(1, 3, figsize=(12, 3))

    for (label, k), ax in zip(loglog_comps, axes):
        # CPU and GPU data points
        ax.errorbar(n_array, cpu_med[:, k], yerr=cpu_iqr_vals[:, k] / 2,
                    fmt='o', label='CPU', color='C0')
        ax.errorbar(n_array, gpu_med[:, k], yerr=gpu_iqr_vals[:, k] / 2,
                    fmt='o', label='GPU', color='C2')

        # Power-law fits (50k+ events)
        m_cpu, c_cpu = np.polyfit(log_n_fit, np.log(cpu_med[mask_fit, k]), 1)
        m_gpu, c_gpu = np.polyfit(log_n_fit, np.log(gpu_med[mask_fit, k]), 1)

        ax.plot(n_smooth, np.exp(m_cpu * np.log(n_smooth) + c_cpu),
                '--', color='C0', alpha=0.5, linewidth=1.5)
        ax.plot(n_smooth, np.exp(m_gpu * np.log(n_smooth) + c_gpu),
                '--', color='C2', alpha=0.5, linewidth=1.5)

        # Annotate the fitted gradients (power-law exponents)
        ax.text(0.05, 0.95, f'$p_{{\\mathrm{{CPU}}}} = {m_cpu:.2f}$',
                transform=ax.transAxes, va='top', ha='left',
                color='C0', fontsize=9)
        ax.text(0.05, 0.85, f'$p_{{\\mathrm{{GPU}}}} = {m_gpu:.2f}$',
                transform=ax.transAxes, va='top', ha='left',
                color='C2', fontsize=9)

        ax.set_title(label)
        ax.set_xlabel('Number of events')
        ax.set_ylabel('Time (s)')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.grid(True, alpha=0.2, which='both')

    axes[0].legend(loc='lower right', fontsize=8)

    fig2.tight_layout()
    fig2.savefig("time-hadronisation-loglog.pdf")
    print("Plot saved as time-hadronisation-loglog.pdf")
