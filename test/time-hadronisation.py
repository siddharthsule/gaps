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

# 10 points from 10^4 to 10^6
n_list = np.logspace(4, 6, 10, dtype=int).tolist()

# GPU Tuning - best value for Threads per Block
thr = 128

time_cpu = 'cpu-time.dat'
time_gpu = 'gpu-time.dat'

# Column order in the timing data files
time_cols = ['me_time', 'shower_time', 'hadronisation_time',
             'decay_time', 'analysis_time', 'total_time']

# Run gaps -r compare once to compile both codes
if option != "plot":
    os.chdir("..")
    subprocess.run("./rungaps -r compare", shell=True)
    os.chdir("test")

# ------------------------------------------------------------------------------
# Run simulations

if option == 'cpu' or option == 'gpu':
    os.chdir("..")

    time_file = time_gpu if option == 'gpu' else time_cpu
    if os.path.exists(time_file):
        os.remove(time_file)

    for n in n_list:
        command = f"./rungaps -p LEP -nlo -hadronise -n {n} --skip-analysis --no-compile"
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

        t = data.columns.get_loc("total_time")
        print(
            f"Subset {i}: median total={median[i][t]:.4f}s, IQR={iqr_values[i][t]:.4f}s")

    return median, iqr_values


if option == 'plot':
    # Every timing file present, as (title, colour, columns, median, IQR)
    candidates = [('CPU', time_cpu, nreps * ncpu, 'C0'),
                  ('GPU', time_gpu, nreps, 'C2')]
    runs = []

    for title, time_file, n_per_group, color in candidates:
        if not os.path.exists("../" + time_file):
            print(f"{time_file} not found, skipping {title}")
            continue

        data = pd.read_csv("../" + time_file, header=None,
                           delimiter=',', names=time_cols)

        # Hadronisation and hadron decays as one component
        data["had_decay_time"] = (data["hadronisation_time"]
                                  + data["decay_time"])

        med, iqr_vals = median_and_iqr(data, n_per_group)
        runs.append((title, color, data.columns, med, iqr_vals))

    if not runs:
        print(f"Error: neither {time_cpu} nor {time_gpu} found.")
        exit(1)

    n_array = np.array(n_list)
    n_groups = len(n_list)
    x = np.arange(n_groups)

    comp_labels = ['ME', 'Shower', 'Hadronisation', 'Decays', 'Analysis']
    comp_colors = ['C0', 'C1', 'C2', 'C3', 'C4']

    tick_labels = [f'$10^{{{np.log10(n):.1f}}}$' for n in n_array]

    fig, axes = plt.subplots(1, len(runs), figsize=(5 * len(runs), 4),
                             sharey=True, squeeze=False)
    axes = axes[0]

    for ax, (title, _, _, med, iqr_vals) in zip(axes, runs):
        for k, (label, color) in enumerate(zip(comp_labels, comp_colors)):
            ax.bar(x, med[:, k], bottom=med[:, :k].sum(axis=1), color=color,
                   label=label)

        # # Error bars on the total time
        # ax.errorbar(x, med[:, 5], yerr=iqr_vals[:, 5] / 2,
        #             fmt='none', color='black', capsize=2, linewidth=0.8)

        ax.set_title(title)
        ax.set_xlabel('Number of events')
        ax.set_xticks(x)
        ax.set_xticklabels(tick_labels, rotation=45, ha='right', fontsize=7)
        ax.grid(True, alpha=0.2, axis='y')
        ax.set_yscale('log')
        # ax.set_ylim(0, None)

    axes[0].set_ylabel('Time (s)')
    axes[0].legend(loc='upper left', fontsize=8)

    fig.tight_layout()
    fig.savefig("time-hadronisation.pdf")
    print("Plot saved as time-hadronisation.pdf")

    # --------------------------------------------------------------------------
    # Log-log plots: time for the shower and for hadronisation with hadron
    # decays, one series per timing file

    # Fit power law in log space (log(y) = m*log(x) + c => y = exp(c) * x^m)
    # using only 50k+ events, as in time-and-energy.py
    mask_fit = n_array >= 50000
    log_n_fit = np.log(n_array[mask_fit])
    n_smooth = np.linspace(
        n_array[mask_fit].min(), n_array[mask_fit].max(), 100)

    # Components to show, ('label', column); Shower on top, Hadronisation and
    # Decays below
    loglog_comps = [('Shower', 'shower_time'),
                    ('Hadronisation + Decays', 'had_decay_time')]

    fig2, axes = plt.subplots(1, 2, figsize=(10, 3), sharey=True)

    for (label, column), ax in zip(loglog_comps, axes):
        fit_lines = []

        for title, color, columns, med, iqr_vals in runs:
            k = columns.get_loc(column)

            # Data points
            ax.errorbar(n_array, med[:, k], yerr=iqr_vals[:, k] / 2,
                        fmt='o', label=title, color=color)

            # Power-law fit (50k+ events)
            m, c = np.polyfit(log_n_fit, np.log(med[mask_fit, k]), 1)
            ax.plot(n_smooth, np.exp(m * np.log(n_smooth) + c),
                    '--', color=color, alpha=0.5, linewidth=1.5)
            fit_lines.append(f'$p_{{\\mathrm{{{title}}}}} = {m:.2f}$')

        # Annotate the fitted gradients (power-law exponents), bottom right
        ax.text(0.95, 0.05, '\n'.join(fit_lines),
                transform=ax.transAxes, va='bottom', ha='right', fontsize=9,
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

        ax.set_title(label)
        ax.set_xlabel('Number of events')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.grid(True, alpha=0.2, which='both')

    axes[0].legend(loc='upper left', fontsize=8)
    axes[0].set_ylabel('Time (s)')

    fig2.tight_layout()
    fig2.savefig("time-hadronisation-loglog.pdf")
    print("Plot saved as time-hadronisation-loglog.pdf")
