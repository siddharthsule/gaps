import os
import argparse
import subprocess
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import iqr

parser = argparse.ArgumentParser(description='Tune and Partitioning Tests')
parser.add_argument('--nreps', type=int, default=10,
                    help='Number of repetitions')
parser.add_argument('--plotting-only', action='store_true',
                    help='Only generate plots')
parser.add_argument('--no-plotting', action='store_true',
                    help='Run tests but do not generate plots')
args = parser.parse_args()

# ------------------------------------------------------------------------------
# Consistent Settings as rungaps

# rows = nev values
# cols = threads per block
nev = [10000, 100000, 1000000]
nev_str = ["10^4", "10^5", "10^6"]
threads_per_block = [32, 64, 128, 256, 512]

# ------------------------------------------------------------------------------
# Run Tune with partitioning

if not args.plotting_only:

    # Move to base directory, ../
    os.chdir('../')

    # Run with Partitioning
    subprocess.run(
        "./rungaps -p LHC -nlo -codecarbon -r tune -do_partitioning yes", shell=True)
    os.rename("gpu-time.dat", "gpu-time-part.dat")

    # Return to original directory
    os.chdir('test')

# ------------------------------------------------------------------------------
# Analysis and Plotting Code


def get_data_and_iqr(df):
    """
    Calculate the median and interquartile range (IQR) of the data.
    """
    # Make an empty dataframe to store the values
    dat = pd.DataFrame(np.zeros((len(nev), len(threads_per_block))),
                       index=nev, columns=threads_per_block)

    dat_iqr = pd.DataFrame(np.zeros(
        (len(nev), len(threads_per_block))), index=nev, columns=threads_per_block)

    # For every nreps values, calculate the median and iqr
    for i in range(0, len(df), args.nreps):
        # Get the values
        values = df[i:i+args.nreps]

        # Remove outliers using IQR method (Tukey's fences)
        q1 = values.quantile(0.25)
        q3 = values.quantile(0.75)
        iqr_val = q3 - q1
        lower_bound = q1 - 1.5 * iqr_val
        upper_bound = q3 + 1.5 * iqr_val

        values_filtered = values[(values >= lower_bound)
                                 & (values <= upper_bound)]

        # Use filtered data if enough points remain, otherwise use original
        if len(values_filtered) >= 3:
            median = values_filtered.median()
            iqr_value = iqr(values_filtered)
        else:
            median = values.median()
            iqr_value = iqr(values)

        # Determine the row and column indices
        # Total combinations = len(nev) * len(threads_per_block)
        combination_index = i // args.nreps
        row = combination_index // len(threads_per_block)
        col = combination_index % len(threads_per_block)

        # Assign the calculated values to the dataframes
        dat.iloc[row, col] = median
        dat_iqr.iloc[row, col] = iqr_value

    # Print the results
    print(dat)
    print(dat_iqr/2)

    # Return
    return dat, dat_iqr/2


def plot_data(dat, dat_iqr, ax1, ax2, color='C2', label_prefix=''):

    for i in range(len(nev)):
        # Define different format styles for each series
        fmt_styles = ['o-', 's-', '^-']  # circle, square, triangle

        # Create label with prefix if provided
        label = f'{label_prefix}, ${nev_str[i]}$ Events'

        # Plot on both axes
        ax1.errorbar(threads_per_block, dat.iloc[i, :],
                     yerr=dat_iqr.iloc[i, :],
                     fmt=fmt_styles[i],
                     label=label,
                     color=color)

        ax2.errorbar(threads_per_block, dat.iloc[i, :],
                     yerr=dat_iqr.iloc[i, :],
                     fmt=fmt_styles[i],
                     label=label,
                     color=color)


if args.no_plotting:
    exit(0)

# Set the style
mpl.rc_file("mplstyleerc")

# Get the data for GPU
part = pd.read_csv('../gpu-time-part.dat', header=None, delimiter=',')

# Isolate the last column (assuming that's the time data we want)
part_data = part.iloc[:, -1]

# Plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 5), sharex=True)

# Process and plot data
dat, dat_iqr = get_data_and_iqr(part_data)
plot_data(dat, dat_iqr, ax1, ax2, color='C2', label_prefix='Partitioning')

# Set fixed y-limits for both plots
ax1.set_ylim(53, 63)
ax1.set_yticks(np.arange(53, 64, 2))

ax2.set_ylim(2, 12)
ax2.set_yticks(np.arange(2, 13, 2))

# Hide the spines between ax1 and ax2
ax1.spines.bottom.set_visible(False)
ax2.spines.top.set_visible(False)
ax1.tick_params(bottom=False, top=False, labeltop=False, labelbottom=False)
ax2.xaxis.tick_bottom()

# Add break lines
d = .5  # proportion of vertical to horizontal extent of the slanted line
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
              linestyle="none", color='k', mec='k', mew=1, clip_on=False)
ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)

# Labels and formatting
ax2.set_xlabel('Threads per block, $N_T$')
fig.supylabel('Execution time (s)', fontsize=11)

ax1.legend(loc='upper right')

ax1.grid(True, alpha=0.3)
ax2.grid(True, alpha=0.3)

# Set x-axis to log base 2 scale
ax2.set_xscale('log', base=2)
ax2.set_xticks(threads_per_block)
ax2.set_xticklabels(threads_per_block)

fig.tight_layout()
fig.savefig('gpu-tune.pdf')
