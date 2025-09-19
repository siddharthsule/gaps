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
# Run Tune with and without partitioning

if not args.plotting_only:

    # Move to base directory, ../
    os.chdir('../')

    # Run with Partitioning
    subprocess.run(
        "./rungaps -p LHC -nlo -codecarbon -r tune -do_partitioning yes", shell=True)
    os.rename("gpu-time.dat", "gpu-time-part.dat")

    # Run without Partitioning
    subprocess.run(
        "./rungaps -p LHC -nlo -codecarbon -r tune -do_partitioning no", shell=True)
    os.rename("gpu-time.dat", "gpu-time-nopart.dat")

    # Return to original directory
    os.chdir('test')

# ------------------------------------------------------------------------------
# Analysis and Plotting Code


def get_data_and_iqr(df):

    # Make an empty dataframe to store the values
    dat = pd.DataFrame(np.zeros((len(nev), len(threads_per_block))),
                       index=nev, columns=threads_per_block)

    dat_iqr = pd.DataFrame(np.zeros(
        (len(nev), len(threads_per_block))), index=nev, columns=threads_per_block)

    # For every nreps values, calculate the median and iqr
    for i in range(0, len(df), args.nreps):
        # Get the values
        values = df[i:i+args.nreps]

        # Remove values orders of magnitude larger/smaller than the median
        median = values.median()
        values = values[values < 10 * median]
        values = values[values > 0.1 * median]

        # Remove extremely differing values that might affect the statistics
        values = values[values < values.quantile(0.95)]
        values = values[values > values.quantile(0.05)]

        # Remove any negative values
        values = values[values > 0]

        # Calculate the median and iqr
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
    print(dat_iqr)

    # Return
    return dat, dat_iqr


def plot_data(dat, dat_iqr, ax1, ax2, color='C2', label_prefix=''):

    for i in range(len(nev)):
        # Define different format styles for each series
        fmt_styles = ['o-', 's-', '^-']  # circle, square, triangle

        # Create label with prefix if provided
        label = f'{label_prefix}, ${nev_str[i]}$'

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
nopart = pd.read_csv('../gpu-time-nopart.dat', header=None, delimiter=',')
part = pd.read_csv('../gpu-time-part.dat', header=None, delimiter=',')

# Isolate the last column (assuming that's the time data we want)
nopart_data = nopart.iloc[:, -1]
part_data = part.iloc[:, -1]

# Plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 6), sharex=True)

# Process and plot both datasets
datasets = [
    (nopart_data, 'C2', 'No Partitioning'),
    (part_data, '#92D050', 'Partitioning')
]

# Store all data for setting consistent y-limits
all_data = []

for df, color, label_prefix in datasets:
    dat, dat_iqr = get_data_and_iqr(df)
    all_data.append(dat)
    plot_data(dat, dat_iqr, ax1, ax2,
              color=color, label_prefix=label_prefix)

# Set consistent y-limits based on all data
if all_data:
    # Upper part - for 1,000,000 events (row 2)
    all_row2_min = min(np.min(data.iloc[2, :]) for data in all_data)
    all_row2_max = max(np.max(data.iloc[2, :]) for data in all_data)
    ax1.set_ylim(all_row2_min - 2, all_row2_max + 5)

    # Lower part - for 10,000 and 100,000 events (rows 0 and 1)
    all_lower_min = min(min(np.min(data.iloc[0, :]), np.min(
        data.iloc[1, :])) for data in all_data)
    all_lower_max = max(max(np.max(data.iloc[0, :]), np.max(
        data.iloc[1, :])) for data in all_data)
    ax2.set_ylim(all_lower_min - 1, all_lower_max + 1)

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

# Get handles and labels from ax1 for legend
handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles, labels, loc='upper right', ncol=2)

ax1.grid(True, alpha=0.3)
ax2.grid(True, alpha=0.3)

# Set x-axis to log base 2 scale
ax2.set_xscale('log', base=2)
ax2.set_xticks(threads_per_block)
ax2.set_xticklabels(threads_per_block)

plt.tight_layout()
plt.savefig('gpu-tune.pdf')
