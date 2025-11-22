import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
import pickle
import argparse

mpl.rc_file("mplstyleerc")

# ------------------------------------------------------------------------------
# Run Gaps


def run(nev, cores):
    os.chdir('../')
    subprocess.run(
        f"./rungaps -p LHC -nlo -r cpu-cluster -ncpu {cores} -n {nev}", shell=True)
    os.chdir('test')

# ------------------------------------------------------------------------------
# Data Reading


def read_data(filename):
    """
    Given a filename, read the data and return a pandas DataFrame.
    """

    print("Reading data from:", filename)

    # Get the total number of lines in the file
    with open(filename, 'r') as file:
        total_lines = sum(1 for line in file)

    # Importing the dataset, skipping the first and last 100 lines
    dataset = pd.read_csv(filename, skiprows=100, nrows=total_lines - 200)
    dataset.columns = ["pid_initial", "pid_final", "eta", "ratio"]

    # Convert columns to appropriate data types
    dataset['pid_initial'] = dataset['pid_initial'].astype(int)
    dataset['pid_final'] = dataset['pid_final'].astype(int)

    # Convert the 'ratio' column to numeric values, forcing errors to NaN
    dataset.loc[:, 'ratio'] = pd.to_numeric(dataset['ratio'], errors='coerce')

    return dataset

# ------------------------------------------------------------------------------
# Physics Processing


def get_data_subset(dataset, pid_i, pid_f):
    """
    Given a dataset, return a subset of the data where pid_initial and pid_final
    match the given values.
    """

    return dataset[(dataset['pid_initial'] == pid_i)
                   & (dataset['pid_final'] == pid_f)]


def apply_cuts_to_subset(subset, fl_i, fl_f):

    # g -> bb known to be large
    if abs(fl_i) == 5 and fl_f == 21:
        subset = subset[subset['ratio'] < 100000]
    # similar case for g -> cc
    elif abs(fl_i) == 4 and fl_f == 21:
        subset = subset[subset['ratio'] < 20000]
    # apply strict upper limit for the rest
    elif abs(fl_i) in [1, 2, 3] and fl_f == 21:
        subset = subset[subset['ratio'] < 10000]
    # apply strict upper limit for the rest
    else:
        subset = subset[subset['ratio'] < 50]

    # Remove any rows with NaN values
    subset = subset.dropna()

    # Remove any vals with inf values
    subset = subset[~subset['ratio'].isin([np.inf, -np.inf])]

    # Remove any vals with eta < 1e-5
    subset = subset[subset['eta'] > 1e-5]

    # Remove any vals with ratio < 1e-5
    subset = subset[subset['ratio'] > 1e-5]

    return subset


def get_subsets(dataset):
    """
    Get all valid subsets from the dataset for different flavour combinations.
    """

    print("Generating subsets from dataset...")

    subsets = []
    subset_flavours = []

    # For all valid flavour combinations
    flavours = np.array([21, 1, 2, 3, 4, 5, -1, -2, -3, -4, -5])
    for i in range(len(flavours)):
        for f in range(len(flavours)):

            fl_i = flavours[i]
            fl_f = flavours[f]

            if (fl_i != fl_f) and not ((fl_i == 21) or (fl_f == 21)):
                continue

            # Get the data subset and apply cuts
            subset = get_data_subset(dataset, fl_i, fl_f)
            subset = apply_cuts_to_subset(subset, fl_i, fl_f)

            # Store the subset if it has valid data
            if len(subset) > 0:
                subsets.append(subset)
                subset_flavours.append((fl_i, fl_f))

    return subsets, subset_flavours


# ------------------------------------------------------------------------------
# Heatmap Generation and Storage

def generate_and_store_heatmap(subsets, subset_flavours):
    """
    For all the subsets from a dataset of file_id, generate and store the
    log10-log10 heatmap data ONLY
    """

    for i in range(len(subsets)):
        subset = subsets[i]
        fl_i, fl_f = subset_flavours[i]

        # Generate the log10-log10 heatmap data
        h, xe, ye = np.histogram2d(
            np.log10(subset['eta']), np.log10(subset['ratio']), bins=100)

        # Also, store y99 and y9999
        y99 = np.percentile(subset['ratio'], 99)
        y9999 = np.percentile(subset['ratio'], 99.99)

        # Map the heatmap data
        heatmap_data = {
            'h': h,
            'xe': xe,
            'ye': ye,
            'fl_i': fl_i,
            'fl_f': fl_f,
            'y99': y99,
            'y9999': y9999
        }

        # Create directory for heatmap data
        os.makedirs('heatmap-data', exist_ok=True)

        # Save the data as a pkl
        filename = f'heatmap-data/heatmap_{fl_i}_{fl_f}.pkl'
        with open(filename, 'wb') as f:
            pickle.dump(heatmap_data, f)

# ------------------------------------------------------------------------------
# Heatmap Loading and Combining


def read_heatmap(fl_i, fl_f):
    """
    Load heatmap data for the given flavor combination.
    """

    print(f"Loading heatmap for flavours: {fl_i}, {fl_f}")

    # Load the single heatmap file for this flavor combination
    filename = f'heatmap-data/heatmap_{fl_i}_{fl_f}.pkl'
    try:
        with open(filename, 'rb') as f:
            heatmap_data = pickle.load(f)
        return heatmap_data
    except FileNotFoundError:
        print(f"Heatmap data not found: {filename}")
        return None

# ------------------------------------------------------------------------------
# Heatmap Plotting


def pids_to_names(fl_i, fl_f):

    pid_names = {
        21: 'g',
        1: 'd',
        2: 'u',
        3: 's',
        4: 'c',
        5: 'b',
        -1: r'\bar{d}',
        -2: r'\bar{u}',
        -3: r'\bar{s}',
        -4: r'\bar{c}',
        -5: r'\bar{b}'
    }

    name_i = pid_names.get(fl_i, str(fl_i))
    name_f = pid_names.get(fl_f, str(fl_f))

    return f'${name_i} \\to {name_f}$'


def plot_heatmap(heatmap_data):
    """
    Generate a combined heatmap from individual heatmap files.
    """

    print(
        f"Plotting heatmap for flavours: {heatmap_data['fl_i']}, {heatmap_data['fl_f']}")

    # Extract the data
    h, xe, ye = heatmap_data['h'], heatmap_data['xe'], heatmap_data['ye']
    fl_i, fl_f = heatmap_data['fl_i'], heatmap_data['fl_f']

    fig, ax = plt.subplots(figsize=(4, 3))

    # Plot heatmap of log10(ratio) vs log10(eta)
    extent = [xe[0], xe[-1], ye[0], ye[-1]]
    ax.imshow(h.T, extent=extent, origin='lower',
              aspect='auto', cmap='cividis')

    ax.text(0.05, 0.15, pids_to_names(fl_i, fl_f), transform=ax.transAxes,
            color='white', fontsize=12, fontweight='bold',
            verticalalignment='top', horizontalalignment='left')

    if fl_i == 21 and fl_f in [1, 2]:

        y9999 = heatmap_data['y9999']
        y9999 = int(np.ceil(y9999 / 5) * 5)

        xvals = np.linspace(xe[0], xe[-1], num=100)
        fit = y9999 * np.sqrt(10 ** xvals)  # Use 10^x instead of exp(x)
        log10fit = np.log10(fit)
        ax.plot(xvals, log10fit, color="white")

        ax.text(0.05, 0.08, f'$f(\\eta) = {y9999}\\sqrt{{\\eta}}$', transform=ax.transAxes,
                color='white', fontsize=8, fontweight='bold',
                verticalalignment='top', horizontalalignment='left')

    else:

        y99 = heatmap_data['y99']
        y9999 = heatmap_data['y9999']

        # Round them up to the nearest 1
        y99 = int(np.ceil(y99))
        y9999 = int(np.ceil(y9999))

        if abs(fl_i) > 2 and fl_i != 21 and fl_f == 21:
            # Round them up to the nearest 50
            y99 = int(np.ceil(y99 / 50) * 50)
            y9999 = int(np.ceil(y9999 / 50) * 50)

        # Convert bin_midpoints from ln to log10 scale
        xvals = np.linspace(xe[0], xe[-1], num=100)

        fit = y9999 * (10 ** xvals)  # Use 10^x instead of exp(x)
        fit[fit < y99] = y99
        log10fit = np.log10(fit)
        ax.plot(xvals, log10fit, color="white")

        ax.text(0.05, 0.08, f'$f(\\eta) = \\max({y99}, {y9999}\\eta)$', transform=ax.transAxes,
                color='white', fontsize=8, fontweight='bold',
                verticalalignment='top', horizontalalignment='left')

    # Set labels for the axes
    ax.set_xlabel(r'$\eta$')
    ax.set_ylabel(r'$\text{PDF Ratio}$')

    # Only Keep whole number x ticks
    xticks = ax.get_xticks()
    xticks = xticks[xticks == np.floor(xticks)]
    ax.set_xticks(xticks)
    # Swap out x tick labels
    xtick_labels = [10**tick for tick in xticks]
    ax.set_xticklabels(xtick_labels)

    # Only Keep whole number y ticks
    yticks = ax.get_yticks()
    yticks = yticks[yticks == np.floor(yticks)]
    ax.set_yticks(yticks)
    # Swap out y tick labels
    ytick_labels = [10**tick for tick in yticks]
    ax.set_yticklabels(ytick_labels)

    # Tighten up xlim and ylim
    ax.set_xlim(left=xe[0], right=np.log10(1.))
    ax.set_ylim(bottom=ye[0], top=ye[-1])

    # Tighten layout to maximize space
    fig.tight_layout()

    # Make a directory to store the plots
    os.makedirs('pdfratio-plots', exist_ok=True)

    # Save the figure
    filename = f'pdfratio-plots/ratio_{fl_i}_{fl_f}.pdf'
    fig.savefig(filename)

    # Close the figure to free up memory
    plt.close(fig)


# ------------------------------------------------------------------------------
# Main execution section

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Run Gaps and plot PDF ratio heatmaps.")
    parser.add_argument('--nev', type=int, default=10000,
                        help='Number of events to generate (default: 10000)')
    parser.add_argument('--cores', type=int, default=4,
                        help='Number of CPU cores to use (default: 4)')
    args = parser.parse_args()

    # Run Gaps to generate data
    nev = args.nev
    cores = args.cores
    run(nev, cores)

    # Files
    files = os.listdir('../')
    in_files = [f for f in files if f.startswith(
        'temp-') and f.endswith('.dat')]

    # Big Dataset - collect all datasets first
    datasets = []

    # Process Each File
    for filename in in_files:

        # Read the data
        file_dataset = read_data(os.path.join('../', filename))

        # Add to the list of datasets
        datasets.append(file_dataset)

    # Concatenate all datasets at once
    dataset_full = pd.concat(datasets, ignore_index=True)

    # Get the subsets
    subsets, subset_flavours = get_subsets(dataset_full)

    # Gen and Store Heatmap
    generate_and_store_heatmap(subsets, subset_flavours)

    # Generate Combined Heatmaps
    flavours = np.array([21, 1, 2, 3, 4, 5, -1, -2, -3, -4, -5])
    for i in range(len(flavours)):
        for f in range(len(flavours)):

            fl_i = flavours[i]
            fl_f = flavours[f]

            if (fl_i != fl_f) and not ((fl_i == 21) or (fl_f == 21)):
                continue

            # Load heatmap
            heatmap_data = read_heatmap(fl_i, fl_f)

            # Plot the heatmap
            plot_heatmap(heatmap_data)
