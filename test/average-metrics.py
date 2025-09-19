import os
import argparse
import subprocess
import numpy as np
import pandas as pd
from scipy.stats import iqr

parser = argparse.ArgumentParser(description='Tune and Partitioning Tests')
parser.add_argument('--nreps', type=int, default=10,
                    help='Number of repetitions')
args = parser.parse_args()

nev = [10000, 100000, 1000000]
do_partitioning = 'no'

# ------------------------------------------------------------------------------
# Run Mutliple Times

os.chdir('../')
if os.path.exists('gpu-time.dat'):
    os.remove('gpu-time.dat')

for n in nev:
    for i in range(args.nreps):
        subprocess.run(
            f"./rungaps -p LHC -nlo -do_partitioning {do_partitioning} -n {n}", shell=True)

os.chdir('test')

# ------------------------------------------------------------------------------
# Calculate Median and IQR

data = pd.read_csv('../gpu-time.dat', header=None, delimiter=',')
data.columns = ["NLO", "Shower", "Observables", "Total"]

# the dataset is len(nev) sets of nreps repetitions. We want the median and IQR for each set.
median = np.zeros((len(data) // args.nreps, data.shape[1]))
iqr_values = np.zeros((len(data) // args.nreps, data.shape[1]))

# Calculate the Medians and IQRs
for i in range(len(data) // args.nreps):
    start = i * args.nreps
    end = start + args.nreps

    # Get subset of data
    subset = data[start:end]

    median[i] = np.median(subset, axis=0)
    iqr_values[i] = iqr(subset, axis=0)

# Print Neatly in Transposed Table
# Create row headers for the transposed table (nev values)
nev_headers = [f"N_EV={nev[i]}" for i in range(len(median))]

# Generate LaTeX table string
latex_table = ""

# Print header row with nev values
header_row = "Metric".ljust(12) + " & " + " & ".join(header.ljust(15)
                                                     for header in nev_headers) + " \\\\\n"
latex_table += header_row

# Print each metric as a row
for j in range(median.shape[1]):
    metric_name = data.columns[j]

    # Add horizontal line before Total row
    if metric_name == "Total":
        latex_table += "\\hline\n"

    # Make Total row bold
    if metric_name == "Total":
        row_values = " & ".join(f"\\textbf{{{median[i, j]:.2f}}}".ljust(
            15) for i in range(len(median)))
        latex_table += f"\\textbf{{{metric_name.ljust(12)}}} & {row_values} \\\\\n"
    else:
        row_values = " & ".join(f"{median[i, j]:.2f}".ljust(
            15) for i in range(len(median)))
        latex_table += f"{metric_name.ljust(12)} & {row_values} \\\\\n"

# Print the table
print(latex_table)

# Write to file
with open('average-metrics-table.dat', 'w') as f:
    f.write(latex_table)
