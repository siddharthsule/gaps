import os
import subprocess
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Profiling Results to LaTeX')
parser.add_argument('--formatting-only', action='store_true',
                    help='Only format the pre-existing profile.dat')
args = parser.parse_args()

# ------------------------------------------------------------------------------
# Run the simulation

nev = 1000000
thr = 128
filename = 'profile.dat'

if not args.formatting_only:
    os.chdir('../')
    with open(filename, 'w') as f:
        subprocess.run(
            f"./rungaps -p LHC -nlo -t {thr} -n {nev} -nsys", shell=True, stdout=f)
    os.chdir('test')

# ------------------------------------------------------------------------------
# Parse CUDA Kernel Statistics

# Read File
with open('../' + filename, 'r') as file:
    content = file.read()

# Find the CUDA Kernel Statistics section (Edit here if issues!)
start_markers = ["CUDA Kernel Statistics:", "[6/8] Executing"]
end_marker = "[7/8] Executing"

# Try to find any of the start markers
start_idx = -1
for marker in start_markers:
    start_idx = content.find(marker)
    if start_idx != -1:
        break

if start_idx == -1:
    raise ValueError(
        f"Could not find any of the start markers: {start_markers}")

end_idx = content.find(end_marker)

if end_idx == -1:
    raise ValueError(f"Could not find end marker: {end_marker}")

# Extract the section
cuda_section = content[start_idx:end_idx]

# Split into lines and find the data lines
lines = cuda_section.split('\n')

# Find the header line (contains "Time (%)")
header_idx = None
data_start_idx = None

for i, line in enumerate(lines):
    if "Time (%)" in line and "Total Time (ns)" in line:
        header_idx = i
        # Data starts after the separator line (next line with dashes)
        for j in range(i + 1, len(lines)):
            if "--------" in lines[j]:
                data_start_idx = j + 1
                break
        break

if header_idx is None or data_start_idx is None:
    raise ValueError("Could not find data table structure")

# Extract column headers
columns = ['Time_Percent', 'Total_Time_ns', 'Instances',
           'Avg_ns', 'Med_ns', 'Min_ns', 'Max_ns', 'StdDev_ns', 'Name']

# Parse data lines
data_rows = []

# Parse each data line
for line in lines[data_start_idx:]:
    line = line.strip()
    if not line or line.startswith('['):
        break

    # Split the line by whitespace, but be careful with the name field
    # The name field can contain spaces and special characters
    parts = line.split()

    # Ensure valid Line
    if len(parts) < 9:
        continue

    try:
        # Extract numeric fields (first 8 columns)
        time_percent = float(parts[0])
        total_time = int(parts[1].replace(',', ''))
        instances = int(parts[2].replace(',', ''))
        avg_ns = float(parts[3].replace(',', ''))
        med_ns = float(parts[4].replace(',', ''))
        min_ns = int(parts[5].replace(',', ''))
        max_ns = int(parts[6].replace(',', ''))
        stddev_ns = float(parts[7].replace(',', ''))

        # The name is everything after the 8th field
        name_parts = parts[8:]
        name = ' '.join(name_parts)

        data_rows.append([
            time_percent, total_time, instances, avg_ns, med_ns,
            min_ns, max_ns, stddev_ns, name
        ])

    except (ValueError, IndexError) as e:
        print(f"Skipping line due to parsing error: {line}")
        continue

# Create DataFrame
df = pd.DataFrame(data_rows, columns=columns)

# Only keep name, instances, time, time percentage
df = df[['Name', 'Instances', 'Total_Time_ns', 'Time_Percent']]

# ------------------------------------------------------------------------------
# Replace long kernel names with shorter versions

# Dictionary to switch long names with short names
name_mapping = {
    'select_winner_split_func': 'Generate the Trial Emission',
    'veto_alg': 'Vetoing Process',
    'xfxQ2': 'PDF 1',
    'setup_pdfratio': 'Setup PDF Ratio Calculation',
    'do_splitting': 'Perform the Emission',
    'check_cutoff': 'Checking Shower Cutoff',
    'check_too_many_particles': 'Check for Max Particles Per Event',
    'select_flavour': 'PDF 2',
    'cluster_genkt': 'Anti-$k_T$ Jet Clustering',
    'fill_histos': 'Fill Histograms',
    'validate_events': 'Validate Events',
    'lhc_lo': 'LHC LO',
    'lhc_nlo': 'LHC NLO',
    'prep_shower': 'Prepare Shower',
    'calculate_mczinc': 'Calculate Z Observables',
    'lo_event': 'NLO 1',
    'h_event': 'NLO 2',
    'c_terms': 'NLO 3',
    'bvic_terms': 'NLO 4'
}

# Replace long names with short names
for old_name, new_name in name_mapping.items():
    df.loc[df['Name'].str.contains(
        old_name, case=False, na=False), 'Name'] = new_name

# There are many __parallel_for but the one which has only one Instance
# Call this Device Prep
df.loc[df['Name'].str.contains('__parallel_for', case=False, na=False) & (
    df['Instances'] == 1), 'Name'] = 'Device Prep'

# There is one _parallel_for Device Prep Task for PDF Evaluation
# it should roughly be twice the instances of "Vetoing Process"
runs = df.loc[df['Name'] == "Vetoing Process"]["Instances"].values[0]
below = (runs * 2) - 50
above = (runs * 2) + 50
df.loc[df['Name'].str.contains('__parallel_for', case=False, na=False) & (
    df['Instances'] > below) & (df['Instances'] < above), 'Name'] = 'PDF 3'

# ------------------------------------------------------------------------------
# Combine everything with 16 instances into Event Record Partitioning

# Get entries with 16 instances
mini_df = df.loc[df['Instances'] == 16]

# Create combined row
total_row = mini_df.sum(numeric_only=True)
total_row['Name'] = 'Event Record Partitioning'
total_row['Instances'] = 16

# Remove original rows and add combined row
df = df[df['Instances'] != 16]
total_row_df = pd.DataFrame([total_row])
df = pd.concat([df, total_row_df], ignore_index=True)

# ------------------------------------------------------------------------------
# Combine PDF 1, 2 and 3 rows into a single row

# Get PDF entries
pdf1_df = df.loc[df['Name'] == 'PDF 1']
pdf2_df = df.loc[df['Name'] == 'PDF 2']
pdf3_df = df.loc[df['Name'] == 'PDF 3']

# Combine the rows
combined_row = pdf1_df.sum(
    numeric_only=True) + pdf2_df.sum(numeric_only=True) + pdf3_df.sum(numeric_only=True)
combined_row['Name'] = 'PDF (11 Evaluations + Selection)'
combined_row['Instances'] = pdf2_df['Instances'].values[0]

# Remove original rows and add combined row
df = df[~df['Name'].isin(['PDF 1', 'PDF 2', 'PDF 3'])]
df = pd.concat([df, pd.DataFrame([combined_row])], ignore_index=True)

# ------------------------------------------------------------------------------
# Combine NLO 1, 2, 3, and 4 rows into a single row

nlo1_df = df.loc[df['Name'] == 'NLO 1']
nlo2_df = df.loc[df['Name'] == 'NLO 2']
nlo3_df = df.loc[df['Name'] == 'NLO 3']
nlo4_df = df.loc[df['Name'] == 'NLO 4']

# Combine the rows
combined_row = nlo1_df.sum(numeric_only=True) + nlo2_df.sum(numeric_only=True) + \
    nlo3_df.sum(numeric_only=True) + nlo4_df.sum(numeric_only=True)
combined_row['Name'] = 'NLO Event Generation'
combined_row['Instances'] = nlo1_df['Instances'].values[0]

# Remove original rows and add combined row
df = df[~df['Name'].isin(['NLO 1', 'NLO 2', 'NLO 3', 'NLO 4'])]
df = pd.concat([df, pd.DataFrame([combined_row])], ignore_index=True)

# ------------------------------------------------------------------------------
# Generate LaTeX table

print(df.sort_values(by='Time_Percent', ascending=False))

# Filter and sort kernels by time percentage
top_kernels = df[df['Time_Percent'] >= 0.1].sort_values(
    by='Time_Percent', ascending=False)

# Format the data
top_kernels['Time_Percent'] = top_kernels['Time_Percent'].round(2)
top_kernels['Short_Name'] = top_kernels['Name'].apply(
    lambda x: x[:60] + '...' if len(x) > 60 else x
)

# Select columns for the table
table_data = top_kernels[['Short_Name',
                          'Instances', 'Total_Time_ns', 'Time_Percent']]

# Generate LaTeX table
latex_table = """\\begin{table}[h!]
\\centering
\\label{tab:cuda_kernels}
\\begin{tabular}{l|r|r|r}
\\hline
\\textbf{Name} & \\textbf{Instances} & \\textbf{Total Time (ns)} & \\textbf{Time (\\%)} \\\\
\\hline
"""

for _, row in table_data.iterrows():
    # Apply italics to specific names
    italic_names = ['Event Record Partitioning', 'Device Prep']
    kernel_name = f"\\textit{{{row['Short_Name']}}}" if row['Short_Name'] in italic_names else row['Short_Name']
    latex_table += f"{kernel_name} & {row['Instances']} & {row['Total_Time_ns']:,} & {row['Time_Percent']:.1f} \\\\\n"

latex_table += """\\hline\n"""
latex_table += """\\end{tabular}
\\end{table}
"""

# Print to Terminal
print("\nLaTeX table preview:")
print(latex_table)

# Write to file
with open("cuda-kernels-table.dat", "w") as f:
    f.write(latex_table)
