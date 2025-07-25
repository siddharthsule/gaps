import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# from scipy.stats import iqr

# Load the data
cpu0 = pd.read_csv('cpu-time-0.dat')
cud0 = pd.read_csv('gpu-time-0.dat')
cpu1 = pd.read_csv('cpu-time-1.dat')
cud1 = pd.read_csv('gpu-time-1.dat')

# Name the Columns Matrix, Shower, Observables, and Total
cpu0.columns = ['Matrix', 'Shower', 'Observables', 'Total']
cud0.columns = ['Matrix', 'Shower', 'Observables', 'Total']
cpu1.columns = ['Matrix', 'Shower', 'Observables', 'Total']
cud1.columns = ['Matrix', 'Shower', 'Observables', 'Total']

# Get the median speedup between cpu and cud for the two versions
speedup0 = cpu0.median() / cud0.median()
speedup1 = cpu1.median() / cud1.median()

# Get the median speedup between the two versions for cpu and cud
speedup_cpu = cpu0.median() / cpu1.median() * 100
speedup_cud = cud0.median() / cud1.median() * 100

# Concatenate all results
speedup0 = pd.concat([speedup0, speedup1], axis=1)
speedup1 = pd.concat([speedup_cpu, speedup_cud], axis=1)
speedup = pd.concat([speedup0, speedup1], axis=1)

# Name the columns
speedup.columns = ['Old Speedup', 'New Speedup',
                   'CPU Old/New (%)', 'GPU Old/New (%)']

# Print the results to 3 decimal places
print(speedup.round(3))
