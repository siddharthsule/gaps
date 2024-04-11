#------------------------------------------------------------------------------
#!/usr/bin/env python
# ------------------------------------------------------------------------------

# GAPS - Run Script
# -----------------
# This script is used to compile and run the GAPS and C++ Shower codes. It
# provides a number of options to control the number of events, the number of
# cores to use, and the type of run to perform. The run types are:
#   - gaps: Run the GAPS simulation
#   - cpp: Run the C++ Shower simulation
#   - compare: Run both the GAPS and C++ Shower and compare the results
#   - full: Run both the GAPS and C++ Shower for a range of event numbers

# ------------------------------------------------------------------------------

import argparse
import os
import shutil
import subprocess

# Set up argument parser
parser = argparse.ArgumentParser(description='Run GAPS or C++ Shower')
parser.add_argument('-n', type=int, default=10000, help='set the number of events (default: 10000)')
parser.add_argument('-e', type=float, default=91.2, help='set the CoM energy of the system (default: 91.2)')
parser.add_argument('-c', type=int, default=1, help='set the number of cores (default: 1)')
parser.add_argument('-r', type=str, default='gaps', help='set the run type (default: gaps, options: gaps, cpp, compare, full)')

args = parser.parse_args()

# Compile code
def compile(dir):
    print(f'Compiling {dir}')
    os.chdir(dir)
    os.makedirs('build', exist_ok=True)
    os.chdir('build')
    subprocess.run(['cmake', '..'])
    subprocess.run(['make', '-j', str(args.c)])

# Run GAPS or C++ Shower
def run(runtype, events, energy):
    print(f'Running {runtype}')
    subprocess.run([f'./{runtype}/bin/{runtype}', str(events), str(energy)])

# Compile and run based on runtype
if args.r in ['gaps', 'compare']:
    compile('gaps')
    run('gaps', args.n, args.e)

if args.r in ['cpp', 'compare']:
    compile('cpp-shower')
    run('cpp-shower', args.n, args.e)


if args.r == 'full':
    # Remove previous results and make folder
    if os.path.exists('results'):
        shutil.rmtree('results')
    os.makedirs('results', exist_ok=True)

    # Clear previous log files
    if os.path.exists('cpp-time.dat'):
        os.remove('cpp-time.dat')
    if os.path.exists('gaps-time.dat'):
        os.remove('gaps-time.dat')

    # Run the comparison 100 times, for different number of events
    neventsarray = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000]
    for n in neventsarray:
        # Run and store the output in a log file
        for i in range(1, 101):
            print(f"Running GAPS with {n} events")
            subprocess.run(['./gaps/bin/gaps', str(n), str(args.e)])
            print(f"Running C++ Shower with {n} events")
            subprocess.run(['./cpp-shower/bin/cpp-shower', str(n), str(args.e)])

    # Move the log files to the results directory
    shutil.move('cpp-time.dat', 'results/')
    shutil.move('gaps-time.dat', 'results/')
    shutil.move('cpp.yoda', 'results/')
    shutil.move('gaps.yoda', 'results/')