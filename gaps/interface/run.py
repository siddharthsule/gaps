#!/usr/bin/env python3

import os
import subprocess


def prepare_runparams(runtype, args):
    """
    @brief Prepare the arguements to be fed to the code.

    Param List:
    - NLO: 0 for no NLO, 1 for NLO
    - Root s: The center of mass energy in GeV
    - asmz: The strong coupling constant at the Z mass
    - t_c: The Shower cutoff in GeV
    - n_emissions_max: The maximum number of emissions
    - nevents: The number of events to generate
    - Event Number Offset: The offset for the event number (default 0)
    - Output filename: The name of the output file
    - Threads: The number of threads to use (only for GPU)

    @param runtype: 'cpu' or 'gpu'
    @param args: The command line arguments parsed by argparse
    @return: A list of parameters to be passed to the code
    """

    # Start with an empty list
    params = []

    # Add 0 if no NLO, 1 if NLO
    params.append('1' if args.nlo else '0')

    # Add the root s value
    params.append(str(args.root_s))

    # Add the strong coupling constant at the Z mass
    params.append(str(args.asmz))

    # Add the shower cutoff in GeV
    params.append(str(args.t_c))

    # Add the maximum number of emissions
    params.append(str(args.n_emissions_max))

    # Add the number of events
    params.append(str(args.nevents))

    # Add Event Number Offset - will be adjusted in run_cpu_cluster
    params.append('0')

    # Add output filename
    if runtype == 'cpu':
        params.append('cpu.yoda')
    elif runtype == 'gpu':
        params.append('gpu.yoda')

    # If GPU and Threads Given
    if runtype == 'gpu':
        params.append(str(args.threads))

    return params


def run(runtype, args):

    print(f'Running {runtype}-shower...')

    # Prepare the arguments based on the runtype - cpu or gpu
    params = prepare_runparams(runtype, args)

    # Define the run command
    run = [f'./gaps/{runtype}-shower/bin/{runtype}-shower']
    command = run + params

    # If Nsys Profiling
    if args.nsysprofile:
        profile = ['nsys', 'profile', '--stats=true']
        command = profile + command

    # Run the command
    subprocess.run(command)
