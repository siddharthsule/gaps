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

    # Add Event Number Offset (use args.offset if present, else 0)
    params.append(str(getattr(args, 'offset', 0)))

    # Add output filename (use args.output if present, else default)
    if hasattr(args, 'output'):
        params.append(args.output)
    elif runtype == 'cpu':
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


def run_gpu_repeated(args, do_zip=True, zipname="gpu-yodas.zip"):
    """
    Mainly for NLL Tests, not part of the main GAPS workflow.
    """

    print(f'Running gpu-shower on repeat...')

    n_events_total = args.nevents
    batch_size = 1000000

    n_batches = n_events_total // batch_size
    remainder = n_events_total % batch_size

    offset = 0
    batch_outputs = []

    for i in range(n_batches):
        this_proc_events = batch_size
        output_filename = f"gpu-{i+1}.yoda"
        print(f"Running batch {i+1}/{n_batches + (1 if remainder > 0 else 0)}: ")

        args.nevents = this_proc_events
        args.offset = offset
        args.output = output_filename

        run('gpu', args)  # Use your run() function, which waits for completion
        batch_outputs.append(output_filename)
        offset += this_proc_events

    if remainder > 0:
        output_filename = f"gpu-{n_batches+1}.yoda"
        print(f"Running batch {n_batches+1}/{n_batches+1}: ")
        args.nevents = remainder
        args.offset = offset
        args.output = output_filename

        run('gpu', args)
        batch_outputs.append(output_filename)

    print("All GPU-based Monte Carlo simulations have completed.")

    if do_zip:
        print(f"Zipping GPU results into {zipname}...")
        os.makedirs(zipname.replace('.zip', ''), exist_ok=True)
        for f in batch_outputs:
            subprocess.run(f"mv {f} {zipname.replace('.zip', '')}/", shell=True, check=True)
        subprocess.run(["zip", "-r", zipname, zipname.replace('.zip', '')], check=True)
        subprocess.run(["rm", "-rf", zipname.replace('.zip', '')], check=True)
        print("GPU .yoda files have been zipped.\n")
