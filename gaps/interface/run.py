#!/usr/bin/env python3

import os
import subprocess


def prepare_runparams(runtype, args):
    """
    @brief Prepare the arguements to be fed to the code.

    Param List:
    - Process: 1 for LEP, 2 for LHC
    - NLO: 0 for no NLO, 1 for NLO
    - Root s: The center of mass energy in GeV
    - asmz: The strong coupling constant at the Z mass
    - fixas: 1 if fixed asmz, 0 if running
    - noshower: 1 if skip shower, 0 if run shower
    - t_c: The Shower cutoff in GeV
    - n_emissions_max: The maximum number of emissions
    - nevents: The number of events to generate
    - Event Number Offset: The offset for the event number (default 0)
    - Output filename: The name of the output file
    - Partitioning: Do Event Record Partitioning (GPU only)
    - Threads: The number of threads to use (GPU only)

    @param runtype: 'cpu' or 'gpu'
    @param args: The command line arguments parsed by argparse
    @return: A list of parameters to be passed to the code
    """

    # Start with an empty list
    params = []

    # 0 - Add 1 for LEP or 2 for LHC
    params.append('1' if args.process.lower() == 'lep' else '2')

    # 1 - Add 0 if no NLO, 1 if NLO
    params.append('1' if args.nlo else '0')

    # 2 - Add the root s value
    params.append(str(args.root_s))

    # 3 - Add the strong coupling constant at the Z mass
    params.append(str(args.asmz))

    # 4 - If fixas is set, use the fixed asmz value
    params.append('1' if args.fixas else '0')

    # 5 - If noshower is set, skip the shower section
    params.append('1' if args.noshower else '0')

    # 6 - Add the shower cutoff in GeV
    params.append(str(args.t_c))

    # 7 - Add the maximum number of emissions
    params.append(str(args.n_emissions_max))

    # 8 - Add the number of events
    params.append(str(args.nevents))

    # 9 - Add Event Number Offset - will be adjusted in run_cpu_cluster
    params.append('0')

    # 10 - Add Storage file name
    if runtype == 'cpu':
        params.append('cpu.yoda')
    elif runtype == 'gpu':
        params.append('gpu.yoda')

    # 11 - If GPU and Threads Given
    if runtype == 'gpu':

        # Add partitioning flag
        params.append('1' if args.do_partitioning.lower() == 'yes' else '0')

        # Add the number of threads
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


def run_cpu_cluster(ncpu, args):

    print(f'Running cpu-shower on cpu cluster...')

    # Number of CPU processes to launch
    num_procs = ncpu
    n_events_total = args.nevents

    # Simple even split + handle remainder
    events_per_proc = n_events_total // num_procs
    remainder = n_events_total % num_procs

    processes = []
    offset = 0

    for i in range(num_procs):

        # This process gets either events_per_proc or +1 if remainder not exhausted
        this_proc_events = events_per_proc + (1 if i < remainder else 0)

        # Output file so each process writes to a different .yoda
        output_filename = f"cpu-{i+1}.yoda"

        print(
            f"[CPU {i}] Launching {this_proc_events} events (offset={offset}) → {output_filename}")

        # Get the code params as if we were running on a single CPU
        params = prepare_runparams('cpu', args)

        # adjust output filename and offset
        params[8] = str(this_proc_events)
        params[9] = str(offset)
        params[10] = output_filename

        # Path to the executable
        exe_path = ['./gaps/cpu-shower/bin/cpu-shower']

        # Combine the command with parameters
        command = exe_path + params

        # Spawn the process
        p = subprocess.Popen(command)
        processes.append(p)

        # Optional - store output from each process
        # p = subprocess.Popen(command, stdout=open(f"temp-{i+1}.dat", 'w'))
        # processes.append(p)

        # Update offset for the next process
        offset += this_proc_events

    # Wait for all processes
    for p in processes:
        ret_code = p.wait()
        if ret_code != 0:
            print(f"Process {p.pid} exited with code {ret_code}")

    print("All CPU-based runs have completed.")

    # Zipping is OFF
    # Don't want zipping time to be included in the profiling!
    # if not (args.runtype == 'full'):
    #     print("Zipping CPU results into cpu-yodas.zip...")
    #     os.makedirs('cpu-yodas', exist_ok=True)
    #     subprocess.run("mv cpu-*.yoda cpu-yodas/", shell=True, check=True)
    #     subprocess.run(["zip", "-r", "cpu-yodas.zip", "cpu-yodas"], check=True)
    #     subprocess.run(["rm", "-rf", "cpu-yodas"], check=True)
    #     print("CPU .yoda files have been zipped.\n")
