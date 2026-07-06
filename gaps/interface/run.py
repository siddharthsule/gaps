#!/usr/bin/env python3

import subprocess


def prepare_runparams(runtype, args):
    """
    @brief Prepare the arguments to be fed to the code.

    @param runtype: 'cpu' or 'gpu'
    @param args: The command line arguments parsed by argparse
    @return: A dictionary of parameters
    """

    # Build parameter dictionary with named keys
    params = {
        'process': '1' if args.process.lower() == 'lep' else '2',
        'nlo': '1' if args.nlo else '0',
        'root_s': str(args.root_s),
        'asmz': str(args.asmz),
        'fixas': '1' if args.fixas else '0',
        'use_cmw': '1' if args.use_cmw else '0',
        'noshower': '1' if args.noshower else '0',
        't_c': str(args.t_c),
        'n_emissions_max': str(args.n_emissions_max),
        'me2pdf': args.me2pdf,
        'showerpdf': args.showerpdf,
        'hadronise': '1' if args.hadronise else '0',
        'clmax': ','.join(str(x) for x in args.clmax),
        'clpow': ','.join(str(x) for x in args.clpow),
        'psplit': ','.join(str(x) for x in args.psplit),
        'pwt': ','.join(str(x) for x in args.pwt),
        'nevents': str(args.nevents),
        'id_offset': '0',
        'output_file': f'{runtype}.yoda',
    }

    # Add GPU-specific parameters
    if runtype == 'gpu':
        params['do_partitioning'] = '1' if args.do_partitioning.lower(
        ) == 'yes' else '0'
        params['threads'] = str(args.threads)

    return params


def params_to_list(params, runtype):
    """
    @brief Convert parameter dictionary to ordered list for C++/CUDA program.

    Order matches the expected argv order in main.cpp/main.cu

    @param params: Dictionary of parameters
    @param runtype: 'cpu' or 'gpu'
    @return: List of parameters in correct order
    """

    param_list = [
        params['process'],
        params['nlo'],
        params['root_s'],
        params['asmz'],
        params['fixas'],
        params['use_cmw'],
        params['noshower'],
        params['t_c'],
        params['n_emissions_max'],
        params['me2pdf'],
        params['showerpdf'],
        params['hadronise'],
        params['clmax'],
        params['clpow'],
        params['psplit'],
        params['pwt'],
        params['nevents'],
        params['id_offset'],
        params['output_file'],
    ]

    # Add GPU-specific parameters
    if runtype == 'gpu':
        param_list.extend([
            params['do_partitioning'],
            params['threads']
        ])

    return param_list


def run(runtype, args):

    print(f'Running {runtype}-shower...')

    # Prepare the arguments based on the runtype - cpu or gpu
    params = prepare_runparams(runtype, args)
    param_list = params_to_list(params, runtype)

    exe_path = f'{args.base_dir}/gaps/{runtype}-shower/bin/{runtype}-shower'

    command = [exe_path] + param_list

    # If Nsys Profiling
    if args.nsysprofile:
        profile = ['nsys', 'profile', '--trace=cuda,nvtx,osrt', '--stats=true']
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

        # Adjust output filename, event count, and offset
        params['nevents'] = str(this_proc_events)
        params['id_offset'] = str(offset)
        params['output_file'] = output_filename

        # Path to the executable using absolute path
        exe_path = f'{args.base_dir}/gaps/cpu-shower/bin/cpu-shower'

        # Combine the command with parameters
        command = [exe_path] + params_to_list(params, 'cpu')

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


def run_gpu_large_sample(args):

    print(f'Running gpu-shower with large sample in batches...')

    # GPU can handle max 1,000,000 events at once
    max_events_per_batch = 1000000
    n_events_total = args.nevents

    # Calculate number of batches needed
    num_batches = (n_events_total + max_events_per_batch -
                   1) // max_events_per_batch

    print(f' - Total events: {n_events_total}')
    print(f' - Events per batch: {max_events_per_batch}')
    print(f' - Number of batches: {num_batches}')
    print('')

    offset = 0

    for i in range(num_batches):

        # Calculate events for this batch
        events_remaining = n_events_total - offset
        this_batch_events = min(max_events_per_batch, events_remaining)

        # Output file for this batch
        output_filename = f"gpu-{i+1}.yoda"

        print(
            f"[Batch {i+1}/{num_batches}] Running {this_batch_events} events (offset={offset}) → {output_filename}")

        # Get the code params as if we were running a single GPU run
        params = prepare_runparams('gpu', args)

        # Adjust number of events, offset, and output filename
        params['nevents'] = str(this_batch_events)
        params['id_offset'] = str(offset)
        params['output_file'] = output_filename

        # Path to the executable using absolute path
        exe_path = f'{args.base_dir}/gaps/gpu-shower/bin/gpu-shower'

        # Combine the command with parameters
        command = [exe_path] + params_to_list(params, 'gpu')

        # If Nsys Profiling (apply to all batches)
        if args.nsysprofile:
            profile = ['nsys', 'profile',
                       '--trace=cuda,nvtx,osrt', '--stats=true']
            command = profile + command

        # Run the batch sequentially (GPU runs one batch at a time)
        ret_code = subprocess.run(command).returncode

        if ret_code != 0:
            print(f"Batch {i+1} exited with code {ret_code}")
            break

        # Update offset for the next batch
        offset += this_batch_events

    print("")
    print("All GPU batches have completed.")


