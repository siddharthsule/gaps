#!/usr/bin/env python3

import os
import subprocess


def compile_code(runtype, args):
    """
    Compile the cpu-shower or gpu-shower with cmake and make.

    @param runtype: 'cpu' or 'gpu'
    """
    current_dir = os.getcwd()

    # Change to base directory for compilation
    os.chdir(args.base_dir)

    os.chdir(os.path.join('gaps', f'{runtype}-shower'))
    os.makedirs('build', exist_ok=True)
    os.chdir('build')

    subprocess.run(
        ['cmake', '..', f'-DGPROF={"ON" if args.gprof else "OFF"}'], check=True)
    subprocess.run(['make', '-j'], check=True)

    # Return to original directory
    os.chdir(current_dir)
