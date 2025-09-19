#!/usr/bin/env python3

import os
import subprocess


def compile_code(runtype, args):
    """
    @brief Compile the code with cmake and make.

    Common compile flow:
      cd gaps/[runtype]-shower
      mkdir -p build
      cd build
      cmake ..
      make -j8

    @param runtype: 'cpu' or 'gpu'
    """

    # Save the current directory
    home_dir = os.getcwd()

    # Change to the directory where the code is
    folder = f"{runtype}-shower"  # 'cpu-shower' or 'gpu-shower'
    dir_path = os.path.join('gaps', folder)
    os.chdir(dir_path)

    # Compile the code
    os.makedirs('build', exist_ok=True)
    os.chdir('build')

    # CMAKE
    if args.gprof:
        # If gprof profiling is enabled, set the CMake variable
        subprocess.run(['cmake', '..', '-DGPROF=ON'], check=True)
    else:
        subprocess.run(['cmake', '..', '-DGPROF=OFF'], check=True)

    # MAKE
    subprocess.run(['make', '-j'], check=True)

    # Return to the original directory
    os.chdir(home_dir)
