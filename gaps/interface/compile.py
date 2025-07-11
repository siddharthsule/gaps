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

    # Print the current working directory
    print(f"Compiling {folder}")

    # Compile the code
    os.makedirs('build', exist_ok=True)
    os.chdir('build')

    is_compiled = False
    try:
        with open('build_output.log', 'a') as log_file:
            # CMAKE
            if args.gprof:
                # If gprof profiling is enabled, set the CMake variable
                subprocess.run(['cmake', '..', '-DGPROF=ON'], check=True,
                               stdout=log_file, stderr=log_file)
            else:
                subprocess.run(['cmake', '..', '-DGPROF=OFF'], check=True,
                               stdout=log_file, stderr=log_file)

            # MAKE
            subprocess.run(['make', '-j'], check=True,
                           stdout=log_file, stderr=log_file)

        # IF Compilation is successful
        is_compiled = True

    except subprocess.CalledProcessError:
        red = "\033[31m"
        res = "\033[0m"
        print(f"{red}Build failed! Check {os.path.abspath('build_output.log')}{res}")
        is_compiled = False

    # Return to the original directory
    finally:
        os.chdir(home_dir)

    return is_compiled
