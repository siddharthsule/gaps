#!/usr/bin/env python3

import os
import subprocess


def install_lhapdf(base_dir):
    """
    @brief Install LHAPDF.

    Install LHAPDF from the source code.
    @param base_dir: Base directory where rungaps is located
    """

    # Boolean to check if LHAPDF is installed
    lhapdf_installed = False

    # Current Directory
    current_dir = os.getcwd()

    # install Directory
    lhapdf_dir = os.path.join(base_dir, 'gaps/lhapdf-gpu/')

    # Create the install directory
    if not os.path.exists(lhapdf_dir + 'install'):
        os.makedirs(lhapdf_dir + 'install')

    # Change to the install directory
    os.chdir(lhapdf_dir + 'install')

    # Set this directory as $INSTALL_DIR
    install_dir = os.getcwd()
    print(f"INSTALL_DIR is set to: {install_dir}")

    # Make an src folder to avoid confusion
    if not os.path.exists('src'):
        os.makedirs('src')

    # Change to the src directory
    os.chdir('src')

    # Get GPU LHAPDF - From GITLAB
    subprocess.run(['git', 'clone', '-b', 'kokkos_version',
                   'https://gitlab.com/hepcedar/lhapdf.git'])
    if not os.path.exists('lhapdf'):
        print("------------------------------------------------")
        print("Error: Failed to clone LHAPDF repository!")
        print("Sometimes running ./rungaps.sh again can fix this issue.")
        print("If the issue persists, please manually clone the repository:")
        print("\t git clone -b kokkos_version https://gitlab.com/hepcedar/lhapdf.git")
        print("------------------------------------------------")
        return False
    os.chdir('lhapdf')

    # Checkout specific commit for stability
    subprocess.run(
        ['git', 'checkout', 'ced1118b3a288117c68631df7bad784aae2a6f25'])

    # Reconf and Configure
    subprocess.run(['autoreconf', '-vi'])
    subprocess.run(
        ['./configure', f'--prefix={install_dir}', '--disable-python'])

    subprocess.run(['make', '-j'])
    subprocess.run(['make', 'install'])

    # Change to the install directory
    os.chdir(install_dir)

    # Place PDF tar.gz files in share/LHAPDF and extract them
    os.chdir('share/LHAPDF')
    pdfs = ['CT14lo', 'NNPDF40MC_lo_as_01180', 'NNPDF40MC_nlo_as_01180']
    for pdf in pdfs:
        subprocess.run(['cp', f'../../../{pdf}.tar.gz', '.'])
        subprocess.run(['tar', '-xzf', f'{pdf}.tar.gz'])
    os.chdir(install_dir)

    # Create a path.txt file
    os.chdir("..")
    with open('path.txt', 'w') as file:
        file.write(install_dir)

    # Return to the current directory
    os.chdir(current_dir)

    # Check if LHAPDF is installed
    if os.path.exists(lhapdf_dir + 'path.txt'):
        with open(lhapdf_dir + 'path.txt', 'r') as file:
            lhapdf_path = file.read().strip()

        if os.path.exists(lhapdf_path):
            lhapdf_installed = True
            print("------------------------------------------------")
            print("LHAPDF (kokkos_version) installed successfully!")
            print("------------------------------------------------")

    return lhapdf_installed


def check_lhapdf_path(base_dir):
    """
    @brief Check if the LHAPDF path is correct.

    If the file lhapdf_path.txt exists, check if the path inside is correct.
    If the path is correct, return it. If not, install LHAPDF within GAPS.
    @param base_dir: Base directory where rungaps is located
    """

    lhapdf_sorted = False

    while not lhapdf_sorted:
        # Check if the file lhapdf_path.txt exists
        lhapdf_txt = os.path.join(base_dir, 'gaps/lhapdf-gpu/path.txt')

        # Make a file lhapdf_txt if it does not exist
        if not os.path.exists(lhapdf_txt):
            with open(lhapdf_txt, 'w') as file:
                file.write('')

        # Check if the path inside is correct
        with open(lhapdf_txt, 'r') as file:
            lhapdf_path = file.read().strip()

        if os.path.exists(lhapdf_path):
            print(f"Using LHAPDF (kokkos_version) from {lhapdf_path}")
            return lhapdf_path

        # If path is not correct, install LHAPDF within GAPS
        print("------------------------------------------------")
        print("GAPS uses LHAPDF (kokkos_version) for GPU PDFs.")
        print(f"It will be installed in {os.path.join(base_dir, 'gaps/lhapdf-gpu')}.")
        print("")
        print("To manually state LHAPDF path, Ctrl+C and do:")
        print(f"\t echo '[PATH-2-LHAPDF]' > {lhapdf_txt}")
        print("------------------------------------------------")

        # Wait for 5 seconds before installing LHAPDF, to allow user to cancel
        subprocess.run(['sleep', '5'])

        # Install LHAPDF
        installed_lhapdf = install_lhapdf(base_dir)

        # Install issue: break the loop
        if not installed_lhapdf:
            break

        # Read the path again
        with open(lhapdf_txt, 'r') as file:
            lhapdf_path = file.read().strip()

        # Verify the path
        if os.path.exists(lhapdf_path):
            lhapdf_sorted = True
            print("--------------")
            print("LHAPDF sorted!")
            print("--------------")
        else:
            print("--------------------------------------")
            print("Error: The LHAPDF path does not exist!")
            print("--------------------------------------")

    return lhapdf_path
