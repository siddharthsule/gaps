# Event Generators on GPU

This code runs a simple LEP Event Generator on NVIDIA GPUs using CUDA.

The current progress is:
- ee -> qq Matrix Element: Fully on GPU
- Parton Shower: Prepared for GPU use, not implemented
- Hadronisation: No attempts made here
- Jet + Event Shapes: In Simple C++, CUDA not neccessary

This code is based on S. Hoeche's Tutorial on Parton Showers [arxiv:1411.4085]

## Building the Program

To build this program, you will need CMake and a valid CUDA installation (the code was tested on CUDA 11.7). To automate the building of the program, simply use the command

```shell
source build-project.sh # single core
source build-project.sh N # N cores