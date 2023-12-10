# Event Generators on GPU

This project runs a simple LEP Event Generator on NVIDIA GPUs using CUDA. It is based on S. Höche's Tutorial on Parton Showers ([arxiv:1411.4085](https://arxiv.org/abs/1411.4085)).

## Current Progress

- **ee -> qq Matrix Element**: Fully implemented on GPU.
- **Parton Shower**: Prepared for GPU use, not yet implemented.
- **Hadronisation**: No attempts made yet.
- **Jet + Event Shapes**: Implemented in simple C++, CUDA not necessary (for now).

## Building the Program

Ensure you have CMake and a valid CUDA installation (tested on CUDA 11.7). To build the program, use the following command:

```shell
source build-project.sh # single core
source build-project.sh N # N cores
```

## Running the Program

To run the event generator, use:

```shell
./bin/GPUEvGen # Default
./bin/GPUEvGen N # N Events
```
The program updates the user after every thousand events. At the end of the run, the output is stored as `output.yoda`.

To generate the plots, use:

```shell
rivet-mkhtml output.yoda:"Results" -s --mc-errs -c plots.conf
```
### Siddharth Sule, December 2023
