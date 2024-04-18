# Getting Started

The code simulates just one experiment, so this should take a little time.

The file `gaps.sh` can be used to build and operate both the C++ and CUDA generators. It has been coded with all the routines (including the one for paper results).

To run GAPS, you will need the following:

- An NVIDIA V100, A100 or above: These are the only GPUs with the required features
- CMake: To create the makefile
- NVCC and GCC: to build the two generators
- Python: To make plots of the results

Simply execute the command:

```bash
./gaps.sh
```

NB: If you get a permission denied error, please run ```chmod +x gaps.sh```.

This should build the program and generate 10000 events on the GPU. The output should look something like this:

```bash
-------------------------------------------------
| GAPS: a GPU-Amplified Parton Shower |
-------------------------------------------------
Process: e+ e- --> q qbar
Number of Events: 10000

Initialising...
Generating Matrix Elements...
Showering Partons...
Analysing Events...

EVENT GENERATION COMPLETE

ME Time: 0.000666688 s
Sh Time: 0.0235208 s
An Time: 0.00896093 s

Total Time: 0.0331484 s

Histograms written to gaps.yoda
Timing data written to gaps-time.dat
------------------------------------------------
```

Then you have free reign over what you need. Like [README.md](../../README.md), here are all the possible ways you can run gaps:

```bash
# Simulate different numbers of events and build the code using multiple CPU cores
./gaps.sh -n nevents -c ncores

# Run C++ Simulation
./gaps.sh -n nevents -c ncores -r cpp

# Run the same number of events on C++ and CUDA and compare times
./gaps.sh -n nevents -c ncores -r compare

# Run a multitude of number of events 100 times, as seen in the paper
./gaps.sh -c ncores -r full
```

And that's all there is to it! In upcoming versions, we'll add features like increasing the number of events, different centres of mass energies, new analyses, and potentially some matching!
