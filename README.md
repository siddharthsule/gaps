# GAPS: a GPU-Amplified Parton Shower

> **Version 1.1.0**: Some extra features added on demand, folders restructured.

Code for "An Algorithm to Parallelise Parton Showers on a GPU" [[arxiv:2403.08692](https://arxiv.org/abs/2403.08692)]

The aim of this project is to demonstrate how a Parton Shower Veto Algorithm can be written to run in parallel on a GPU. The code runs a simple LEP Event Generator on NVIDIA GPUs using CUDA. It is based on S. Höche's Tutorial on Parton Showers [[arxiv:1411.4085](https://arxiv.org/abs/1411.4085)].

## What can the code do on the GPU?

- Calculate the Matrix Element for $e^+ e^- \to q \bar{q}$ at 91.2 GeV
- Simulate a Final State Dipole Shower
- Calculate Jet Rates and Event Shapes

## Requirements

You will need an NVIDIA GPU, desgined for data centres (this code is verified to run on the NVIDIA Tesla V100 and A100 Devices). To build the code, you will need CMake, G++, Python and the NVIDIA Development Toolkit, which contains the NVCC compiler.

## Running the Code

The executable ```rungaps``` is written to simplify the use of the code. One can simply execute the command:

```bash
./rungaps
```

NB: If you get a permission denied error, please run ```chmod +x rungaps```.

This should build the program and generate 10000 events on the GPU. More customisation options are available, and are listed below:

```bash
# Simulate different number of events and build the code using multiple cpu cores
./rungaps -n nevents -c ncores

# Run C++ Simulation
./rungaps -n nevents -c ncores -r cpp

# Run the same number of events on C++ and CUDA and compare times
./rungaps -n nevents -c ncores -r compare

# Run a multitude of number of events 100 times, as seen in paper
./rungaps -c ncores -r full
```

The histograms are saved as yoda files [[arxiv:2312.15070](https://arxiv.org/abs/2312.15070)]. To generate the plots, use Rivet [[arxiv:1912.05451](https://arxiv.org/abs/1912.05451)] as follows:

```shell
rivet-mkhtml my-output.yoda:"Results" -s --mc-errs -c plots.conf
```

## Modifying Parameters

**New**: You can now adjust the Centre of Mass Energy using the ```-e``` flag

To focus on the computational aspects and make it simple to replicate the results in the paper, we don't allow direct access to the physics parameters (yet!). For now, please use the ```base.cuh``` file to adjust parameters like $\alpha_s(m_Z)$, $t_{C}$ and $n_{Bins}$.

***

### Sid Sule + Mike Seymour, March 2024

For issues and queries, email: [siddharth.sule@manchester.ac.uk](mailto:siddharth.sule@manchester.ac.uk)
