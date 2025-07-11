# GAPS: a GPU-Amplified Parton Shower

> **Version 1.3.0**: Improved code structure, more parameter options, MC@NLO Matching

Codebase for [M. H. Seymour and S. Sule, _An Algorithm to Parallelise Parton Showers on a GPU_, SciPost Phys. Codebases 33 (2024)](https://scipost.org/SciPostPhysCodeb.33)

The aim of this project is to demonstrate how a Parton Shower Veto Algorithm can be written to run in parallel on a GPU. The code runs a simple LEP Event Generator on NVIDIA GPUs using CUDA. It is based on S. Höche's Tutorial on Parton Showers [[arxiv:1411.4085](https://arxiv.org/abs/1411.4085)].

## What can the code do on the GPU?

- Generate the process $e^+ e^- \to q \bar{q}$ at 91.2 GeV
- Convert to NLO using Catani-Seymour Subtraction
- Generate Emissions using a Dipole Shower
- Calculate Jet Rates and Event Shapes

## Requirements

You will need an NVIDIA GPU, with CUDA Compatibility 7.0 [[Guide](https://developer.nvidia.com/cuda-gpus)]. To build the code, you will need `cmake`, `g++`, `python` and NVIDIA's `nvcc` compiler.

## Running the Code

The executable ```rungaps``` is written to simplify the use of the code. One can simply execute the command:

```bash
./rungaps
```

NB: If you get a permission denied error, please run ```chmod +x rungaps```.

This should build the program and generate 10000 events on the GPU. More customisation options are available, and are listed below:

```bash
# Simulate N Events on GPU (Default)
./rungaps -n nevents

# Simulate N Events on CPU
./rungaps -n nevents -c ncores -r cpu

# Simulate N Events on CPU and GPU and compare times and results
./rungaps -n nevents -c ncores -r compare

# Simulate a range of N on both CPU and GPU and compare (in paper)
./rungaps -c ncores -r full
```

The histograms are saved as `yoda` files [[arxiv:2312.15070](https://arxiv.org/abs/2312.15070)]. To generate the plots, use `rivet` [[arxiv:1912.05451](https://arxiv.org/abs/1912.05451)] as follows:

```shell
rivet-mkhtml --mc-errs -c test/plots.conf <yoda file>
```

## Going Further

**New!** You can now adjust the following paramters:

- `-nlo`: Generate NLO Events
- `-e, --root_s`: Adjust the center of mass energy
- `-asmz, -t_c`,: Adjust $\alpha_s(m_Z)$ and the shower cutoff $t_{C}$
- `-n_em_max`: Limit the number of emissions
- `-nsys, -codecarbon, -gprof`: Profiling Tools
- `-t`: Number of threads per block on the GPU

To learn more about the code and how it all works, see the [documentation](doc/README.md).

***

### Sid Sule + Mike Seymour, July 2025

For issues and queries, email: [siddharth.sule@manchester.ac.uk](mailto:siddharth.sule@manchester.ac.uk)
