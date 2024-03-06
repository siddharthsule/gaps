# GAPS: a GPU-Amplified Parton Shower

The aim of this project is to demonstrate how a Parton Shower Veto Algorithm can be written to run in parallel on a GPU. The code runs a simple LEP Event Generator on NVIDIA GPUs using CUDA. It is based on S. Höche's Tutorial on Parton Showers [[arxiv:1411.4085](https://arxiv.org/abs/1411.4085)].

## What can the code do on the GPU?

- Calculate the Matrix Element for $e^+ e^- \to q \bar{q}$
- Simulate a Dipole Shower
- Calculate Jet Rates and Event Shapes

## Running the Code

```bash
# Run CUDA/GPU Simulation
source runcode.sh cuda 10000 8 # 10000 events, 8 cores to compile

# Run C++ Simulation
source runcode.sh cpp 10000 8

# Run the same number of events and compare times
source runcode.sh compare 10000 8

# Run a multitude of number of events 10 times, as seen in paper
source runcode.sh full 8
```

The histograms are saved as yoda files [[arxiv:2312.15070](https://arxiv.org/abs/2312.15070)]. To generate the plots, use Rivet [[arxiv:1912.05451](https://arxiv.org/abs/1912.05451)] as follows:

```shell
rivet-mkhtml my-output.yoda:"Results" -s --mc-errs -c plots.conf
```

## Modifying Parameters

To make it simple to replicate the results in the paper, we don't allow direct access to the physics paramters (yet!). For now, please use the ```base.cuh``` file to adjust paramters like $\alpha_s(m_Z)$, $t_{C}$ and $n_{Bins}$.

***
### Mike Seymour + Sid Sule, February 2024
