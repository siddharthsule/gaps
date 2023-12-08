# Event Generators on GPU

This code runs a simple LEP Event Generator on NVIDIA GPUs using CUDA.

The current progress is:
- ee -> qq Matrix Element: Fully on GPU
- Parton Shower: Prepared for GPU use, not implemented
- Hadronisation: No attempts made here
- Jet + Event Shapes: In Simple C++, CUDA not neccessary

This code is based on S. Hoeche's Tutorial on Parton Showers [ArXiV:1411.4085]