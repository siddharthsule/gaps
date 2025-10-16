# Tests and Analyses for GAPS Papers

This folder contains useful analyses and results for understanding and playing around with the GAPS code.

- **Herwig**: Input files, Rivet Analyses, Yoda Files and a Shell Script to run the simulation. Also includes the PDF Ratio diff to align Herwig and GAPS further. Used for Physics Cross Check in Paper 2.

- **SH-Tutorial**: Results from Stefan Hoche's "Introduction to Parton Showers and Matching" Tutorial [[1411.4085](https://arxiv.org/abs/1411.4085),[shoeche/tutorials](https://gitlab.com/shoeche/tutorials/)]. Useful for validating LEP results (the results should be identical). Used for Physics Cross Check in Paper 1.

- **plots.conf**: When plotting with Rivet, use this file using the argument `-c plots.conf` to add titles, axes labels and adjusted heights. For Rivet 4, you can adjust the Python file generated for each plot, too.

- **pdf-evaluations.cu**: A simple example of using LHAPDF on CPU and GPU. Compares the time required for 1,000,000 evaluations on both hardware platforms. Designed to test different flavours in each 'event'.

- **Python Scripts**: Used for the analyses in the papers, can be used on all hardware to generate the same(?) results. Useful to compare GPUs.

- **Results in Papers**: Although results can be generated from the Python Scripts, results from the papers are stored here, as raw data (the plots are available on ArXiV as TeX Source).
