# Matrix Element Calculation

This part of the code resembles the classic example of axpy (a*x + p = y). The matrix element is in charge of doing some calculations of pre-computed matrix elements and applying random numbers to create authentic events.

Here are the steps involved:

- The code generates a random number for the flavour of the quark-antiquark pair and for the angles theta and Phi.

- Then it does an analytical calculation, which might be complicated for us but is simple for the device as it is designed for this. This calculation gives the matrix element and the differential cross-section.

- It generates the momenta of the e+e- pair and the quark-antiquark pair and adds all of this information to the event variable

And that's all! The shower and observables sections are a lot more complicated. For more sophisticated calculations, look at the PEPPER [2311.06198](https://arxiv.org/abs/2311.06198.pdf) and MG4GPU [2303.18244](https://arxiv.org/abs/2303.18244) projects!
