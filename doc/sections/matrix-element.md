# Matrix Element Calculation

This part of the code is akin to the classic example of axpy (a*x + p = y). The matrix element is incharge of doing some calculation of pre-computed matrix elements and apply random numbers to create authentic events.

Here are the steps involved:

- The code generates a random number for the flavour of the quark-antiquark pair and for the angles theta and Phi.

- Then it does an analytical calulation, which might be complicated for us, but is simple for the device as it is designed for this. This calculation gives the matrix element, and hence the diff. cross section.

- It generates the momenta of the e+e- pair and the quark antiquark pair, and adds all of this information to the event variable

And that's all! The shower and observables sections are a lot more complicated. For more sophisticated calculations have a look at PEPPER and MG4GPU projects!
