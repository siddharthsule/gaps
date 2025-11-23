# Matrix Element Calculation

This part of the code resembles the classic example of axpy (a*x + p = y). The matrix element is in charge of doing some calculations of pre-computed matrix elements and applying random numbers to create authentic events.

Here are the steps involved:

- The code generates the momentum and layout of the initial and final state particles.

- Then it does an analytical calculation, which might be complicated for us, but is simple for the device, as it is designed for this. This calculation gives the matrix element and the differential cross-section.

The shower and observables sections are a lot more complicated. For more sophisticated calculations, look at the PEPPER [[2311.06198](https://arxiv.org/abs/2311.06198.pdf)] and MG4GPU [[2507.21039](https://arxiv.org/abs/2507.21039)] projects!

---

## Navigation

- [📚 Documentation Home](../README.md)
- [🏠 Repository Home](https://gitlab.com/siddharthsule/gaps)
