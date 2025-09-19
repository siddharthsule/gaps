# Calculating Observables and Histogramming

After the parton shower stage, we are left with multiple events of different sizes (meaning different numbers of partons in the final state). Like the shower, we write our code to treat each event independently and use the as-prescribed algorithm to calculate the observables.

The observables that are currently available for LEP are:

- **Durham Jet Rates**: These are the $y_{23}$, $y_{34}$, $y_{45}$ and $y_{56}$ Jet rates on a log scale. These tell us how many jets a particular event possibly has, along with intricate details like the transverse momentum spectra [[Phys. Lett. B 269 (1991)](https://inspirehep.net/literature/317695)]

- **Thrust**: This tells us how pencil-like or spherical a final state is [[Phys. Rev. Lett. 39 (1977)](https://link.aps.org/doi/10.1103/PhysRevLett.39.1587)]. There are two possible algorithms for calculating thrust: one that is more time-consuming but consistent for all event shapes and a simplified, faster algorithm that incorrectly calculates thrust for highly spherical events [[(Reported in Pythia)](https://pythia.org/latest-manual/EventAnalysis.html)]. We used the more time-consuming algorithm, first seen in [[Z. Physik C 1, 61 (1979)](https://link.springer.com/article/10.1007/BF01450381)]

- **Jet Masses and Broadenings**: These split the momenta of the partons with and against the thrust direction and yield properties of each hemisphere. [[Phys. Lett. B 272 (1991)](https://www.sciencedirect.com/science/article/pii/037026939191845M), [Phys. Lett. B. 295 (1992)](https://www.sciencedirect.com/science/article/pii/037026939291565Q)]

- **Dalitz Plot**: This studies the emitter and spectator momentum in a single emission event [[Phil. Mag. 44 (1953)](https://www.tandfonline.com/doi/abs/10.1080/14786441008520365)] (OFF by default though, as the writing time dominates the GPU analysis time!)

The observables that are currently available for LHC are:

- **Z Boson Observables**: In both LHC and LHCNLO, the momentum of the Z Boson is used to generate the distribution of its mass, rapidity, transverse momentum and azimuth. This is in line with the MC_ZINC analysis in Rivet [[arXiV:2404.15984](https://arxiv.org/abs/2404.15984), [MC_ZINC](https://rivet.hepforge.org/analyses/MC_ZINC.html)].

- **Leading and Subleading Jet Properties**: The final state partons were clustered using the anti-$k_T$ algorithm [[arXiv:0802.15984](https://arxiv.org/abs/0802.1189)]

One exciting aspect of this stage is binning. CUDA has a feature called atomic operations, where a variable can be operated on simultaneously. In our case, this means that once an observable is calculated for all events, every value of that observable is binned at once (Pretty cool, right!).

---

## Navigation

- [📚 Documentation Home](../README.md)
- [🏠 Repository Home](https://gitlab.com/siddharthsule/gaps)
