# Code Structure

> **Version 1.0.0**: All code in `include` and `src` to be split by component later

While the focus of the paper is the parton shower, there are several necessary components that are exciting in their own right. Below is a figure showing the structure of the classes.

**Note**: The code structure is identical for the C++ and CUDA generators. We want to demonstrate that one can get a speedup without altering the code much.

![Structure of the Files](structure.png)

Let's go through the Classes:

- Base: All the necessary definitions and settings are provided here. For now, you can adjust the Centre of Mass Energy and Cutoff and the physical contents here.

- Vec4 contains the definition of the four vectors—the primary tool for all our kinematics. This is one of many header-only files! This is because the code inside a .cu file might not be accessible to another .cu file (weird, I know; it could be me struggling to do things the C++ way, though!).

- Parton and Event: This is the backbone of the Generator. Unlike the original tutorial, we choose to store all and any properties of the code in these objects. The parton class holds the ID, momentum and colour of the parton. The event class stores an array of partons (it has to be static because CUDA won't allow dynamic arrays/unsuitable for speed). It also stores several parameters like the differential cross section and observable values and acts as a temporary store for shower parameters. It's good to have everything in there rather than in separate arrays!

- Histogram: This has a class for Bins and Histograms. Fun feature: CUDA allows you to do Atomic operations, i.e., I can bin all my events at the same time! It makes things unbelievably fast!

- Physics Classes
  - Matrix: Computes the matrix element for $ee \to qq$
  - Shower: Runs the final state shower
  - Observable: Calculates the observable and bins the Histogram
