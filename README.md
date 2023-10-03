# BlaschkeSantalo
Matlab Codes for generating Blaschke Santalo diagrams related to the paper:  [The numerical approximation of Blaschke-Santalo diagrams using Centroidal Voronoi Tessellations](https://hal.science/hal-03966754) by B. Bogosel, G. Buttazzo and E. Oudet. If you use the code in your work, please give credit by citing the paper and this Github repository.

The library Geogram is needed to run the code. If possible, it should be compiled to use more than the default 6 digits of precision when writing to files. The file 'mesh_io.cpp' needs to be modified accordingly. The Matlab routine `fmincon` is used, therefore the Matlab Optimization toolbox is required. 

Make sure the correct path to the routine 'compute_RVD' is set in VoronoiGeogram.m

The main files of the toolbox are as follows:
- VoronoiGeogram.m  computes a restricted Voronoi diagram. A set of points and a polygon are given as inputs.
- BlaschkeSantaloWeb.m is used for the algebraic examples
- RunAPWweb.m is used for the geometric ones

Check the corresponding files for typical examples of usage. 
