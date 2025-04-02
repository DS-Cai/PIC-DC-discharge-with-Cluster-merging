# PIC_DC_discharge

This is the code used in the paper "A Conservative-Constrained Clustering-Merging Algorithm for Particle-in-Cell Codes" 
Note:
Since the merging function calls a Python package, a Python runtime environment is required when using merging (If merging is not used, then this step is not necessary).
Follow the tutorial to install the k-means-constrained Python package: MATLAB Documentation. 
Then, modify the location where the Python package is called in the merging function, specifically in line 20: "py.k_means_constrained.KMeansConstrained".

<1> Case 1
Case 1 uses two groups of particles that follow a Maxwell distribution but have opposite velocities to test the effect of the particle merging algorithm on the particle distribution function (PDF). Both position and velocity spaces are 1D. The merging process is executed cell by cell. “test_minsize” represents the minimum cluster limits for Merge1-4, corresponding to the symbol "m" in the paper. 
Running Case1.m will get four comparison results for Merge1-4, using the same velocity distribution data before merging, as shown in the Figure of Case 1 merging results comparison.
 
<2> Case 2
Case 2 evaluates the quality of the PDF obtained from the merging algorithm using the two-stream instability model. Similarly, "test_minsize" represents the minimum cluster limits for Merge1-3. Running Case3.m will yield the particle distributions and energy change comparisons for four different merging types at various time points. The figure of "PDF for particles without merge" is an example of the results without merging.
<3> Case 3
Case 3 is a 1D3V simulation of direct current discharge between two parallel plate electrodes, considering electron neutrality, ion neutrality, and Coulomb collisions during the ionization process. The value of the variable “merge” is used to enable (1) or disable (0) merging, while “test_minsize” is used to set the minimum number of particles for merging, corresponding to the variable “m” in the paper. Due to the large number of results generated in this case, only the storage description of the main results is provided here.

Variables:
CPU_t: current time step and total CPU time.
“variable”_i represents ion, “variable”_e represents electron.
x: current particle position of ions. 
vx: current particle velocity along x-direction. 
vy: current particle velocity along y-direction.
vz: current particle velocity along z-direction.
weight: current macro particle weight.
pot_xt: variation of electric potential with time and coordinates.
ne_xt: variation of electron density with time and coordinates.
ni_xt: variation of ion density with time and coordinates.
meanei_xt: variation of electron mean energy with time and coordinates.
meanei_xt: variation of electron mean energy with time and coordinates.
Files:
conv.dat: storage time, Current macroscopic electron and ion number.
Picdata.dat: storage current result, including current time, electron and ion number, weight, position and velocity.

Ref:
https://github.com/joshlk/k-means-constrained

Bradley, P. S., K. P. Bennett, and Ayhan Demiriz. "Constrained k-means clustering." Microsoft Research, Redmond (2000): 1-8.

Google's SimpleMinCostFlow C++ implementation

A. Gonoskov, Agnostic conservative down-sampling for optimizing statistical representations and PIC simulations, Computer Physics Communications 271 (2022) 108200. https://doi.org/10.1016/j.cpc.2021.108200.

**Citations**

