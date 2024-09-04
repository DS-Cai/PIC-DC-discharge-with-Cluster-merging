# PIC_DC_discharge
This is a 1D PIC code for DC discharge between parallel plates.

If merging is required, set "size_min" (greater than `constraint_num`) in the parameter_define, and ensure that a callable Python version of MATLAB is installed on your computer. Then, download and install the Constraint_k_means Python code from:  https://joshlk.github.io/k-means-constrained.

If merging is not required, remove the following two lines of code from do_one_cycle: [x_e,vx_e,vy_e,vz_e,WEIGHT_e]=merging(xN_e,xG_e,x_e,vx_e,vy_e,vz_e,WEIGHT_e,N_G,constraint_num,size_min) & [x_i,vx_i,vy_i,vz_i,WEIGHT_i]=merging(xN_i,xG_i,x_i,vx_i,vy_i,vz_i,WEIGHT_i,N_G,constraint_num,size_min)
RUn PIC_1D.
Ref:
https://github.com/joshlk/k-means-constrained
Bradley, P. S., K. P. Bennett, and Ayhan Demiriz. "Constrained k-means clustering." Microsoft Research, Redmond (2000): 1-8.
Google's SimpleMinCostFlow C++ implementation
A. Gonoskov, Agnostic conservative down-sampling for optimizing statistical representations and PIC simulations, Computer Physics Communications 271 (2022) 108200. https://doi.org/10.1016/j.cpc.2021.108200.
**Citations**
If you use this software in your research, please use the following citation:
