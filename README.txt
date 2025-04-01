# PIC_DC_discharge

This is the code used in the paper "A Conservative-Constrained Clustering-Merging Algorithm for Particle-in-Cell Codes" 
Note:
Since the merging function calls a Python package, a Python runtime environment is required when using merging (If merging is not used, then this step is not necessary).
Follow the tutorial to install the k-means-constrained Python package: MATLAB Documentation. 
Then, modify the location where the Python package is called in the merging function, specifically in line 20: "py.k_means_constrained.KMeansConstrained".

<1> Case 1
In Case1.m, test_minsize stores the minimum cluster limits for Merge1-4, corresponding to the symbol "m" in the paper. 
Running Case1.m will get four comparison results for Merge1-4, using the same velocity distribution data before merging.
<2> Case 2




This is a 1D3V PIC code for DC discharge between parallel plates.

If merging is required, set "size_min" (greater than 'constraint_num') in the parameter_define, and ensure that a callable Python version of MATLAB is installed on your computer. Then, download and install the Constraint_k_means Python code from:  https://joshlk.github.io/k-means-constrained. For calling Python code from MATLAB, please refer to：https://ww2.mathworks.cn/help/matlab/matlab_external/ways-to-call-python-from-matlab.html#mw_53c718f8-1bbb-49bd-a17a-344979a5878a

If merging is not required, remove the following two lines of code from “do_one_cycle.m": [x_e,vx_e,vy_e,vz_e,WEIGHT_e]=merging(xN_e,xG_e,x_e,vx_e,vy_e,vz_e,WEIGHT_e,N_G,constraint_num,size_min) & [x_i,vx_i,vy_i,vz_i,WEIGHT_i]=merging(xN_i,xG_i,x_i,vx_i,vy_i,vz_i,WEIGHT_i,N_G,constraint_num,size_min)
RUn PIC_1D.
Ref:
https://github.com/joshlk/k-means-constrained

Bradley, P. S., K. P. Bennett, and Ayhan Demiriz. "Constrained k-means clustering." Microsoft Research, Redmond (2000): 1-8.

Google's SimpleMinCostFlow C++ implementation

A. Gonoskov, Agnostic conservative down-sampling for optimizing statistical representations and PIC simulations, Computer Physics Communications 271 (2022) 108200. https://doi.org/10.1016/j.cpc.2021.108200.

**Citations**

