## Side-information-Aided Non-coherent Beam Alignment (SANBA) for MmWave Systems. 
## The SINBA algorithm is designed and prototyped with MATLAB, USRP and phased arrays.

This repository contains the source codes for the following paper:
 
Yi Zhang, Kartik Patel, Sanjay Shakkottai, and Robert W. Heath Jr.. 2019. 
Side-information-aided Non-coherent Beam Alignment Design for Millimeter 
Wave Systems. In MobiHoc '19: The Twentieth ACM International Symposium on
Mobile Ad Hoc Networking and Computing, July 02-05, 2019, Catania, Italy.
ACM, New York, NY, USA, 10 pages.

In particular, this repository contains the experimental validation and 
numerical simulation of the contribution of the above paper. Details of the
method, experimental results, and representative numerical results can be 
also found in the above paper.

This software was developed at the [Department of Electrical and Computer 
Engineering][UT_ECE], [The University of Texas at Austin][UT_Austin] 
and is released under the MIT license. 

The software was developed and tested using Matlab R2017b and uses 
[SparsePR][SparsePR], a Matlab software package for the efficient 
construction and solution of convex optimization problems. A copy of the 
sparsepr package is included under the third party software directory 
at [Numerical_Simulation/3rd_software_component/]. A selection of scripts 
also uses the CVX optimization software [(link)][cvx], the Generalized 
Approximate Message Passing (GAMP) in MRI MATLAB package [(link)][GAMP] and
the CoSaMP/OMP toolbox [(link)][OMP].


## Directory Structure and Contents

The software package is organized under the following directory structure:
- Experimental_Validation    
     - data/
            - cal_data/
              This folder contains the collected data which will be used to
              calibrate the antennas, i.e. calculating the radiation gain 
              and phase error of each antenna element of a phased array.
            - cal_result/ 
              This folder contains the calibration result after running the 
              functions in [Experimental_Validation/main_programs/
              fine_phased_array_calibration/].
            - val_cal_data/ 
              This folder contains the collected data which will be used to
              validate the proposed calibration method by testing the
              RSSI of the directional beam patterns.
            - val_cpr_data/ 
              This folder contains the collected data which will be used to
              validate the proposed two-stage non-coherent beam alignment
              algorithm.

     - main_programs/
            - cal_decoder/
              This folder contains the specific decoder to decode the  
              collected calibration data to check more details of the 
              performance of the proposed calibration methods.
            - fine_phased_array_calibration/
              This folder contains the main scripts for calculating the
              radiation gain and the phase error of the phase shifters of a
              phased array.
            - general_receiver_decoder/
              This folder contains the general receiver and decoder to run
              the testbed.
            - general_transmitter/
              This folder contains the general transmitter to run the 
              testbed.
            - measurement_collection/
              This folder contains the specific receiver to collect 
              measurements for validating the proposed calibration methods 
              and the non-coherent beam alignment algorithm.
            - plot_results/
              This folder contains the scripts which output the 
              experimental results that are shown in the paper.

     - result_figures/
       This folder contains the experimental results that are shown in the
       paper

     - src/
       This folder contains auxiliary functions necessary to implement the 
       proposed calibration methods and non-coherent algorithms.

- Numerical_Simulation  
     - 3rd_software_component/
       This software contains third party software used by the proposed 
       non-coherent algorithm. In particular, it contains GAMP, sparsepr 
       and an installation file of CoSaMP_OMP toolbox. 

     - main_programs/
       This folder contains the main scripts to reproduce the numerical 
       results that are presented in this paper. In particular, the 
       scripts with postfix _par executes for-loop iterations in parallel 
       on workers.

     - result_figures/
       This folder contains the numerical results that are shown in the
       paper   

     - src/
       This folder contains auxiliary functions necessary to implement the 
       proposed non-coherent algorithm and related benchmarking algorithms.


- Paper_Materials
     - camera_ready/
     - rebuttal/
     - submitted_draft/
     - supplementary_figures/

## Instructions

1. Extract the contents of the zip file.
2. Add folder Experimental_Validation, Numerical_Simulation and their 
   subfolders into Path.
3. Execute the scripts in Experimental_Validation/main_programs/ or 
   Numerical_Simulation/main_programs.


## Contact

Bug reports, comments and suggestions are welcome 
yi.zhang.cn@utexas.edu.


[UT_Austin]: https://www.utexas.edu/
[UT_ECE]: http://www.ece.utexas.edu/
[cvx]: http://cvxr.com/cvx/
[GAMP]: https://gampmatlab.svn.sourceforge.net/svnroot/gampmatlab
[OMP]: https://www.mathworks.com/matlabcentral/fileexchange/32402-cosamp-and-omp-for-sparse-recovery
[SparsePR]: https://bitbucket.org/charms/sparsepr/
[tfocs]: http://cvxr.com/tfocs/

