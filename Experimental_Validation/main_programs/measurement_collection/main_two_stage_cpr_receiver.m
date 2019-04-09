%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2019 Yi Zhang and The University of Texas at Austin 
%  
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including 
% without limitation the rights to use, copy, modify, merge, publish, 
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the 
% following conditions:
% 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you use this code or any (modified) part of it in any publication,
% please cite:
%
% Yi Zhang, Kartik Patel, Sanjay Shakkottai, and Robert W. Heath Jr.. 2019. 
% Side-information-aided Non-coherent Beam Alignment Design for Millimeter 
% Wave Systems. In MobiHoc '19: The Twentieth ACM International Symposium 
% on Mobile Ad Hoc Networking and Computing, July 02-05, 2019, Catania, 
% Italy. ACM, New York, NY, USA, 10 pages.
%
% Author: Yi Zhang
% Contact email: yi.zhang.cn@utexas.edu 
% Last modified: Apr. 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script description:
% This script validates the two-stage non-coherent beam alignment algorithm
% proposed in the above paper. 
% Key variables:
% AoD_True: Ground truth of the AoD.
% AoA_True: Ground truth of the AoA.
% num_codebook_to_train: number of codebook used.
% cpr_result: the variable that is saved in the function cal_val_receiver, 
% which stores the estimated AoD and AoA via the proposed algorithm.
% The following variable is store for future demonstration 
% (plotted by functions in folder [Experimental_Validation\main_programs
% \plot_results].
% BeamformerQuant_W: quantized best combiner codebook
% BeamformerQuant_F: quantized best beamformer codebook
% cpr_result:
%    cpr_result.AoA_Estimated: estimated AoA by non-coherent algorithm.
%    cpr_result.AoD_Estimated: estimated AoD by non-coherent algorithm.
%    cpr_result.AoA_Estimated_Sweep: estimated AoA by beam sweeping.
%    cpr_result.AoD_Estimated_Sweep: estimated AoD by beam sweeping.
%    cpr_result.Training_measurement = RSSI of beam sweeping results.
%    cpr_result.Testing_measurement = RSSI of non-coherent result.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Add path
folder_name = ...
[
"/Numerical_Simulation";
"/Experimental_Validation";
];
for i=1:length(folder_name)
    Lib_path = char(strcat(pwd,folder_name(i)));
    addpath(genpath(Lib_path));
end

%%
clear all
t_start = tic;
compileIt = true;
AoD_True = 20;
AoA_True = 0;
calibration_method_id = 4;
for num_codebook_to_train = [8]
    switch num_codebook_to_train
        case 6
            G = 50;
        case 7
            G = 55;
        case 8
            G = 60;
    end
    cpr_receiver(AoD_True, AoA_True, compileIt, calibration_method_id, num_codebook_to_train, G);
    pause();
end
