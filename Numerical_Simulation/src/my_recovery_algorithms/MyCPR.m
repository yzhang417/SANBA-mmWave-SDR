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
% Function description:
% This function manages all the recovery (compressive phase retrieval) 
% algorithms and the related benchmarking algorithms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% measurements: a structure array that groups the different types of  
% received noised measurements. Please refer to Generate_Measurement.m for 
% the details of its fields.
% measurementMat: sensing matrix, it is stored in the field.
% "measurementMat" of the output of the function Generate_Sensing_Matrix
% z: vectorized and quantized mmWave channel. it is stored in the field "z" 
% of the output of the function Sparse_Channel_Formulation.
% s: level of sparsity, i.e. number of large dominant paths. 
% plot_flag: whether show the recovery performance of the sparse channel 
% vector.
% noise_power: power of noise.
% Method: a structure array that groups related parameters on the recovery
% algorithms that needs to be tested.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output arguments:
% recoveredSig_Set: a structure array that groups all the recovered sparse
% signal (estimated z) by different methods. Each field of recoveredSig_Set
% is a column vector.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function recoveredSig_Set = MyCPR(measurements, measurementMat, z, s, plot_flag, noise_power, Method)
%% Problem parameters
measurements_norm_square = measurements.measurements_norm_square;
measurements_with_perfect_phase = measurements.measurements_with_perfect_phase;
measurements_with_noisy_phase = measurements.measurements_with_noisy_phase;
TrueSig = z;
[m, n] = size(measurementMat);
if plot_flag % Print out the key parameters of the problem to be solved
    fprintf('Problem dimension - %d\n', n);
    fprintf('Sparsity - %d\n', s );
    fprintf('No. of measurements_norm_square - %d\n', m);
end
recoveredSig_Set.PerfectPhaseCS = zeros(n,1);


%% PhaseLift
if Method.PhaseLift
    recoveredSig_Set.PhaseLift = MyPhaseLift(measurements_norm_square, measurementMat);
    if plot_flag
        Plot_Recovery_Performance(TrueSig,recoveredSig_Set.PhaseLift,'PhaseLift');
    end
end


%% CPRL
if Method.CPRL
    recoveredSig_Set.CPRL = MyCPRL(measurements_norm_square, measurementMat);
    if plot_flag
        Plot_Recovery_Performance(TrueSig,recoveredSig_Set.CPRL,'CPRL');
    end
end


%% PR-GAMP
if Method.PRGAMP
    opt_pr.xreal = 0;
    opt_pr.xnonneg = 0;
    opt_pr.tuneOn = 1;
    opt_pr.sparseRat = s/n;
    opt_pr.SNRdBinit = 20;
    opt_pr.SNRdB = 25;
    opt_pr.verbose = 0;
    opt_pr.maxTry = 100;
    opt_pr.nitEM = 50;     
    recoveredSig_Set.PRGAMP = MyPRGAMP(measurements_norm_square, measurementMat, opt_pr);
    if plot_flag
        Plot_Recovery_Performance(TrueSig,recoveredSig_Set.PRGAMP,'PRGAMP');
    end
end


%% SparsePhaseLift
if Method.MySparsePL
    opt_AltMin.MaxIte = 1e4;
    opt_AltMin.s = s;
    opt_AltMin.opt_epslion = 0.001;
    recoveredSig_Set.MySparsePL = ...
        MySparsePL(measurements_norm_square, measurementMat, opt_AltMin);
    if plot_flag
        Plot_Recovery_Performance(TrueSig,recoveredSig_Set.MySparsePL,'MySparsePL');
    end
end


%% Two_Stage_Recovery (Our proposed algorithm)
if Method.PLOMP || Method.PLGAMP
    [recoveredSig_Set.PLOMP, recoveredSig_Set.PLGAMP] = ...
        My_TwoStage_Recovery(measurements_norm_square, measurementMat, s, noise_power, plot_flag, Method);
    if plot_flag && Method.PLOMP
        Plot_Recovery_Performance(TrueSig,recoveredSig_Set.PLOMP,'PLOMP');
    end
    if plot_flag && Method.PLGAMP
        Plot_Recovery_Performance(TrueSig,recoveredSig_Set.PLGAMP,'PLGAMP');
    end
end


%% Conventional CS via measurement with perfect phase information
if Method.PerfectPhaseCS
    recoveredSig_Set.PerfectPhaseCS = My_Conventional_CS(measurements_with_perfect_phase, measurementMat, s, noise_power);
    if plot_flag
        Plot_Recovery_Performance(TrueSig,recoveredSig_Set.PerfectPhaseCS,'PerfectPhaseCS');
    end
end


%% Conventional CS via measurement with noisy phase information
if Method.NoisyPhaseCS
    recoveredSig_Set.NoisyPhaseCS = My_Conventional_CS(measurements_with_noisy_phase, measurementMat, s, noise_power);
    if plot_flag
        Plot_Recovery_Performance(TrueSig,recoveredSig_Set.NoisyPhaseCS,'NoiseyPhaseCS');
    end
end
