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
% This script numerically demonstrates the spectrum efficiency (SE) and 
% mean angle estimation error (MAEE) performance of the proposed 
% non-coherent beam alignment with different number of measurements, in 
% comparison with several benchmarking algorithms. Figure 8 (a) and (b) in 
% the paper are produced by this script.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization 
clc; clear all; close all
format long;
profile on;

%% Add paths
% In this section, three 3rd party software components are included
% 1. SparsePR: Matlab Software for Sparse Phase Retrieval
% 2. Generalized Approximate Message Passing (GAMP) in MRI MATLAB package
% 3. Orthogonal Matching Pursuit (OMP) and Compressive Sampling Matched
% Pursuit (CoSaMP) should be installed
folder_name = ...
[
"Numerical_Simulation/3rd_software_component";
"Numerical_Simulation/src";
"/3rd_software_component";
"/src";
];
for i=1:length(folder_name)
    Lib_path = char(strcat(pwd,folder_name(i)));
    addpath(genpath(Lib_path));
end

%% Uniform linear array (ULA) parameter
ULA_main.lambda = 3*10^8/(60.48*10^9);  % Wavelength of 60.48 GHz signal
ULA_main.d = 3.055*10^(-3);             % Spacing distance between two neighboring antenna elements
ULA_main.Phase_Bit = 2;                 % Number of bits of phase shifters
ULA_main.Nt = 12;                       % Number of transmitter antenna
ULA_main.Nr = 12;                       % Number of receiver antenna
ULA_main.NQt = 4*ULA_main.Nt;           % Initial value of number of AoD quantization
ULA_main.NQr = 4*ULA_main.Nr;           % Initial value of number of AoA quantization
Mtr_main.Mt = 0;                        % Initial value of number of measurements at Tx side
Mtr_main.Mr = 0;                        % Initial value of number of measurements at Rx side

%% Number of measurement range 
M = [ 5  6  7  8  9];               % Number of measurements at Tx and Rx side
G = [25 40 55 60 70];               % Number of AoD and AoA quantization

%% The SNR before beam training in dB
SNR = 0;

%% Setting of searching range in degree
Searching_Area = 40; 

%% Other parameters
Beampattern_Mode = 'Directional_Beam_Angular'; % Available options shown in Generate_Sensing_Matrix.m
Add_Noise_Flag = 1;         % Whether the gaussian noise is introduced in the Rx side
on_grid = 0;                % Whether the generated AoD and AoA is on the quantized grid
L = 1;                      % Number of the dominant path
Rician_K = 5;               % Number of the other paths
Fix_angle_flag = 0;         % Whether fix the AoD and AoA for algorithm debug and test 
Show_Beam_Pattern_Flag = 0; % Whether show the exploited beam patterns in algorithm
Show_leakeage_Flag = 0;     % Whether show the quantization of the generated channel
plot_flag = 0;              % Whether show the recovery performance of the sparse channel vector

%% The different methods (algorithms) to test
Method.PhaseLift = 0;       % High computation and memory requirement
Method.CPRL = 0;            % High computation and memory requirement
Method.PRGAMP = 0;          % Bad performance by directly using PRGAMP algorithm
Method.MySparsePL = 0;      % Bad performance by directly using alternative minimization
Method.PLOMP = 0;           % Proposed algorithm by using OMP
Method.PLGAMP = 1;          % Proposed algorithm by using GAMP
Method.PerfectPhaseCS = 1;  % Perfect phase benchmark
Method.NoisyPhaseCS = 1;    % Noisy phase benchmark
Method.Number = Method.PhaseLift + ...
                Method.CPRL + ...
                Method.PRGAMP + ...
                Method.MySparsePL + ...
                Method.PLOMP + ...
                Method.PLGAMP + ...
                Method.PerfectPhaseCS + ...
                Method.NoisyPhaseCS;
Method.State = [Method.PhaseLift ...
                Method.CPRL ...
                Method.PRGAMP ...   
                Method.MySparsePL ...
                Method.PLOMP ...
                Method.PLGAMP ...
                Method.PerfectPhaseCS ... 
                Method.NoisyPhaseCS];
Method.Name_of_method = ["PhaseLift","CPRL","PRGAMP","MySparsePL","PLOMP",...
                         "PLGAMP","PerfectPhaseCS","NoisyPhaseCS"]; 
Method.Color_Line = ["gd-", "g+-", "yx-", "yo-", "r^-",...
                     "r*-","bo-","mx-"];

%% Number of simulation 
LOOP = 10;
problem_dimension = zeros(LOOP,length(M));
Evaluation = zeros(LOOP,length(M),12,length(Method.Name_of_method));
AoD = zeros(LOOP,1);
AoA = zeros(LOOP,1);

%% Loop of channel generation
tic
parfor loop = 1:LOOP
    loop
    Evaluation_parfor = zeros(length(M),12,length(Method.Name_of_method));
    %% Loop for number of measurements
    for i = 1:length(M)
        
        % Number of measurement    
        Mt = M(i);
        Mr = Mt;
        N_Measurement = Mt*Mr;
        Mtr = Mtr_main;
        Mtr.Mt = Mt;
        Mtr.Mr = Mr;
        ULA = ULA_main;
        ULA.NQt = G(i);
        ULA.NQr = G(i);
        
        % Generate channel with number of path and AoA/AoD setting 
        H = Generate_Channel(ULA, L, Searching_Area, Rician_K, on_grid, Fix_angle_flag);

        % Generate the sparse formulation of channel
        Sparse_Channel_Representation = Sparse_Channel_Formulation(ULA,L,H,Show_leakeage_Flag);
        z = Sparse_Channel_Representation.z;

        % Number of measurement required
        n = length(z);
        s = L;
        
        % Generate Sensing Matrix 
        Sensing_Matrix =...
            Generate_Sensing_Matrix(Beampattern_Mode,Mt,Mr,ULA,...
            Sparse_Channel_Representation,Show_Beam_Pattern_Flag,L,H);
        measurementMat = Sensing_Matrix.measurementMat;

        % Generate Measurement
        [measurements, noise_power] = Generate_Measurement(Sensing_Matrix, SNR, H, Add_Noise_Flag);
        
        % Recovery
        recoveredSig_Set = MyCPR(measurements, measurementMat, z, s, plot_flag, noise_power, Method);
        
        % Evaluation
        ind = 0;
        for name = Method.Name_of_method
            ind = ind + 1;
            if Method.(name)
                Evaluation_parfor(i,:,ind) = ...
                    Evaluation_Recovery(recoveredSig_Set.(name),L,H,...
                    Sparse_Channel_Representation,ULA,SNR,Mtr,noise_power);
            end
        end                       
    end
    Evaluation(loop,:,:,:) = Evaluation_parfor;
end
Running_Time = toc;
profile viewer

%% Save data
Simulation_result_Vs_M.Range = M.^2;
Simulation_result_Vs_M.G = G;
Simulation_result_Vs_M.ULA = ULA_main;
Simulation_result_Vs_M.SNR = SNR;
Simulation_result_Vs_M.L = L;
Simulation_result_Vs_M.Beampattern_Mode = Beampattern_Mode;
Simulation_result_Vs_M.Searching_Area = Searching_Area;
Simulation_result_Vs_M.problem_dimension = problem_dimension;
Simulation_result_Vs_M.Running_Time = Running_Time;
Simulation_result_Vs_M.Mean_Evaluation = mean(Evaluation);
Simulation_result_Vs_M.Method = Method;
Simulation_result_Vs_M.Num_Quantization_Error = Evaluation;
save(['Simulation_result_Vs_M.mat'],'Simulation_result_Vs_M');

%% Plot result
load(['Simulation_result_Vs_M.mat'],'Simulation_result_Vs_M')
Label_Name = 'Number of Measurements $M^2$';
Plot_result(Label_Name, Simulation_result_Vs_M);


