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
% This script provides the theoretically estimated AoA and AoD by numerical
% simulation with the experimental setting to confirm the correctness of 
% the algorithm implementation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization 
clc; clear all; close all
format long;
profile on;

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

%% ULA parameter
ULA.lambda = 3*10^8/(60.48*10^9);
ULA.d = 3.055*10^(-3);
ULA.Phase_Bit = 2;      % Number of bits of phase shifters
ULA.Nt = 12;             % Number of transmitter antenna
ULA.Nr = 12;           % Number of receiver antenna
num_training = 8;
Mt = num_training;
Mr = num_training;
N_Measurement = Mt*Mr;
Mtr.Mt = Mt;
Mtr.Mr = Mr;
switch num_training
    case 6
        G = 50;
        Searching_Area = 45; 
    case 7
        G = 55;
        Searching_Area = 47.5; 
    case 8
        G = 60;
        Searching_Area = 50; 
end
ULA.NQt = G;       % Number of Quantization of AoA and AoD
ULA.NQr = G;

%% SNR
SNR = 9;

%% Testing angles
Angle_Set = [-20:5:20];

%% Other parameters
Beampattern_Mode = 'ProtoType_Pattern';
Add_Noise_Flag = 1;
on_grid = 0;
L = 1;
Rician_K =0;
Fix_angle_flag = 1;
Show_Beam_Pattern_Flag = 0;
Show_leakeage_Flag = 0;
plot_flag = 0;
problem_dimension = 1;

%% The methods to test
Method.PhaseLift = 0;
Method.CPRL = 0;
Method.PRGAMP = 0;
Method.MySparsePL = 0;
Method.PLOMP = 0;
Method.PLGAMP = 1;
Method.PerfectPhaseCS = 0;
Method.NoisyPhaseCS = 0;
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
Method.Color_Line = ["gd-", "g+-", "cx-", "co-", "r^-",...
                     "r*-","bo-.","b*-."];

%% Number of Simulation 
LOOP = 10;
Evaluation = zeros(LOOP,length(Angle_Set),12,length(Method.Name_of_method),2);

%% Loop of Channel Generation
tic
parfor loop = 1:LOOP
    for side = 1:2
        Evaluation_parfor = zeros(length(Angle_Set),12,length(Method.Name_of_method));
        for i = 1:length(Angle_Set)
            if side == 1
                AoD = 0;
                AoA = Angle_Set(i);
            else
                AoA = 0;
                AoD = Angle_Set(i);
            end   
            % Generate channel with number of path and AoA/AoD setting 
            Rician_K = 0;
            H = Generate_Channel_fix_angle(ULA, L, Searching_Area, Rician_K, on_grid, Fix_angle_flag, AoD, AoA);

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
        Evaluation(loop,:,:,:,side) = Evaluation_parfor;
    end
end
Running_Time = toc;
profile viewer

%% Save data
Simulation_result_AoD.Range = Angle_Set;
Simulation_result_AoD.ULA = ULA;
Simulation_result_AoD.SNR = SNR;
Simulation_result_AoD.L = L;
Simulation_result_AoD.Beampattern_Mode = Beampattern_Mode;
Simulation_result_AoD.Searching_Area = Searching_Area;
Simulation_result_AoD.problem_dimension = problem_dimension;
Simulation_result_AoD.N_Measurement = N_Measurement;
Simulation_result_AoD.Running_Time = Running_Time;
Simulation_result_AoD.Mean_Evaluation = mean(Evaluation(:,:,:,:,2));
Simulation_result_AoD.Method = Method;
Simulation_result_AoD.Num_Quantization_Error = Evaluation;
Simulation_result_AoA = Simulation_result_AoD;
Simulation_result_AoA.Mean_Evaluation = mean(Evaluation(:,:,:,:,1));

%% Plot the theoretical estimated AoD and AoA
Label_Name = 'Testing angle of AoD in degree';
Plot_result(Label_Name, Simulation_result_AoD);
Label_Name = 'Testing angle of AoA in degree';
Plot_result(Label_Name, Simulation_result_AoA);

%% Save the theoretically estimated AoD and AoA
AoD_estimated_theory = Simulation_result_AoD.Mean_Evaluation(1,:,2,6);
AoA_estimated_theory = Simulation_result_AoA.Mean_Evaluation(1,:,2,6);
cpr_data_folder = char(strcat(pwd,'/Experimental_Validation/data/val_cpr_data'));
save([fullfile(cpr_data_folder,['AoD_estimated_theory_M' num2str(num_training) '.mat'])],'AoD_estimated_theory');
save([fullfile(cpr_data_folder,['AoA_estimated_theory_M' num2str(num_training) '.mat'])],'AoA_estimated_theory');

