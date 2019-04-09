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
% This function performs three major tasks illustrated as below:
% Task I: running the Rx USRP to receive the signal over the air and 
% fetching the digital data stream from the Rx USRP. To this end the
% necessary Tx/Rx phased array configurations are performed. 
% Task II: after collecting the measurements, this function runs the 
% proposed two-stage non-coherent beam alignment algorithm to estimate the 
% AoA and AoD. 
% Task III: to further evaluate and justify the correct implementation. 
% This function test the estimated best beamformer and combiner pair
% over the air and compare the RSSI with the beam sweeping strategy.
% The collected data is saved in 
% [Experimental_Validation\data\plot_results].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% MeasurementsPrototype: the non-coherent measurements collected over the
% air.
% num_codebook_to_train: number of beam patterns to be probed at each side.
% prmCode: 
%   prmCode.Shift: the number of measurement to be omitted.
%   prmCode.calibration_method_id: the calibration method.
% Searching_Area: searching range provided by side-information.
% G: quantization level.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output arguments:
% AoA_Estimated: estimated AoA by non-coherent algorithm.
% AoD_Estimated: estimated AoD by non-coherent algorithm.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [AoD_Estimated, AoA_Estimated] = run_cpr(MeasurementsPrototype, num_codebook_to_train, prmCode, Searching_Area, G)
    %% ULA parameter
    ULA.lambda = 3*10^8/(60.48*10^9);
    ULA.d = 3.055*10^(-3);
    ULA.Phase_Bit = 2;      % Number of bits of phase shifters
    ULA.Nt = 12;            % Number of transmitter antenna
    ULA.Nr = 12;            % Number of receiver antenna
    ULA.NQt = G;            % Quantization of AoA and AoD
    ULA.NQr = ULA.NQt;

    %% SNR
    SNR = 25;

    %% Other parameters
    Beampattern_Mode = 'ProtoType_Pattern';
    Add_Noise_Flag = 0;
    on_grid = 0;
    L = 1;
    Rician_K = 0;
    Fix_angle_flag = 0;
    Show_Beam_Pattern_Flag = 0;
    Show_leakeage_Flag = 0;
    plot_flag = 0;

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


    %% Generate Channel with number of path and AoA/AoD setting 
    H = Generate_Channel(ULA, L, Searching_Area, Rician_K, on_grid, Fix_angle_flag);

    %% Generate the sparse formulation of channel
    Sparse_Channel_Representation = ...
        Sparse_Channel_Formulation(ULA,L,H,Show_leakeage_Flag);
    z = Sparse_Channel_Representation.z_leakage_reduced;

    %% Number of measurement required
    n = length(z);
    s = L;
    Mt = num_codebook_to_train;
    Mr = num_codebook_to_train;
    Mtr.Mt = Mt;
    Mtr.Mr = Mr;

    %% Generate Sensing Matrix 
    Sensing_Matrix =...
        Generate_Sensing_Matrix(Beampattern_Mode,Mt,Mr,ULA,...
        Sparse_Channel_Representation,Show_Beam_Pattern_Flag,L,H);
    measurementMat = Sensing_Matrix.measurementMat;

    %% Generate Measurement
    [measurements, noise_power] = ...
        Generate_Measurement(Sensing_Matrix, SNR, H, Add_Noise_Flag);
    if strcmp(Beampattern_Mode,"ProtoType_Pattern")
        %bar(measurements.measurements_norm_square); hold on;
        measurements_prototype = extract_measurement(prmCode,MeasurementsPrototype,Mt,Mr);
        %measurements_prototype = measurements_prototype - min(measurements_prototype);
        measurements.measurements_norm_square = measurements_prototype/norm(measurements_prototype)*norm(measurements.measurements_norm_square);
        measurements.measurements_with_noisy_phase = sqrt(measurements.measurements_norm_square);
        %figure
        %stem(measurements.measurements_norm_square,'r*');
        %legend('Measured RSSI')
        %xlabel('Measurement Index')
        %ylabel('RSSI')
    end

    %% Recovery
    recoveredSig_Set = MyCPR(measurements, measurementMat, z, s, plot_flag, noise_power, Method);

    %% Metric Evaluation
    Evaluation = zeros(12,length(Method.Name_of_method));
    ind = 0;
    for name = Method.Name_of_method
        ind = ind + 1;
        if Method.(name)
            Evaluation(:,ind) = ...
                Evaluation_Recovery(recoveredSig_Set.(name),L,H,Sparse_Channel_Representation,ULA,SNR,Mtr,noise_power);
        end
    end      

    %% display the estiamtion result if required
    load('AoD_Estimated.mat')
    AoD_Estimated = AoD_Estimated;
    load('AoA_Estimated.mat')
    AoA_Estimated = AoA_Estimated;
    load('AoDA_Err.mat')
    AoDA_Err = AoDA_Err;
end
