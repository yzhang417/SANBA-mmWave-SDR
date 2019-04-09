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
% This script loads the collected data and the estimated results to provide 
% a clear presentation of the performance of the proposed algorithm. It 
% shows the amplitude calibration, which corresponds to Figure 4 in the 
% paper. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization
clc; clear all; close all;


%% Add path
folder_name = ...
[
"/Experimental_Validation";
"/Numerical_Simulation";
];
for i=1:length(folder_name)
    Lib_path = char(strcat(pwd,folder_name(i)));
    addpath(genpath(Lib_path));
end


%% Loading data
load('Element_Gain.mat')
load('Noise_Gain.mat')
Fixed_Side_Antenna_ID = [5];
Fixed_Side_Antenna_Phase = [1];
%RSSI
RSSI = Element_Gain.^2+Noise_Gain.^2;
RSSI_TX = RSSI(:,:,Fixed_Side_Antenna_ID,Fixed_Side_Antenna_Phase);
RSSI_TX_Normalized = RSSI_TX/(max(max(RSSI_TX)));
RSSI_RX = squeeze(RSSI(Fixed_Side_Antenna_ID,Fixed_Side_Antenna_Phase,:,:));
RSSI_RX_Normalized = RSSI_RX/(max(max(RSSI_RX)));
%SNR
SNR_dB= 10*log10((Element_Gain.^2)./Noise_Gain.^2);
SNR_dB_Tx = SNR_dB(:,:,Fixed_Side_Antenna_ID,Fixed_Side_Antenna_Phase);
SNR_dB_Rx = squeeze(SNR_dB(Fixed_Side_Antenna_ID,Fixed_Side_Antenna_Phase,:,:));


%% Plot SNR in dB      
figure
subplot(2,1,1)
bar(SNR_dB_Tx)
xlabel('Index of Tx Antenna Element')
ylabel('SNR (dB)')
title('Receiver SNR with Individual Tx Antenna Element Activated'); grid on;
L1 = legend('Phase State 1', 'Phase State 2', 'Phase State 3', 'Phase State 4');
L1.Orientation = 'vertical';
axis([-inf inf 25 32])

subplot(2,1,2)
bar(SNR_dB_Rx);
xlabel('Index of Rx Antenna Element')
ylabel('SNR (dB)')
title('Receiver SNR with Individual Rx Antenna Element Activated'); grid on;
L2 = legend('Phase State 1', 'Phase State 2', 'Phase State 3', 'Phase State 4');
L2.Orientation = 'vertical';
axis([-inf inf 25 32])


%% Plot RSSI
figure
subplot(2,1,1)
bar(RSSI_TX_Normalized)
xlabel('Index of Tx Antenna Element')
ylabel('RSSI')
title('Normalized RSSI of Each Tx Antenna Element'); grid on;
L1 = legend('Phase State 1', 'Phase State 2', 'Phase State 3', 'Phase State 4');
L1.Orientation = 'vertical';

subplot(2,1,2)
bar(RSSI_RX_Normalized);
xlabel('Index of Rx Antenna Element')
ylabel('RSSI')
title('Normalized RSSI of Each Rx Antenna Element'); grid on;
L2 = legend('Phase State 1', 'Phase State 2', 'Phase State 3', 'Phase State 4');
L2.Orientation = 'vertical';


