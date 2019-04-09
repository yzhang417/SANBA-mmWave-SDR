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
% This script validates the fine phased array calibration proposed in the
% above paper. The related coarse calibration method is also evaluated. We 
% validate all the calibration methods by demonstrating the their RSSI 
% performance of the directional beams.
% Key variables:
% AoD_True: Ground truth of the AoD.
% AoA_True: Ground truth of the AoA.
% cpr_result: a structure array detailed as below:
%    cal_val_result.AoA: true AoA.
%    cal_val_result.AoD: true AoD.
%    cal_val_result.RSSI: RSSI received by different calibration methods,
%    including the major baseline where no calibration is performed.
% Four involved calibration methods are:
% 1. without any calibration (baseline).
% 2. coarse phased array calibration: in this method, the four phase states
%    are considered to be exactly 0, 90, 180, 270. The results are stored
%    in coarse_calibrate_Tx.mat and coarse_calibrate_Rx.mat.
%    For this method, the phase error is stored in coarse_calibrate_Tx.mat 
%    and coarse_calibrate_Rx.mat in folder [Experimental_Validation/data/
%    cal_result].
%    For this method, coarse_calibrate_Tx.mat/coarse_calibrate_Rx.mat only 
%    have four possibilities of values which are four types of error 0, 90, 
%    180 and 270 degrees.
% 3. fine phased array calibration method 1: in this method, the four phase
%    states are considered to be 0+epsilon_n, 90+epsilon_n, 180+epsilon_n,
%    and 270+epsilon_n for the n-th antenna element. Thus, the phase error
%    is antenna dependent but independent of the phase state.
%    For this method, the phase error is stored in Ave_Error_Phase_Tx.mat 
%    and Ave_Error_Phase_Rx.mat in folder [Experimental_Validation/data/
%    cal_result]. For this method, the Ave_Error_Phase_Tx.mat and 
%    Ave_Error_Phase_Rx.mat store the estimated phase shift error epsilon_n 
%    which needs to be added to the nominal phase 0, 90, 180 or 270 to be 
%    the real phase corresponding to the nominal four types 00, 01, 10, 11 
%    phase states.
% 4. fine phased array calibration method 2: in this method, the four phase
%    states are considered to be 0+epsilon_n,1, 90+epsilon_n,2, 
%    180+epsilon_n,3, and 270+epsilon_n,4. This means the phase error is
%    not only antenna dependent but also phase state dependent, which is
%    a more realistic model.
%    For this method, the phase error is stored in Fine_calibrate_Tx.mat 
%    and Fine_calibrate_Rx.mat in folder [Experimental_Validation/data/
%    cal_result]. In particular, the Fine_calibrate_Tx.mat and 
%    Fine_calibrate_Rx.mat store the estimated phase of the nominal four 
%    types 00, 01, 10, 11 phase states.
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


%% Validate calibrated directional beam
clear all
t_start = tic;
compileIt = true;
AoD_True = 5; 
AoA_True = 0; 
cal_val_result = cal_val_receiver(AoD_True, AoA_True, compileIt);
toc(t_start);


%% Noise floor estimation
% clear all
% compileIt = false;
% Range_Set = [-18.7500 -11.2500 -3.7500 3.7500 11.2500 18.7500];
% rssi_nosie = zeros(1,length(Range_Set));
% 
% %
% AoD_True = 0;
% for i = 1:length(Range_Set)
%     AoA_True = Range_Set(i);
%     cal_val_result = cal_val_receiver(AoD_True, AoA_True, compileIt);
%     rssi_nosie(i) = mean(cal_val_result.RSSI);
% end
% var_noise_rx_all_on = mean(rssi_nosie)^2;
% save('var_noise_rx_all_on.mat','var_noise_rx_all_on');
% 
% %
% AoA_True = 0;
% for i = 1:length(Range_Set)
%     AoD_True = Range_Set(i);
%     cal_val_result = cal_val_receiver(AoD_True, AoA_True, compileIt);
%     rssi_nosie(i) = mean(cal_val_result.RSSI);
% end
% var_noise_rx_single_on = mean(rssi_nosie)^2;
% save('var_noise_rx_single_on.mat','var_noise_rx_single_on');

