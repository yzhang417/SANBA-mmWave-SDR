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
% shows the PMBC, which corresponds to Figure 5 in the paper. 
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


%% Parameters
angles = [-45:5:45];
ntest = length(angles);
cal_data_folder = char(strcat(pwd,'/Experimental_Validation/data/val_cal_data'));
rssi = zeros(length(angles),4,2);
snr = rssi;
load(fullfile(cal_data_folder,'var_noise_rx_single_on.mat'));
load(fullfile(cal_data_folder,'var_noise_rx_all_on.mat'));


%% read rssi of the calibration result
for side_id = [1 2]
    if side_id == 1
        AoA = 0*ones(1,length(angles));
        AoD = angles;
    else
        AoD = 0*ones(1,length(angles));
        AoA = angles;
    end
    for i=1:1:length(angles)
        mat_name = ['/val_cal_experiment_result_' num2str(AoD(i)) '_' num2str(AoA(i)) '.mat'];
        load(fullfile(cal_data_folder,mat_name));
        rssi(i,1,side_id) = mean(val_cal_result.RSSI(1:3));
        rssi(i,2,side_id) = mean(val_cal_result.RSSI(4:6));
        rssi(i,3,side_id) = mean(val_cal_result.RSSI(7:9));
        rssi(i,4,side_id) = mean(val_cal_result.RSSI(10:12));
    end
end


%% Plot Original
plot_init;
xticklabels_var = {};
for i=1:1:length(angles)
    xticklabels_var{i} = strcat(num2str(angles(i)),'\circ');
end
after_cal_color = 'r*-';
before_cal_color = 'k^-';
mean_after_cal_color = 'r-.';
mean_before_cal_color = 'k--';


%% SNR
snr(:,:,1) = (rssi(:,:,1).^2-var_noise_rx_single_on)/var_noise_rx_single_on;
snr(:,:,2) = (rssi(:,:,2).^2-var_noise_rx_all_on)/var_noise_rx_all_on;
snr_db = 10*log10(snr);

figure
subplot(2,1,1)
b1 = bar(snr_db(:,[1 4],1)); hold on;
b1(1).FaceColor = 'flat';
b1(1).CData(:,:) = repmat([0 0 0],length(snr_db(:,[1 4],1)),1);
b1(2).FaceColor = 'flat';
b1(2).CData(:,:) = repmat([255 0 0],length(snr_db(:,[1 4],1)),1);
mean_snr_db_after_cal = mean(snr_db(:,4,1));
mean_snr_db_before_cal = mean(snr_db(:,1,1));
mean_snr_db_improvement_tx = mean_snr_db_after_cal - mean_snr_db_before_cal
plot(mean_snr_db_before_cal*ones(1,length(angles)),mean_before_cal_color); hold on;
plot(mean_snr_db_after_cal*ones(1,length(angles)),mean_after_cal_color); hold on; grid on;
title('Receiver SNR Before and After Tx Calibration');
xlabel('AoD and Directional Beam Applied at Tx (degree)');
ylabel('SNR (dB)');
legend('Before Calibration','After Calibration','Before Calibration (Mean)', 'After Calibration (Mean)');
xticks([1:length(angles)])
xticklabels(xticklabels_var);
axis([-inf inf 10 35]);
arrow = annotation('arrow');  
arrow.Parent = gca;           
arrow.X = [2.5,2.5]; 
arrow.Y = [mean_snr_db_before_cal,mean_snr_db_after_cal];
arrow.LineStyle = ':';
arrow.LineWidth  = 2;      
arrow.HeadWidth  = 7;
arrow.HeadLength = 7;
arrow.Color = 'blue';
str = ['Averaged SNR Improvement: ' num2str(mean_snr_db_improvement_tx) ' dB'];
prop = text(arrow.X(2)-1,arrow.Y(2)+2.5,str,'FontSize',16,'Color','blue');
prop.LineStyle = ':';
prop.LineWidth = arrow.LineWidth+1;
prop.EdgeColor = 'blue';

subplot(2,1,2)
b2 = bar(snr_db(:,[1 4],2)); hold on;
b2(1).FaceColor = 'flat';
b2(1).CData(:,:) = repmat([0 0 0],length(snr_db(:,[1 4],2)),1);
b2(2).FaceColor = 'flat';
b2(2).CData(:,:) = repmat([255 0 0],length(snr_db(:,[1 4],2)),1);
mean_snr_db_after_cal = mean(snr_db(:,4,2));
mean_snr_db_before_cal = mean(snr_db(:,1,2));
mean_snr_db_improvement_rx = mean_snr_db_after_cal - mean_snr_db_before_cal
plot(mean_snr_db_before_cal*ones(1,length(angles)),mean_before_cal_color); hold on;
plot(mean_snr_db_after_cal*ones(1,length(angles)),mean_after_cal_color); hold on; grid on;
title('Receiver SNR Before and After Rx Calibration');
xlabel('AoA and Directional Beam Applied at Rx');
ylabel('SNR (dB)');
legend('Before Calibration','After Calibration','Before Calibration (Mean)', 'After Calibration (Mean)');
xticks([1:length(angles)])
xticklabels(xticklabels_var);
axis([-inf inf 10 35]);
arrow = annotation('arrow');  
arrow.Parent = gca;           
arrow.X = [2.5,2.5]; 
arrow.Y = [mean_snr_db_before_cal,mean_snr_db_after_cal];
arrow.LineWidth  = 2;
arrow.LineStyle = ':';
arrow.HeadWidth  = 7;
arrow.HeadLength = 7;
arrow.Color = 'blue';
str = ['Averaged SNR Improvement: ' num2str(mean_snr_db_improvement_rx) ' dB'];
prop = text(arrow.X(2)-1,arrow.Y(2)+2.5,str,'FontSize',16,'Color','blue');
prop.LineStyle = ':';
prop.LineWidth = arrow.LineWidth+1;
prop.EdgeColor = 'blue';

