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
% a clear presentation of the proposed two-stage non-coherent beam 
% alignment algorithm, which corresponds to Figure 7 in the paper. There 
% are two options for variable M: 6 and 8.
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
calibration_method_id = 4;
M = 6;
testing_angle_set = -20:5:20;
AoD_Fix = 0;
AoA_Fix = 0;
ntest = length(testing_angle_set);
cpr_data_folder = char(strcat(pwd,'/Experimental_Validation/data/val_cpr_data'));
MAEE = zeros(length(testing_angle_set),2,2);
SNR = zeros(length(testing_angle_set),2,2);


%% Noise floor
load(fullfile(cpr_data_folder,'var_noise_rx_all_on.mat'));
noise_floor = var_noise_rx_all_on;


%% Calculate the MAEE and SNR improvement
for i=1:1:ntest
    for j= 1:2
        if j==2
            mat_name = ['/CPR_experiment_result_' num2str(AoD_Fix) '_' num2str(testing_angle_set(i)) '_calM' num2str(calibration_method_id) '_M' num2str(M) '.mat'];
        else
            mat_name = ['/CPR_experiment_result_' num2str(testing_angle_set(i)) '_' num2str(AoA_Fix) '_calM' num2str(calibration_method_id) '_M' num2str(M) '.mat'];
        end
        matfile = fullfile(cpr_data_folder,mat_name);
        load(matfile)
        if j==2
            AME_CPR = (abs(cpr_result.AoA_Estimated-testing_angle_set(i))+abs(cpr_result.AoD_Estimated-AoD_Fix))/2;
            AME_Sweep = (abs(cpr_result.AoA_Estimated_Sweep-testing_angle_set(i))+abs(cpr_result.AoD_Estimated_Sweep-AoD_Fix))/2;
        else
            AME_CPR = (abs(cpr_result.AoA_Estimated-AoA_Fix)+abs(cpr_result.AoD_Estimated-testing_angle_set(i)))/2;
            AME_Sweep = (abs(cpr_result.AoA_Estimated_Sweep-AoA_Fix)+abs(cpr_result.AoD_Estimated_Sweep-testing_angle_set(i)))/2;
        end
        MAEE(i,j,1) = AME_Sweep;
        MAEE(i,j,2) = AME_CPR;
        [v1,p1] = max(Measurement_train_test(1:M,:));
        [v2, p2] = max(v1);
        r = p1(p2);
        c = p2;
        SNR(i,j,1) = ((max(max(Measurement_train_test(1:M,:))).^2-noise_floor))/noise_floor;
        [v1,p1] = max(Measurement_train_test(M+1:M+2,:));
        [v2, p2] = max(v1);
        r = p1(p2);
        c = p2;
        SNR(i,j,2) = ((max(max(Measurement_train_test(M+1:M+2,:))).^2-noise_floor))/noise_floor;
    end
end


%% Visualization
plot_init;
xticklabels_var = {};
for i=1:1:length(testing_angle_set)
    xticklabels_var{i} = strcat(num2str(testing_angle_set(i)),'\circ');
end
cpr_color_theory = 'r*--';
sweep_color_theory = 'k^--';
cpr_color = 'ro-';
sweep_color = 'ks-';
mean_cpr_color = 'r-.';
mean_sweep_color = 'k--';

%MAEE
load(['AoD_estimated_theory_M' num2str(M) '.mat']);
load(['AoA_estimated_theory_M' num2str(M) '.mat']);
    
M = num_codebook_to_train;
switch M
    case 6
        Searching_Area = 45;
    case 7
        Searching_Area = 47.5;
    case 8
        Searching_Area = 50;
end
Angle_range = [-Searching_Area/2 Searching_Area/2];
Partition_Angle_range = linspace(Angle_range(1),Angle_range(2),M+1);
Training_Angles = (Partition_Angle_range(1:end-1)+Partition_Angle_range(2:end))/2;
Beam_sweep_estimated_theory = zeros(size(testing_angle_set));
for i=1:length(testing_angle_set)
    [Beam_sweep_estimated_theory(i), p] = min(abs(Training_Angles-testing_angle_set(i)));
    Beam_sweep_estimated_theory(i) = (Beam_sweep_estimated_theory(i) + min(abs(Training_Angles-0)))/2;
end

figure
subplot(1,2,1)
plot(MAEE(:,1,1),sweep_color); hold on;
plot(Beam_sweep_estimated_theory,sweep_color_theory);hold on;
plot(MAEE(:,1,2),cpr_color); hold on;
plot(AoA_estimated_theory,cpr_color_theory);hold on; 
grid on;

title(['MAEE with AoA = $0^\circ$ and $M = ' num2str(M) '$'],'Interpreter','latex')
xlabel('AoD (degree)');
ylabel('MAEE (degree)')
legend('Beam Sweeping (Experiment)',...
       'Beam Sweeping (Simulation)',...
       'Non-Coherent Estimation (Experiment)',...
       'Non-Coherent Estimation (Simulation)');
xticks([1:length(testing_angle_set)])
xticklabels(xticklabels_var);

subplot(1,2,2)
plot(MAEE(:,2,1),sweep_color); hold on;
plot(Beam_sweep_estimated_theory,sweep_color_theory);hold on;
plot(MAEE(:,2,2),cpr_color); hold on;
plot(AoD_estimated_theory,cpr_color_theory);hold on;grid on;
title(['MAEE with AoD = $0^\circ$ and $M = ' num2str(M) '$'],'Interpreter','latex')
xlabel('AoA (degree)');
ylabel('MAEE (degree)')
legend('Beam Sweeping (Experiment)',...
       'Beam Sweeping (Simulation)',...
       'Non-Coherent Estimation (Experiment)',...
       'Non-Coherent Estimation (Simulation)');   
xticks([1:length(testing_angle_set)])
xticklabels(xticklabels_var);
axis([-inf inf 0 5])

%SNR
SNRdB = 10*log10(SNR);
figure
subplot(1,2,1)
b1 = bar(squeeze(SNRdB(:,1,:))); hold on;
b1(1).FaceColor = 'flat';
b1(1).CData(:,:) = repmat([0 0 0],length(SNRdB(:,1,:)),1);
b1(2).FaceColor = 'flat';
b1(2).CData(:,:) = repmat([255 0 0],length(SNRdB(:,1,:)),1);
Mean_Sweep = mean(SNRdB(:,1,1));
plot(Mean_Sweep*ones(1,length(testing_angle_set)),mean_sweep_color);hold on;
Mean_CPR = mean(SNRdB(:,1,2));
plot(Mean_CPR*ones(1,length(testing_angle_set)),mean_cpr_color);hold on; grid on;
title(['Receiver SNR After Beam Training (AoA = $0^\circ$ and $M = ' num2str(M) '$)'],'Interpreter','latex')
xlabel('AoD (degree)');
ylabel('SNR (dB)');
legend('Beam Sweeping',...
       'Non-Coherent Estimation', ...
       'Beam Sweeping (Mean)',...
       'Non-Coherent Estimation (Mean)');
xticks([1:length(testing_angle_set)])
xticklabels(xticklabels_var);
axis([-inf inf -inf inf]);
arrow = annotation('arrow');  
arrow.Parent = gca;           
arrow.X = [3.5,3.5]; 
arrow.Y = [Mean_Sweep,Mean_CPR];
arrow.LineWidth  = 3;      
arrow.HeadWidth  = 15;
arrow.HeadLength = 15;
arrow.Color = 'blue';
arrow.LineStyle = ':';
str = ['Averaged SNR Improvement: ' num2str(Mean_CPR-Mean_Sweep) ' dB'];
prop = text(arrow.X(2)-1,arrow.Y(2)+1,str,'FontSize',16,'Color','blue');
prop.LineStyle = ':';
prop.LineWidth = arrow.LineWidth+1;
prop.EdgeColor = 'blue';
axis([-inf inf 25 32])

subplot(1,2,2)
b2 = bar(squeeze(SNRdB(:,1,:))); hold on;
b2(1).FaceColor = 'flat';
b2(1).CData(:,:) = repmat([0 0 0],length(SNRdB(:,1,:)),1);
b2(2).FaceColor = 'flat';
b2(2).CData(:,:) = repmat([255 0 0],length(SNRdB(:,1,:)),1);
Mean_Sweep = mean(SNRdB(:,2,1));
plot(Mean_Sweep*ones(1,length(testing_angle_set)),mean_sweep_color);hold on;
Mean_CPR = mean(SNRdB(:,2,2));
plot(Mean_CPR*ones(1,length(testing_angle_set)),mean_cpr_color);hold on; grid on;
title(['Receiver SNR After Beam Training (AoD = $0^\circ$ and $M = ' num2str(M) '$)'],'Interpreter','latex')
xlabel('AoA (degree)');
ylabel('SNR (dB)');
legend('Beam Sweeping',...
       'Non-Coherent Estimation', ...
       'Beam Sweeping (Mean)',...
       'Non-Coherent Estimation (Mean)');
xticks([1:length(testing_angle_set)])
xticklabels(xticklabels_var);
xticks([1:length(testing_angle_set)])
xticklabels(xticklabels_var);
axis([-inf inf -inf inf]);
arrow = annotation('arrow');  
arrow.Parent = gca;           
arrow.X = [2.5,2.5]; 
arrow.Y = [Mean_Sweep,Mean_CPR];
arrow.LineWidth  = 3;      
arrow.HeadWidth  = 15;
arrow.HeadLength = 15;
arrow.Color = 'blue';
arrow.LineStyle = ':';
str = ['Averaged SNR Improvement: ' num2str(Mean_CPR-Mean_Sweep) ' dB'];
prop = text(arrow.X(2)-1,arrow.Y(2)+1,str,'FontSize',16,'Color','blue');
prop.LineStyle = ':';
prop.LineWidth = arrow.LineWidth+1;
prop.EdgeColor = 'blue';
axis([-inf inf 25 32])

% %SNR
% SNRdB = 10*log10(SNR);
% figure
% subplot(1,2,1)
% plot(SNRdB(:,1,1),sweep_color); hold on;
% Mean_Sweep = mean(SNRdB(:,1,1));
% plot(Mean_Sweep*ones(1,length(testing_angle_set)),mean_sweep_color);hold on;
% plot(SNRdB(:,1,2),cpr_color); hold on;
% Mean_CPR = mean(SNRdB(:,1,2));
% plot(Mean_CPR*ones(1,length(testing_angle_set)),mean_cpr_color);hold on; grid on;
% title(['Receiver SNR After Beam Training (AoA = $0^\circ$ and $M = ' num2str(M) '$)'],'Interpreter','latex')
% xlabel('AoD (degree)');
% ylabel('SNR (dB)');
% legend('Beam Sweeping',...
%        'Beam Sweeping (Mean)',...
%        'Non-Coherent Estimation', ...
%        'Non-Coherent Estimation (Mean)');
% xticks([1:length(testing_angle_set)])
% xticklabels(xticklabels_var);
% axis([-inf inf -inf inf]);
% arrow = annotation('arrow');  
% arrow.Parent = gca;           
% arrow.X = [3.5,3.5]; 
% arrow.Y = [Mean_Sweep,Mean_CPR];
% arrow.LineWidth  = 3;      
% arrow.HeadWidth  = 15;
% arrow.HeadLength = 15;
% arrow.Color = 'blue';
% arrow.LineStyle = ':';
% str = ['Averaged SNR Improvement: ' num2str(Mean_CPR-Mean_Sweep) ' dB'];
% prop = text(arrow.X(2)-1,arrow.Y(1)-0.5,str,'FontSize',16,'Color','blue');
% prop.LineStyle = ':';
% prop.LineWidth = arrow.LineWidth+1;
% prop.EdgeColor = 'blue';
% axis([-inf inf 26 30])
% 
% subplot(1,2,2)
% plot(SNRdB(:,2,1),sweep_color); hold on;
% Mean_Sweep = mean(SNRdB(:,2,1));
% plot(Mean_Sweep*ones(1,length(testing_angle_set)),mean_sweep_color);hold on;
% plot(SNRdB(:,2,2),cpr_color); hold on;
% Mean_CPR = mean(SNRdB(:,2,2));
% plot(Mean_CPR*ones(1,length(testing_angle_set)),mean_cpr_color);hold on; grid on;
% title(['Receiver SNR After Beam Training (AoD = $0^\circ$ and $M = ' num2str(M) '$)'],'Interpreter','latex')
% xlabel('AoA (degree)');
% ylabel('SNR (dB)');
% legend('Beam Sweeping',...
%        'Beam Sweeping (Mean)',...
%        'Non-Coherent Estimation', ...
%        'Non-Coherent Estimation (Mean)');
% xticks([1:length(testing_angle_set)])
% xticklabels(xticklabels_var);
% xticks([1:length(testing_angle_set)])
% xticklabels(xticklabels_var);
% axis([-inf inf -inf inf]);
% arrow = annotation('arrow');  
% arrow.Parent = gca;           
% arrow.X = [2.5,2.5]; 
% arrow.Y = [Mean_Sweep,Mean_CPR];
% arrow.LineWidth  = 3;      
% arrow.HeadWidth  = 15;
% arrow.HeadLength = 15;
% arrow.Color = 'blue';
% arrow.LineStyle = ':';
% str = ['Averaged SNR Improvement: ' num2str(Mean_CPR-Mean_Sweep) ' dB'];
% prop = text(arrow.X(2)-1,arrow.Y(1)-0.5,str,'FontSize',16,'Color','blue');
% prop.LineStyle = ':';
% prop.LineWidth = arrow.LineWidth+1;
% prop.EdgeColor = 'blue';
% axis([-inf inf 26 30])

% %Spectrum Efficiency
% rate = log2(1+SNR);
% figure
% subplot(2,1,1)
% plot(rate(:,1,1),sweep_color); hold on;
% plot(mean(rate(:,1,1))*ones(1,length(testing_angle_set)),mean_sweep_color);hold on;
% plot(rate(:,1,2),cpr_color); hold on;
% plot(mean(rate(:,1,2))*ones(1,length(testing_angle_set)),mean_cpr_color);hold on; grid on;
% title(['Spectrum Efficiency After Beam Training (AoA = $0^\circ$ and $M = ' num2str(M) '$)'],'Interpreter','latex')
% xlabel('AoD (degree)');
% ylabel('SE (bits/s/Hz)');
% legend('Beam Sweeping',...
%        'Beam Sweeping (Mean)',...
%        'Non-Coherent Estimation', ...
%        'Non-Coherent Estimation (Mean)');
% xticks([1:length(testing_angle_set)])
% xticklabels(xticklabels_var);
% 
% subplot(2,1,2)
% plot(rate(:,2,1),sweep_color); hold on;
% plot(mean(rate(:,2,1))*ones(1,length(testing_angle_set)),mean_sweep_color);hold on;
% plot(rate(:,2,2),cpr_color); hold on;
% plot(mean(rate(:,2,2))*ones(1,length(testing_angle_set)),mean_cpr_color);hold on; grid on;
% title(['Spectrum Efficiency After Beam Training (AoD = $0^\circ$ and $M = ' num2str(M) '$)'],'Interpreter','latex')
% xlabel('AoA (degree)');
% ylabel('SE (bits/s/Hz)');
% legend('Beam Sweeping',...
%        'Beam Sweeping (Mean)',...
%        'Non-Coherent Estimation', ...
%        'Non-Coherent Estimation (Mean)');
% xticks([1:length(testing_angle_set)])
% xticklabels(xticklabels_var);



