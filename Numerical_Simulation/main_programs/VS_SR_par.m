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
% This script numerically demonstrates how side information helps in
% reducing the beam training overhead. Figure 9 in the paper is produced 
% by this script.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization 
clc; clear all; close all
format long;
LOOP = 20;

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

%% Simulation settings
SNR = 0; % The SNR before beam training in dBB
Searching_Area_Set = 20:10:80;
for i = 1:length(Searching_Area_Set)
    Searching_Area = Searching_Area_Set(i) % Setting of searching range in degree
    switch Searching_Area
        case 20
            M = [ 2  3  4  5];
            G = [25 35 45 55];
        case 30
            M = [ 4  5  6  7];
            G = [25 40 55 60];
        case 40
            M = [ 5  6  7  8  9];
            G = [25 40 55 60 70];
        case 50
            M = [ 6  7  8  9 10 11];
            G = [25 40 45 55 65 70];
        case 60
            M = [ 8  9 10 11 12];
            G = [40 50 55 60 70];
        case 70
            M = [ 9 10 11 12 13];
            G = [40 55 60 70 75];
        case 80
            M = [10 11 12 13 14];
            G = [45 55 60 70 75];
    end
    plot_result = 0;
    sub_VS_SR_par; % Call script sub_VS_SR_par.m
end

%% Load all simulated data and plot the result
MAEE_Set = [0.6, 0.8, 1.0];
M_SR = zeros(length(Searching_Area_Set),length(MAEE_Set),2);
for j = 1:1:length(Searching_Area_Set)
    Searching_Area = Searching_Area_Set(j);
    load(['Simulation_result_Vs_SA_' num2str(Searching_Area) '.mat'],'Simulation_result_Vs_SA');
    %Simulation_result_Vs_SA.G;
    PLGAMP = squeeze(Simulation_result_Vs_SA.Mean_Evaluation(1,:,2,6));
    PerfectPhase = squeeze(Simulation_result_Vs_SA.Mean_Evaluation(1,:,2,7));    
    for i = 1:1:length(MAEE_Set)
        MAEE = MAEE_Set(i);
        [v p] = min(abs(PLGAMP-MAEE));
        M_CPR = Simulation_result_Vs_SA.M(p);
        [v p] = min(abs(PerfectPhase-MAEE));
        M_Perfect = Simulation_result_Vs_SA.M(p);
        M_SR(j,i,1) = M_CPR^2;
        M_SR(j,i,2) = M_Perfect^2;
    end
end

%% Plot result
set(0,'DefaultAxesFontSize',18);
set(0,'DefaultLineMarkerSize',13);
set(0,'DefaultLineLineWidth', 2);

figure
subplot(2,1,2)
plot(Searching_Area_Set,M_SR(:,1,1),'b*-'); hold on;
plot(Searching_Area_Set,M_SR(:,2,1),'ms-'); hold on;
plot(Searching_Area_Set,M_SR(:,3,1),'ko-'); grid on;
set(gca, 'XDir','reverse')
legend_name = ["MAEE $\delta \approx 0.6^{\circ}$","MAEE $\delta \approx 0.8^{\circ}$","MAEE $\delta \approx 1.0^{\circ}$"];
legend(legend_name','Interpreter','latex')
xlabel('Searching Range $\Delta\theta$ in Degree','Interpreter','latex')
ylabel('Number of Measurements $M^2$','Interpreter','latex')

subplot(2,1,1)
for j = 1:1:length(Searching_Area_Set)
    Searching_Area = Searching_Area_Set(j);
    load(['Simulation_result_Vs_SA_' num2str(Searching_Area) '.mat'],'Simulation_result_Vs_SA');
    PLGAMP = squeeze(Simulation_result_Vs_SA.Mean_Evaluation(1,:,2,6));
    PerfectPhase = squeeze(Simulation_result_Vs_SA.Mean_Evaluation(1,:,2,7));
    plot(Simulation_result_Vs_SA.Range,PLGAMP,'ro-'); hold on;
    %plot(Simulation_result_Vs_SA.Range,PerfectPhase,'bo--'); hold on;
end
xrange = [2:1:14].^2;
h(1) = plot(xrange,MAEE_Set(1)*ones(1,length(xrange)),'b:'); hold on;
h(2) = plot(xrange,MAEE_Set(2)*ones(1,length(xrange)),'m:'); hold on;
h(3) = plot(xrange,MAEE_Set(3)*ones(1,length(xrange)),'k:'); grid on;
legend_name = ["MAEE $\delta = 0.6^{\circ}$","MAEE $\delta = 0.8^{\circ}$","MAEE $\delta = 1.0^{\circ}$"];
legend(h,legend_name','Interpreter','latex');
arrow = annotation('arrow');  
arrow.Parent = gca;           
arrow.X = [3,144]; 
arrow.Y = [0.75, 1.75];
arrow.LineWidth  = 2;      
arrow.HeadWidth  = 7;
arrow.HeadLength = 7;
arrow.Color = 'black';
str = ['$\Delta\theta=20^\circ,30^\circ,...,80^\circ$'];
prop = text(arrow.X(2)+3,arrow.Y(2),str,'Interpreter','latex');
prop.FontSize = 16;
xticks(xrange);
xlabel('Number of Measurements $M^2$','Interpreter','latex');
ylabel('MAEE in Degree','Interpreter','latex'); 



