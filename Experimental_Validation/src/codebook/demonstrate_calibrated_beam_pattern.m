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
% This script numerically demonstrates the expected radiation beam pattern 
% obtained by different calibration methods.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization
clc; close all

%% ULA Parameters
M = 12;
lambda = 3*10^8/(60.48*10^9);
d = 3.055*10^(-3); %It is measured from the board that the spacing is around 3mm

%% '_d' indicates value in degree and '_r' indicates value in radian
Theta_disired_d = [0]; % Desired Directional Beam
load('Element_Gain.mat'); % load element gain calibration
phase_array_id = 3;
if phase_array_id == 1 || phase_array_id == 3
    Theta_disired_d = -Theta_disired_d;
end
switch  phase_array_id 
    case 1 %% ID 1419-01 Tx
        %calibrate = [1 2 4 3 1 2 3 2 2 3 2 4] - 1;    
    case 2 %% ID 1419-01 Rx
        load('Ave_Error_Phase_Rx.mat');
        load('coarse_calibrate_Rx.mat');
        load('Fine_calibrate_Rx.mat');    
        coarse_calibrate = (coarse_calibrate - 1)';
        phase_error = Ave_Error_Phase;
        element_gain = Element_Gain(5,1,:,:);
        element_gain = squeeze(element_gain);
    case 3 %% ID 0900-02 Tx
        load('Ave_Error_Phase_Tx.mat');
        load('coarse_calibrate_Tx.mat');
        load('Fine_calibrate_Tx.mat');  
        coarse_calibrate = (coarse_calibrate - 1)';
        phase_error = Ave_Error_Phase;
        element_gain = Element_Gain(:,:,5,1);
        element_gain = squeeze(element_gain);
    case 4 %% ID 0900-02 Rx
        %calibrate = [1 2 4 3 1 2 4 2 2 4 2 1] - 1;
end
Theta_r = deg2rad(Theta_disired_d);

%% Storage Vectors
Num_codebook_type = 4;
Codebook = zeros(length(Theta_disired_d),M,Num_codebook_type);
gain_perfect = zeros(length(Theta_disired_d));
gain_quantized = zeros(length(Theta_disired_d));
Theta_r_Channel = -pi/2:0.001:pi/2; 
BeamformerQuant = zeros(M,length(Theta_disired_d),Num_codebook_type);

%% Plot the perfect beam pattern (red) and the quantized beam pattern (blue)
NBeams = length(Theta_r);
r = ceil(sqrt(NBeams));
c = ceil(sqrt(NBeams));
Fig_ID1 = randi([1 10000],1,1);
Fig_ID2 = randi([1 10000],1,1);
Fig_ID3 = randi([1 10000],1,1);
Fig_ID4 = randi([1 10000],1,1);
figure(Fig_ID1);
figure(Fig_ID2);
for nthBeam = 1:NBeams
    Theta_r_Beamformer = Theta_r(nthBeam); %radian of the nth beam pattern
    Beamformer = 1/sqrt(M)*transpose(exp(-1i*2*pi*d/lambda*sin(Theta_r_Beamformer).*(0:1:M-1))); %Perfect beamformer
    Beamformer = 1/sqrt(M)*Beamformer/(Beamformer(1));
    Codebook_Pattern = Quantize(Beamformer, coarse_calibrate); %Quantized beamformer, i.e., in format 012301230123
    [Codebook_Pattern1, Codebook_Pattern2, Codebook_Pattern3, Codebook_Pattern4] = Quantize2(Beamformer, coarse_calibrate, phase_error, Fine_calibrate);
    Codebook(nthBeam,:,1) = Codebook_Pattern1; %Un_calibrated
    Codebook(nthBeam,:,2) = Codebook_Pattern2; %Coarsed_calibrated
    Codebook(nthBeam,:,3) = Codebook_Pattern3; %Fine_calibrated_m1
    Codebook(nthBeam,:,4) = Codebook_Pattern4; %Fine_calibrated_m2
    
    %% Go through all channel to plot the beam patterns
    element_gain1 = zeros(1,12);
    element_gain2 = element_gain1;
    element_gain3 = element_gain1;
    element_gain4 = element_gain1;
    for n = 1:M
        element_gain1(n) = element_gain(n,Codebook_Pattern1(n)+1);
        element_gain2(n) = element_gain(n,Codebook_Pattern2(n)+1);
        element_gain3(n) = element_gain(n,Codebook_Pattern3(n)+1);
        element_gain4(n) = element_gain(n,Codebook_Pattern4(n)+1);
    end
    element_gain1 = element_gain1/norm(element_gain1)*sqrt(M);
    element_gain2 = element_gain2/norm(element_gain2)*sqrt(M);
    element_gain3 = element_gain3/norm(element_gain3)*sqrt(M);
    element_gain4 = element_gain4/norm(element_gain4)*sqrt(M);
    BeamformerQuant(:,nthBeam,1) = 1/sqrt(M)*exp(1j*Codebook_Pattern1*pi/2).*element_gain1'; %Convert quantized format to steering vector
    BeamformerQuant(:,nthBeam,2) = 1/sqrt(M)*exp(1j*Codebook_Pattern2*pi/2)./exp(1j*coarse_calibrate'*pi/2).*element_gain2'; 
    BeamformerQuant(:,nthBeam,3) = 1/sqrt(M)*exp(1j*Codebook_Pattern3*pi/2).*exp(1j*phase_error').*element_gain3';
    phase_error_fine = phase_error;
    for n=1:12
        phase_error_fine(n) = Fine_calibrate(Codebook_Pattern4(n)+1,n);
    end
    BeamformerQuant(:,nthBeam,4) = 1/sqrt(M)*exp(1j*phase_error_fine').*element_gain4';
    
    Gain_Perfect = zeros(1,length(Theta_r_Channel));
    Gain_Quan = zeros(length(Theta_r_Channel),Num_codebook_type);
    for ind_Channel = 1:length(Theta_r_Channel)
        Theta_r_Channel_Realization = Theta_r_Channel(ind_Channel);
        Channel = transpose(exp(-1i*2*pi*d/lambda*sin(Theta_r_Channel_Realization).*(0:1:M-1)));
        Gain_Perfect(ind_Channel) = abs(Beamformer'*Channel);
        Gain_Quan(ind_Channel,1) = abs(BeamformerQuant(:,nthBeam,1)'*Channel);
        Gain_Quan(ind_Channel,2) = abs(BeamformerQuant(:,nthBeam,2)'*Channel);
        Gain_Quan(ind_Channel,3) = abs(BeamformerQuant(:,nthBeam,3)'*Channel);
        Gain_Quan(ind_Channel,4) = abs(BeamformerQuant(:,nthBeam,4)'*Channel);
    end
    Gain_Perfect_dB = 10*log10(Gain_Perfect);
    Gain_Quan_dB = 10*log10(Gain_Quan);
    figure(Fig_ID1);
    subplot(r,c,nthBeam);
    plot(rad2deg(Theta_r_Channel),Gain_Perfect_dB,'r-'); hold on;
    plot(rad2deg(Theta_r_Channel),Gain_Quan_dB(:,1),'b--'); hold on; grid on;
    plot(rad2deg(Theta_r_Channel),Gain_Quan_dB(:,2),'c--'); hold on; grid on;
    plot(rad2deg(Theta_r_Channel),Gain_Quan_dB(:,3),'k--'); hold on; grid on;
    plot(rad2deg(Theta_r_Channel),Gain_Quan_dB(:,4),'m--'); hold on; grid on;
    title(num2str(Theta_disired_d(nthBeam)))
    figure(Fig_ID2);
    plot(rad2deg(Theta_r_Channel),Gain_Perfect_dB,'r-'); hold on;
    plot(rad2deg(Theta_r_Channel),Gain_Quan_dB(:,1),'b--'); hold on; grid on;
    plot(rad2deg(Theta_r_Channel),Gain_Quan_dB(:,2),'c--'); hold on; grid on;
    plot(rad2deg(Theta_r_Channel),Gain_Quan_dB(:,3),'k--'); hold on; grid on;
    plot(rad2deg(Theta_r_Channel),Gain_Quan_dB(:,4),'m--'); hold on; grid on;

end

%% Codebook to Phvec command
%
command=[];
for i=1:1:length(Theta_disired_d)
    command = [command;phvec(Codebook(i,:,1))];
end
Phvec_command_no_calibration = command

%
command=[];
for i=1:1:length(Theta_disired_d)
    command = [command;phvec(Codebook(i,:,2))];
end
Phvec_command_coarse = command
%
command=[];
for i=1:1:length(Theta_disired_d)
    command = [command;phvec(Codebook(i,:,3))];
end
Phvec_command_fine_m1 = command
%
command=[];
for i=1:1:length(Theta_disired_d)
    command = [command;phvec(Codebook(i,:,4))];
end
Phvec_command_fine_m2 = command


%% Save Beam Patterns
switch phase_array_id
    case 3
        save('BeamformerQuant_F.mat','BeamformerQuant')
    case 2
        save('BeamformerQuant_W.mat','BeamformerQuant')
end

%%
Test_Channel_Estimation = 0;
if Test_Channel_Estimation
    load('Ave_Error_Phase_Tx.mat');
    phase_error = Ave_Error_Phase;
    load('coarse_calibrate_Tx.mat');
    coarse_calibrate = (coarse_calibrate - 1)';
    load('Fine_calibrate_Tx.mat');
    load('f_Est_Infinite_Bit.mat')
    f_Est_Infinite_Bit = f_Est_Infinite_Bit/(f_Est_Infinite_Bit(1));
    [CodeBook_F1, CodeBook_F2, CodeBook_F3, CodeBook_F4] = Quantize2(f_Est_Infinite_Bit, coarse_calibrate, phase_error, Fine_calibrate);
    Command_F1 = phvec(CodeBook_F1')
    Command_F2 = phvec(CodeBook_F2')
    Command_F3 = phvec(CodeBook_F3')
    Command_F4 = phvec(CodeBook_F4')
    
    load('Ave_Error_Phase_Rx.mat');
    phase_error = Ave_Error_Phase;
    load('coarse_calibrate_Rx.mat');
    coarse_calibrate = (coarse_calibrate - 1)';
    load('Fine_calibrate_Rx.mat');  
    load('w_Est_Infinite_Bit.mat')
    w_Est_Infinite_Bit = w_Est_Infinite_Bit/(w_Est_Infinite_Bit(1));
    [CodeBook_W1, CodeBook_W2, CodeBook_W3, CodeBook_W4] = Quantize2(w_Est_Infinite_Bit, coarse_calibrate, phase_error, Fine_calibrate);
    Command_W1 = phvec(CodeBook_W1')
    Command_W2 = phvec(CodeBook_W2')
    Command_W3 = phvec(CodeBook_W3')
    Command_W4 = phvec(CodeBook_W4')
end
