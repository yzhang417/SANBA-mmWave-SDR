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
% This script numerically demonstrates the ideal beam patterns.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization
clc
clear all
close all

%% ULA
M = 12;
lambda = 3*10^8/(60.48*10^9);
%lambda = 3*10^8/(62.58*10^9);
d = 3.055*10^(-3); %It is measured from the board that the spacing is around 3mm

%% Calibration Referece
calibrate = [0 0 0 0 0 0 0 0 0 0 0 0];
ind_beamformer = 1;

%% Coverage
Theta_d_Max = 30;

%% AoD, in the following, '_d' indicates value in degree and '_r' indicates value in radian
if  ind_beamformer == 1
    Theta_disired_d = Theta_d_Max:-10:-Theta_d_Max;
else
    Theta_disired_d = -Theta_d_Max:10:Theta_d_Max;
end 
Theta_r = deg2rad(Theta_disired_d); %90 degree corresponds to 000000000000, i.e. 0x000000 phvec command


%% Storage Vectors
Codebook = zeros(length(Theta_disired_d),12);
gain_perfect = zeros(length(Theta_disired_d));
gain_quantized = zeros(length(Theta_disired_d));
Theta_r_Channel = -pi/2:0.001:pi/2; 

%% Plot the perfect beam pattern (red) and the quantized beam pattern (blue)
for nthBeam = 1:length(Theta_r)
    Theta_r_Beamformer = Theta_r(nthBeam); %radian of the nth beam pattern
    Beamformer = 1/sqrt(M)*transpose(exp(-1i*2*pi*d/lambda*sin(Theta_r_Beamformer).*(0:1:M-1))); %Perfect beamformer
    Beamformer_Error = Beamformer.*exp(-1i.*unifrnd(0,45,length(Beamformer),1).*pi/180);
    Codebook_Pattern = Quantize(Beamformer_Error,calibrate);   %Quantized beamformer, i.e., in format 012301230123
    BeamformerQuant = 1/sqrt(M)*exp(1j*Codebook_Pattern*pi/2); %Convert quantized format to steering vector
    Codebook(nthBeam,:) = Codebook_Pattern; %save the codebook
    
    %% Go through all channel to plot the beam patterns
    Gain_Perfect = zeros(1,length(Theta_r_Channel));
    Gain_Quan = zeros(1,length(Theta_r_Channel));
    for ind_Channel = 1:length(Theta_r_Channel)
        Theta_r_Channel_Realization = Theta_r_Channel(ind_Channel);
        Channel = transpose(exp(-1i*2*pi*d/lambda*sin(Theta_r_Channel_Realization).*(0:1:M-1)));
        Gain_Perfect(ind_Channel) = abs(Beamformer'*Channel);
        Gain_Quan(ind_Channel) = abs(BeamformerQuant'*Channel);
    end
    figure
    plot(rad2deg(Theta_r_Channel),10*log10(Gain_Perfect),'r-'); grid on;
    hold on;
    plot(rad2deg(Theta_r_Channel),10*log10(Gain_Quan),'b-'); grid on;
    hold on;
end
