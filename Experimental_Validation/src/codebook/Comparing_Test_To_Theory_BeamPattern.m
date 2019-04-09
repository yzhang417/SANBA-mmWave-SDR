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
% This function plots the comparison between the estimated real beam 
% pattern and the theoretical beam pattern.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% Theta_desired_d: desired direction in degree.
% Codebook_Pattern: codebooks corresponding to Theta_desired_d, which are
% measured manually or calculated. Please refer to
% [Experimental_Validation\src\codebook\measured_beam_pattern.m] for its
% usage example.
% calibrate: coarse calibrate result.
% ULA: structure array of ULA parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Comparing_Test_To_Theory_BeamPattern(Theta_desired_d, Codebook_Pattern, calibrate, ULA)
    M = ULA.M;
    d = ULA.d;
    lambda = ULA.lambda;
    NBeams = length(Codebook_Pattern(1,:));
    Theta_r_Channel = -pi/2:0.001:pi/2; 
    r = ceil(sqrt(NBeams));
    c = floor(sqrt(NBeams));
    figure
    for i=1:NBeams
        Theta_r_Beamformer = deg2rad(Theta_desired_d(i)); %radian of the nth beam pattern
        Beamformer = 1/sqrt(M)*transpose(exp(-1i*2*pi*d/lambda*sin(Theta_r_Beamformer).*(0:1:M-1))); %Perfect beamformer
        BeamformerQuant = 1/sqrt(M)*exp(1j*Codebook_Pattern(:,i)*pi/2)./exp(1j*calibrate'*pi/2); %Convert quantized format to steering vector

        %% Go through all channel to plot the beam patterns
        Gain_Perfect = zeros(1,length(Theta_r_Channel));
        Gain_Quan = zeros(1,length(Theta_r_Channel));
        for ind_Channel = 1:length(Theta_r_Channel)
            Theta_r_Channel_Realization = Theta_r_Channel(ind_Channel);
            Channel = transpose(exp(-1i*2*pi*d/lambda*sin(Theta_r_Channel_Realization).*(0:1:M-1)));
            Gain_Perfect(ind_Channel) = abs(Beamformer'*Channel);
            Gain_Quan(ind_Channel) = abs(BeamformerQuant'*Channel);
        end
        Gain_Perfect_dB = 10*log10(Gain_Perfect);
        Gain_Quan_dB = 10*log10(Gain_Quan);
        subplot(r,c,i);
        plot(rad2deg(Theta_r_Channel),Gain_Perfect_dB,'r-'); hold on;
        plot(rad2deg(Theta_r_Channel),Gain_Quan_dB,'b-'); hold on; grid on;
        title(num2str(Theta_desired_d(i)));
    end
end

