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
% This script loads and demonstrates the data related to the radiation gain 
% of each array element. In particular, the averaged radiation gain of each
% antenna element is stored in Element_Gain.mat in the folder 
% [Experimental_Validation/data/cal_result] alongside with the 
% corresponding noise floor is stored in Noise_Gain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization
clc; clear all; close all

%% Add path
folder_name = ...
[
"/Experimental_Validation";
];
for i=1:length(folder_name)
    Lib_path = char(strcat(pwd,folder_name(i)));
    addpath(genpath(Lib_path));
end


%% Transfer corrupted raw data to element gain
Activate_Tx_ID_range = 1:1:12;
Activate_Tx_Phase_ID_range = 1:1:4;
Activate_Rx_ID_range = 1:1:12;
Activate_Rx_Phase_ID_range = 1:1:4;
Element_Gain = zeros(length(Activate_Tx_ID_range),...
                     length(Activate_Tx_Phase_ID_range),...
                     length(Activate_Rx_ID_range),...
                     length(Activate_Rx_Phase_ID_range));  
Noise_Gain = Element_Gain;
remove_boundary = 0.2;
Fixed_Side_Antenna_ID_Range = [5];
Fixed_Side_Antenna_Phase_Range = [1];


%% Element gain collection
for Activate_Tx_ID = Activate_Tx_ID_range
       for Activate_Tx_Phase_ID = Activate_Tx_Phase_ID_range
           for Activate_Rx_ID = Activate_Rx_ID_range
               for Activate_Rx_Phase_ID = Activate_Rx_Phase_ID_range
                   if (ismember(Activate_Tx_ID,Fixed_Side_Antenna_ID_Range) && ismember(Activate_Tx_Phase_ID,Fixed_Side_Antenna_Phase_Range)) ...
                           || (ismember(Activate_Rx_ID,Fixed_Side_Antenna_ID_Range) && ismember(Activate_Rx_Phase_ID,Fixed_Side_Antenna_Phase_Range) )
                       try
                           load(['Raw_corruptSignal_TX' num2str(Activate_Tx_ID) '_'...
                                                num2str(Activate_Tx_Phase_ID) '_Rx' ...
                                                num2str(Activate_Rx_ID) '_' ...
                                                num2str(Activate_Rx_Phase_ID) '.mat'],'Raw_corruptSignal');
                       catch
                           fprintf("\nNo such files\n");
                       end
                       Raw_corruptSignal_Removed_Boundary = Raw_corruptSignal(:,:,round(44*remove_boundary):round(44*(1-remove_boundary)));
                       rms_temp = rms(rms(rms(Raw_corruptSignal_Removed_Boundary)));
                       Element_Gain(Activate_Tx_ID,Activate_Tx_Phase_ID,Activate_Rx_ID,Activate_Rx_Phase_ID) = rms_temp;
                       try
                           load(['Noise_Raw_corruptSignal_TX' num2str(Activate_Tx_ID) '_'...
                                                num2str(Activate_Tx_Phase_ID) '_Rx' ...
                                                num2str(Activate_Rx_ID) '_' ...
                                                num2str(Activate_Rx_Phase_ID) '.mat'],'Raw_corruptSignal');
                       catch
                           fprintf("\nNo such files\n");
                       end
                       Noise_Raw_corruptSignal_Removed_Boundary = Raw_corruptSignal(:,:,round(44*remove_boundary):round(44*(1-remove_boundary)));
                       Noise_rms_temp = rms(rms(rms(Noise_Raw_corruptSignal_Removed_Boundary)));
                       Noise_Gain(Activate_Tx_ID,Activate_Tx_Phase_ID,Activate_Rx_ID,Activate_Rx_Phase_ID) = Noise_rms_temp;
                   end
               end
           end
       end
end
Element_Gain_With_Noise = Element_Gain;
Element_Gain = sqrt(Element_Gain.^2 - Noise_Gain.^2);
cal_data_path = char(strcat(pwd,'/Experimental_Validation/data/cal_result/'));
save(char(strcat(cal_data_path,'/Element_Gain.mat')),'Element_Gain');
save(char(strcat(cal_data_path,'/Noise_Gain.mat')),'Noise_Gain');


%% Element gain saving and showing
load(char(strcat(cal_data_path,'/Element_Gain.mat')),'Element_Gain');
Element_Gain_With_Noise = Element_Gain;
load(char(strcat(cal_data_path,'/Noise_Gain.mat')),'Noise_Gain');
Fixed_Side_Antenna_ID_Range = [5];
Fixed_Side_Antenna_Phase_Range = [1];

for Fixed_Side_Antenna_ID = Fixed_Side_Antenna_ID_Range
    for Fixed_Side_Antenna_Phase = Fixed_Side_Antenna_Phase_Range
        Element_Gain_Temp_Tx = Element_Gain(:,:,Fixed_Side_Antenna_ID,Fixed_Side_Antenna_Phase);
        figure
        subplot(2,1,1)
        bar(Element_Gain_Temp_Tx);grid on;
        xlabel('Index of Antenna')
        ylabel('RSSI')
        title('RSSI of Each Tx Antenna Element')
        Element_Gain_Temp_Rx = squeeze(Element_Gain(Fixed_Side_Antenna_ID,Fixed_Side_Antenna_Phase,:,:));
        subplot(2,1,2)
        bar(Element_Gain_Temp_Rx);grid on;
        xlabel('RSSI of Each Rx Antenna Element')
        ylabel('RSSI')
        title('Comparisson Between RSSI at Rx')
    end
end

Fixed_Side_Antenna_ID_Range = [5];
Fixed_Side_Antenna_Phase_Range = [1];
Element_Gain_With_Noise;
Noise_Gain;
SNR_dB= 10*log10(Element_Gain.^2./Noise_Gain.^2);
for Fixed_Side_Antenna_ID = Fixed_Side_Antenna_ID_Range
    for Fixed_Side_Antenna_Phase = Fixed_Side_Antenna_Phase_Range
        SNR_dB_Tx = SNR_dB(:,:,Fixed_Side_Antenna_ID,Fixed_Side_Antenna_Phase);
        figure
        subplot(2,1,1)
        bar(SNR_dB_Tx); grid on;
        xlabel('Index of Antenna')
        ylabel('SNR (dB)')
        title('Receiver SNR of Each Tx Antenna Element (dB)')
        SNR_dB_Rx = squeeze(SNR_dB(Fixed_Side_Antenna_ID,Fixed_Side_Antenna_Phase,:,:));
        subplot(2,1,2)
        bar(SNR_dB_Rx); grid on;
        xlabel('Index of Antenna')
        ylabel('SNR (dB)')
        title('Receiver SNR of Each Rx Antenna Element (dB)')
        %legend(['Phase State 1';'Phase State 2';'Phase State 3';'Phase State 4']);
    end
end


