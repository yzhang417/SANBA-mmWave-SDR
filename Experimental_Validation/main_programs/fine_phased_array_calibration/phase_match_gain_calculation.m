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
% This script loads, calculates and demonstrates the data related to the 
% averaged radiation gain of each combination of two antenna elements. The 
% result is stored in REV_Gain9.1.mat, in the folder 
% [Experimental_Validation/data/cal_result], for Tx phased array PMBC and 
% in REV_Gain9.2.mat for Rx phased array PMBC.
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

%% Parameters
prmPAControl.Program_ID = 9.2;
remove_boundary = 0.2;
Ant_To_Be_Calibrated_Ind_range = 2:12;
Ant_To_Be_Calibrated_Phase_range = 1:1:4;
Ant_Ref_Ind_Range = [1];
Ant_Ref_Phase_Range = [1 2 3 4];
Fixed_Side_Antenna_ID_Range = [5];
Fixed_Side_Antenna_Phase_Range = [1];
REV_Gain = zeros(12,4,12,4,12,4);    
Noise_REV_Gain = REV_Gain;
cal_data_path = char(strcat(pwd,'/Experimental_Validation/data/cal_result/'));

%% Transfer corrupted raw data to phase match gain
for Fixed_Antenna_ID = Fixed_Side_Antenna_ID_Range
   for Fixed_Antenna_Phase = Fixed_Side_Antenna_Phase_Range
       for Ant_Ref_Ind = Ant_Ref_Ind_Range
           for Ant_Ref_Ind_Phase = Ant_Ref_Phase_Range
               for Ant_To_Be_Calibrated_Ind = Ant_To_Be_Calibrated_Ind_range            
                   for Ant_To_Be_Calibrated_Phase = Ant_To_Be_Calibrated_Phase_range
                        try
                        load(['Raw_corruptSignal_To_Cal_' num2str(Ant_To_Be_Calibrated_Ind) '_'...
                                                                        num2str(Ant_To_Be_Calibrated_Phase) '_Ref' ...
                                                                        num2str(Ant_Ref_Ind) '_' ...
                                                                        num2str(Ant_Ref_Ind_Phase) '_Fix' ...
                                                                        num2str(Fixed_Antenna_ID) '_' ...
                                                                        num2str(Fixed_Antenna_Phase) ...
                                                                        num2str(prmPAControl.Program_ID) '.mat'],'Raw_corruptSignal');                          
                        Raw_corruptSignal_Removed_Boundary = Raw_corruptSignal(:,:,round(44*remove_boundary)+1:round(44*(1-remove_boundary)));
                        rms_temp = rms(rms(rms(Raw_corruptSignal_Removed_Boundary)));
                        REV_Gain(Ant_To_Be_Calibrated_Ind,...
                                 Ant_To_Be_Calibrated_Phase,...
                                 Ant_Ref_Ind,...
                                 Ant_Ref_Ind_Phase,...
                                 Fixed_Antenna_ID,...
                                 Fixed_Antenna_Phase) = rms_temp;
                        catch
                            fprintf("\nNo such files 1\n");
                        end
                        try
                        load(['Noise_Raw_corruptSignal_To_Cal_' num2str(Ant_To_Be_Calibrated_Ind) '_'...
                                                                        num2str(Ant_To_Be_Calibrated_Phase) '_Ref' ...
                                                                        num2str(Ant_Ref_Ind) '_' ...
                                                                        num2str(Ant_Ref_Ind_Phase) '_Fix' ...
                                                                        num2str(Fixed_Antenna_ID) '_' ...
                                                                        num2str(Fixed_Antenna_Phase) ...
                                                                        num2str(prmPAControl.Program_ID) '.mat'],'Raw_corruptSignal');                                          
                        Noise_Raw_corruptSignal_Removed_Boundary = Raw_corruptSignal(:,:,round(44*remove_boundary)+1:round(44*(1-remove_boundary)));
                        Noise_rms_temp = rms(rms(rms(Noise_Raw_corruptSignal_Removed_Boundary)));
                        Noise_REV_Gain(Ant_To_Be_Calibrated_Ind,...
                        Ant_To_Be_Calibrated_Phase,...
                        Ant_Ref_Ind,...
                        Ant_Ref_Ind_Phase,...
                        Fixed_Antenna_ID,...
                        Fixed_Antenna_Phase) = Noise_rms_temp;                        
                        catch
                            fprintf("\nNo such files 2\n");
                        end
                   end
               end
            end
        end
    end
end
REV_Gain = sqrt(REV_Gain.^2 - Noise_REV_Gain.^2);
REV_Gain = real(REV_Gain);
save(char(strcat(cal_data_path,['REV_Gain' num2str(prmPAControl.Program_ID) '.mat'])),'REV_Gain');

%% Plot coarse calibration result
for Fixed_Antenna_ID = Fixed_Side_Antenna_ID_Range
   for Fixed_Antenna_Phase = Fixed_Side_Antenna_Phase_Range
       for Ant_Ref_Ind = Ant_Ref_Ind_Range
           for Ant_Ref_Ind_Phase = Ant_Ref_Phase_Range
                REV_Cal = REV_Gain(:,:,Ant_Ref_Ind,Ant_Ref_Ind_Phase,Fixed_Antenna_ID,Fixed_Antenna_Phase);
                figure
                subplot(3,1,1);
                bar(REV_Cal)
                [v P] = max(REV_Cal');
                subplot(3,1,2);
                P(1) = Ant_Ref_Ind_Phase;
                bar(P)
                subplot(3,1,3);
                Delta = P(1) - 1;
                for n=1:1:12
                    P(n) = P(n)-Delta;
                    if P(n) <= 0
                        P(n) = P(n) + 4;
                    end
                end
                bar(P)
           end
       end
   end
end



