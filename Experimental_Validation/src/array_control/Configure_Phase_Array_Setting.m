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
% function description:
% This function configures the phased arrays during the calibration data 
% collection process. It is used for Program_ID = 9.1 or 9.2. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input argument:
% prmPAControl: a structure array that groups the related parameters
% on phased arrays including available codebooks, receiving power gain,
% antenna switching on/off settings and etc. With prmPAControl, the RF 
% configuration is set.
% Fixed_Antenna_ID: the index of antenna fixed to transmit (or receive) 
% signal to (or from) the antenna to be calibrated (during the
% calibration data collection process).
% Fixed_Antenna_Phase: the index of phase state of the previously defined 
% fixed antenna.
% Ant_Ref_Ind: the index of antenna fixed at the calibration side for 
% reference (during the calibration data collection process).
% Ant_To_Be_Calibrated_Ind: the index of antenna to be calibrated.
% Ant_To_Be_Calibrated_Phase: the index of the phase state to be calibrated 
% of the antenna to be calibrated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function  Configure_Phase_Array_Setting(prmPAControl,...
                                        Fixed_Antenna_ID,...
                                        Fixed_Antenna_Phase,...
                                        Ant_Ref_Ind,...
                                        Ant_Ref_Ind_Phase,...
                                        Ant_To_Be_Calibrated_Ind,...
                                        Ant_To_Be_Calibrated_Phase) 
    %% Fixed side
    % Configure side
    if prmPAControl.Program_ID == 9.1 
       Antenna_Phase_Setting_Fixed.Configure_Side = 'Rx';
    else
       Antenna_Phase_Setting_Fixed.Configure_Side = 'Tx';
    end
    
    % Antenna setting
    antenna_array_setting = zeros(1,12);
    antenna_array_setting(Fixed_Antenna_ID) = 1;
    Antenna_Phase_Setting_Fixed.antenna_array_setting = antenna_array_setting;

    % Phase setting
    phase_array_setting = ones(1,12);
    phase_array_setting(Fixed_Antenna_ID) = Fixed_Antenna_Phase;
    Antenna_Phase_Setting_Fixed.phase_array_setting = phase_array_setting;

    % Set phase array     
    clear Control_Phase_Array_External_Phase_Error_Calibration_mex %#ok<UNRCH>
    Control_Phase_Array_External_Phase_Error_Calibration_mex(prmPAControl,Antenna_Phase_Setting_Fixed);
    
    %% Calibrated Side
    % Configure Side
    if prmPAControl.Program_ID == 9.1 
        Antenna_Phase_Setting_Cal.Configure_Side = 'Tx';
    else
        Antenna_Phase_Setting_Cal.Configure_Side = 'Rx';
    end
    
    % Antenna setting for the side to be calibrated
    antenna_array_setting = zeros(1,12);
    antenna_array_setting([Ant_Ref_Ind, Ant_To_Be_Calibrated_Ind]) = 1;
    Antenna_Phase_Setting_Cal.antenna_array_setting = antenna_array_setting;

    % Phase setting
    phase_array_setting = ones(1,12);
    phase_array_setting([Ant_Ref_Ind, Ant_To_Be_Calibrated_Ind]) = [Ant_Ref_Ind_Phase, Ant_To_Be_Calibrated_Phase];
    Antenna_Phase_Setting_Cal.phase_array_setting = phase_array_setting;

    % Set phase array            
    clear Control_Phase_Array_External_Phase_Error_Calibration_mex %#ok<UNRCH>
    Control_Phase_Array_External_Phase_Error_Calibration_mex(prmPAControl,Antenna_Phase_Setting_Cal);               
end

