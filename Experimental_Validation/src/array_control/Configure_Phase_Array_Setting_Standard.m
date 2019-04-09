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
% collection process. It is used for Program_ID = 9.1 or 9.2. It is 
% actually a generalized version of Configure_Phase_Array_Setting.m so that
% it allows more general configurations of phased arrays.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input argument:
% prmPAControl: a structure array that groups the related parameters
% on phased arrays including available codebooks, receiving power gain,
% antenna switching on/off settings and etc. With prmPAControl, the RF 
% configuration is set.
% Activate_Tx_ID: the indicator vector of the activation settings for 
% the Tx antenna elements
% Activate_Tx_Phase_ID: the indicator vector of the phase states for the Tx
% antenna elements
% Activate_Rx_ID: the indicator vector of the activation settings for 
% the Rx antenna elements
% Activate_Rx_Phase_ID: the indicator vector of the phase states for the Rx
% antenna elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Configure_Phase_Array_Setting_Standard(prmPAControl,...
                                                Activate_Tx_ID,...
                                                Activate_Tx_Phase_ID,...
                                                Activate_Rx_ID,...
                                                Activate_Rx_Phase_ID)
    %% TX Side
    % Configuring side
    Antenna_Phase_Setting_Standard.Configure_Side = 'Tx';
    
    % Antenna setting for the Tx side
    antenna_array_setting = zeros(1,12);
    antenna_array_setting(Activate_Tx_ID) = 1;
    Antenna_Phase_Setting_Standard.antenna_array_setting = antenna_array_setting;

    % Phase setting
    phase_array_setting = ones(1,12);
    phase_array_setting(Activate_Tx_ID) = Activate_Tx_Phase_ID;
    Antenna_Phase_Setting_Standard.phase_array_setting = phase_array_setting;

    % Set phase array            
    clear Control_Phase_Array_External_Phase_Error_Calibration_mex %#ok<UNRCH>
    Control_Phase_Array_External_Phase_Error_Calibration_mex(prmPAControl,Antenna_Phase_Setting_Standard);               

    %% RX side
    % Configuring side
    Antenna_Phase_Setting_Standard.Configure_Side = 'Rx';
    
    % Antenna setting for the Rx side
    antenna_array_setting = zeros(1,12);
    antenna_array_setting(Activate_Rx_ID) = 1;
    Antenna_Phase_Setting_Standard.antenna_array_setting = antenna_array_setting;

    % Phase setting
    phase_array_setting = ones(1,12);
    phase_array_setting(Activate_Rx_ID) = Activate_Rx_Phase_ID;
    Antenna_Phase_Setting_Standard.phase_array_setting = phase_array_setting;

    % Set phase array            
    clear Control_Phase_Array_External_Phase_Error_Calibration_mex %#ok<UNRCH>
    Control_Phase_Array_External_Phase_Error_Calibration_mex(prmPAControl,Antenna_Phase_Setting_Standard);               
end

