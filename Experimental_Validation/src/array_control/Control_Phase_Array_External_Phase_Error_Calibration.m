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
% This function is configures the phased array with the input configuration
% profile store in Antenna_Phase_Setting. It is in particular used for
% calibration process as the antenna on off and phase state would be
% frequently changed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input argument:
% prmPAControl: a structure array that groups the related parameters
% on phased arrays including available codebooks, receiving power gain,
% antenna switching on/off settings and etc. With prmPAControl, the RF 
% configuration is set.
% Antenna_Phase_Setting: the antenna/phase settings for both Tx and Rx
% phased arrays.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Control_Phase_Array_External_Phase_Error_Calibration(prmPAControl,Antenna_Phase_Setting)
%#codegen
    % Open serial port and Start Receiving Process
    Holding_Time_Open_Serial = prmPAControl.Holding_Time_Open_Serial;
    Holding_Time_For_Writing_Buffer = prmPAControl.Holding_Time_For_Writing_Buffer;
    Holding_Time_Turn_On_Antenna = prmPAControl.Holding_Time_Turn_On_Antenna;
    Holding_Time_Change_Phase_Shifter_State = prmPAControl.Holding_Time_Change_Phase_Shifter_State;
    
    % Extrinsic function
    coder.extrinsic('phvec'); 
    coder.extrinsic('antenna'); 
    command_ant = ['xxxx xxxxx'];
    command_phase = ['xxxxx xxxxxxxx'];
    
    if Antenna_Phase_Setting.Configure_Side == "Tx"
        %Open Tx serial
        coder.ceval('OpenSerialPortTx', prmPAControl.SerialPortTx);
        coder.ceval('Operation_Interval_Time',Holding_Time_Open_Serial);
        %Open the Antenna
        command_ant = antenna(Antenna_Phase_Setting.antenna_array_setting,"Tx");
        coder.ceval('WriteCommandTx', command_ant, prmPAControl.SerialPortTx);
        coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
        coder.ceval('ReadCommandTxPrint');    
        %Apply the phase
        command_phase = phvec(Antenna_Phase_Setting.phase_array_setting-1);
        coder.ceval('WriteCommandTx', command_phase, prmPAControl.SerialPortTx);
        coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
        coder.ceval('ReadCommandTxPrint');
        %Close serial
        coder.ceval('CloseSerialPortTx');
        coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
    else
        %Open Rx serial
        coder.ceval('OpenSerialPortRx', prmPAControl.SerialPortRx);
        coder.ceval('Operation_Interval_Time',Holding_Time_Open_Serial);
        %Open the Antenna
        command_ant = antenna(Antenna_Phase_Setting.antenna_array_setting,"Rx");
        coder.ceval('WriteCommandRx', command_ant, prmPAControl.SerialPortRx);
        coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
        coder.ceval('ReadCommandRxPrint');    
        %Apply the phase
        command_phase = phvec(Antenna_Phase_Setting.phase_array_setting-1);
        coder.ceval('WriteCommandRx', command_phase, prmPAControl.SerialPortRx);
        coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
        coder.ceval('ReadCommandRxPrint');
        %Close serial
        coder.ceval('CloseSerialPortRx');
        coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);             
    end
end
