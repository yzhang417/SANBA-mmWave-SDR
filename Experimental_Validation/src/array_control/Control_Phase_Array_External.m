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
% This function is only used for Program_ID = 8.1, where the overall beam
% training is performed or the measurements collection is performed for
% validating the proposed non-coherent algorithm. The reason of creating
% this function was due to that simultaneous controls of two phased arrays
% in one program is not allowed if the mex file is not cleared. By clearing
% the mex file every time before controlling the different phased array,
% this function could be abandoned actually.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input argument:
% prmPAControl: a structure array that groups the related parameters
% on phased arrays including available codebooks, receiving power gain,
% antenna switching on/off settings and etc. With prmPAControl, the RF 
% configuration is set.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Control_Phase_Array_External(prmPAControl)
%#codegen
    % Open serial port and start receiving process
    Holding_Time_Open_Serial = prmPAControl.Holding_Time_Open_Serial;
    Holding_Time_For_Writing_Buffer = prmPAControl.Holding_Time_For_Writing_Buffer;
    Holding_Time_Turn_On_Antenna = prmPAControl.Holding_Time_Turn_On_Antenna;
    Holding_Time_Change_Phase_Shifter_State = prmPAControl.Holding_Time_Change_Phase_Shifter_State;
    switch prmPAControl.Program_ID 
        case 8.1
            %Open Tx serial
            coder.ceval('OpenSerialPortTx', prmPAControl.SerialPortTx);
            coder.ceval('Operation_Interval_Time',Holding_Time_Open_Serial);
            coder.ceval('WriteCommandTx', prmPAControl.txonall, prmPAControl.SerialPortTx);
            coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
            coder.ceval('ReadCommandTxPrint');    
            %Apply the tx beam to train
            coder.ceval('WriteCommandTx', prmPAControl.CodeBook_Tx(prmPAControl.nthTxBeam,:), prmPAControl.SerialPortTx);
            coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
            coder.ceval('ReadCommandTxPrint');
            %Close serial
            coder.ceval('CloseSerialPortTx');
            coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
        case 9
            if prmPAControl.ControlTx == 1
                %Open Tx serial
                coder.ceval('OpenSerialPortTx', prmPAControl.SerialPortTx);
                coder.ceval('Operation_Interval_Time',Holding_Time_Open_Serial);
                coder.ceval('WriteCommandTx', prmPAControl.txon_single(prmPAControl.Activate_Tx_ID,:), prmPAControl.SerialPortTx);
                coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
                coder.ceval('ReadCommandTxPrint');    
                %Apply the tx phase state to train
                if prmPAControl.Activate_Tx_ID == 1
                    coder.ceval('WriteCommandTx', prmPAControl.phvec_Antenna1(prmPAControl.Activate_Tx_Phase_ID,:), prmPAControl.SerialPortTx);    
                else
                    coder.ceval('WriteCommandTx', prmPAControl.phvec((prmPAControl.Activate_Tx_ID-2)*4+prmPAControl.Activate_Tx_Phase_ID,:), prmPAControl.SerialPortTx);
                end
                coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
                coder.ceval('ReadCommandTxPrint');
                %Close serial
                coder.ceval('CloseSerialPortTx');
                coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
            else
                %Open Rx serial
                coder.ceval('OpenSerialPortRx', prmPAControl.SerialPortRx);
                coder.ceval('Operation_Interval_Time',Holding_Time_Open_Serial);
                coder.ceval('WriteCommandRx', prmPAControl.rxon_single(prmPAControl.Activate_Rx_ID,:), prmPAControl.SerialPortRx);
                coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
                coder.ceval('ReadCommandRxPrint');    
                %Apply the Rx phase state to train
                if prmPAControl.Activate_Rx_ID == 1
                    coder.ceval('WriteCommandRx', prmPAControl.phvec_Antenna1(prmPAControl.Activate_Rx_Phase_ID,:), prmPAControl.SerialPortTx);    
                else
                    coder.ceval('WriteCommandRx', prmPAControl.phvec((prmPAControl.Activate_Rx_ID-2)*4+prmPAControl.Activate_Rx_Phase_ID,:), prmPAControl.SerialPortTx);
                end
                coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
                coder.ceval('ReadCommandRxPrint');
                %Close serial
                coder.ceval('CloseSerialPortRx');
                coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);                
            end
    end
end
