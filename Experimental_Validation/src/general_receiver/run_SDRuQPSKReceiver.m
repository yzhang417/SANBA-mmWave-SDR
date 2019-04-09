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
% This function runs the Rx USRP and feed the digital data stream 
% into the Rx USRP. Also, the necessary phased array configurations are
% performed, which may also involve the control of the Tx phased array
% configuration. For different objectives, different orders of phased array
% configuration are performed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% prmQPSKReceiver: a structure array that groups the related 
% parameters on the QPSK receiver including frame format, modulation 
% scheme and etc. With the generated structure array, the baseband post 
% signal processing is set.
% prmPAControl: a structure array that groups the related parameters
% on phased arrays including available codebooks, receiving power gain,
% antenna switching on/off settings and etc. With prmPAControl, the RF 
% configuration is set.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function run_SDRuQPSKReceiver(prmQPSKReceiver, prmPAControl)
%#codegen
global Raw_corruptSignal
persistent radio

if isempty(radio)
%% Create and configure the SDRu System object. Set the SerialNum for N210
      radio = comm.SDRuReceiver(...
        'Platform',                 prmQPSKReceiver.Platform, ...
        'IPAddress',                prmQPSKReceiver.Address, ...
        'CenterFrequency',          prmQPSKReceiver.USRPCenterFrequency, ...
        'Gain',                     prmQPSKReceiver.USRPGain, ...
        'DecimationFactor',         prmQPSKReceiver.USRPDecimationFactor, ...
        'SamplesPerFrame',          prmQPSKReceiver.USRPFrameLength, ...
        'OutputDataType',           'double', ...
        'LocalOscillatorOffset',    prmQPSKReceiver.USRPLocalOscillatorOffset, ...
        'TransportDataType',        prmQPSKReceiver.USRPTransportDataType, ...
        'EnableBurstMode',          prmQPSKReceiver.USRPEnableBurstMode, ...
        'NumFramesInBurst',         prmQPSKReceiver.USRPNumFramesInBurst, ...
        'ClockSource',              prmQPSKReceiver.USRPClockSource, ...
        'PPSSource',                prmQPSKReceiver.USRPPPSSource);
end

%% Get Information about the USRP
coder.extrinsic('disp');
if coder.target('MATLAB')
    radioinfo = info(radio);
    disp(radioinfo);
end

%% Open serial port and Start Receiving Process
Holding_Time_Open_Serial = prmPAControl.Holding_Time_Open_Serial;
Holding_Time_For_Writing_Buffer = prmPAControl.Holding_Time_For_Writing_Buffer;
Holding_Time_Turn_On_Antenna = prmPAControl.Holding_Time_Turn_On_Antenna;
Holding_Time_Change_Phase_Shifter_State = prmPAControl.Holding_Time_Change_Phase_Shifter_State;

if ~coder.target('MATLAB')         
    %% Program_ID = 9 Element Gain Calibration
    if prmPAControl.Program_ID == 9 || prmPAControl.Program_ID == 9.1 || prmPAControl.Program_ID == 9.2
        disp('-------------------------------------------------------------------------');
        disp('---------- Start Element Gain/Phase Error Calibration Program -----------');
        disp('-------------------------------------------------------------------------');      
        N_State = prmPAControl.Number_State_To_Test;
        for state = 1:1:N_State
            NthFrame = 1;
            len=uint32(0);
            while NthFrame <= prmPAControl.Number_of_Frames_Per_Test_Rx_Buffered
                %Keep accessing the SDRu System object output until it is valid       
                while len <= 0
                    [corruptSignal, len] = step(radio);
                    Raw_corruptSignal(:, NthFrame, state) = corruptSignal;
                end
                len=uint32(0);
                NthFrame = NthFrame + 1;
            end
        end
        disp('------------------------------------------------------------------------');
        disp('---------- Stop Element Gain/Phase Error Calibration Program -----------');
        disp('------------------------------------------------------------------------');
    %% Program_ID = 1 Power Level test
    elseif prmPAControl.Program_ID == 1
        disp('-----------------------------------------------------');
        disp('---------- Start Power Level Test Program -----------');
        disp('-----------------------------------------------------');
        %Open serial
        coder.ceval('OpenSerialPortRx', prmPAControl.SerialPortRx);
        coder.ceval('Operation_Interval_Time',Holding_Time_Open_Serial);
        %Open 1-st receiver antenna
        coder.ceval('WriteCommandRx', prmPAControl.rxon(1,:), prmPAControl.SerialPortRx);
        coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
        coder.ceval('ReadCommandRxPrint');
        None_Zero_List_Index_rfvga_to_test = zeros(1,prmPAControl.Number_rfvga_to_test);
        ind2 = 1;
        for ind1 = 1:1:length(prmPAControl.List_Index_rfvga_to_test)
            if prmPAControl.List_Index_rfvga_to_test(ind1) == 1
                None_Zero_List_Index_rfvga_to_test(ind2) = ind1; 
                ind2 = ind2 + 1;
            end
        end
        None_Zero_List_Index_rxvga_to_test = zeros(1,prmPAControl.Number_rxvga_to_test);
        ind2 = 1;
        for ind1 = 1:1:length(prmPAControl.List_Index_rxvga_to_test)
            if prmPAControl.List_Index_rxvga_to_test(ind1) == 1
                None_Zero_List_Index_rxvga_to_test(ind2) = ind1; 
                ind2 = ind2 + 1;
            end
        end
        disp('********** Start Receiving Signal **********');
        state = 1;
        flag = 0;
        for i = 1:1:length(prmPAControl.List_Index_rfvga_to_test) 
            if i > prmPAControl.Number_rfvga_to_test
                break;
            end
            Current_Index_rfvga = None_Zero_List_Index_rfvga_to_test(i);
            coder.ceval('WriteCommandRx', prmPAControl.rfvga(Current_Index_rfvga,:), prmPAControl.SerialPortRx);
            coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
            coder.ceval('ReadCommandRxPrint');
            for j = 1:1:length(prmPAControl.List_Index_rxvga_to_test)
                if j > prmPAControl.Number_rxvga_to_test
                    break;
                end
                Current_Index_rxvga = None_Zero_List_Index_rxvga_to_test(j);
                coder.ceval('WriteCommandRx', prmPAControl.rxvga(Current_Index_rxvga,:), prmPAControl.SerialPortRx);
                coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
                coder.ceval('ReadCommandRxPrint');
                NthFrame = 1;
                len=uint32(0);    
                while NthFrame <= prmPAControl.Number_of_Frames_Per_Test_Rx_Buffered
                    %Keep accessing the SDRu System object output until it is valid       
                    while len <= 0
                        [corruptSignal, len] = step(radio);
                        Raw_corruptSignal(:, NthFrame, state) = corruptSignal;
                        %Raw_corruptSignal(:, NthFrame, state) = Raw_corruptSignal(:, NthFrame, state);
                        %len = uint32(1);
                    end
                    len = uint32(0);
                    NthFrame = NthFrame + 1;
                end
                state = state + 1;
                if state > prmPAControl.Maximum_Number_State_To_Test
                    flag = 1;
                    break;
                end
            end
            if flag == 1
                break;
            end
        end
        coder.ceval('CloseSerialPortRx');
        disp('Number of state recorded:');
        disp(state-1);
        disp('********** Stop Receiving Signal **********');
        disp('----------------------------------------------------');
        disp('---------- Stop Power Level Test Program -----------');   
        disp('----------------------------------------------------');

        
    %% Program_ID = 2 || 4 || 6 : Transmitter Phase Calibration/Transmit Beam Training
    elseif prmPAControl.Program_ID == 2 || prmPAControl.Program_ID == 4 || prmPAControl.Program_ID == 6
        if prmPAControl.Program_ID == 2
            disp('-------------------------------------------------------------------------');
            disp('---------- Start Transmitter Phase Shifter Calibration Program ----------');
            disp('-------------------------------------------------------------------------');
        elseif prmPAControl.Program_ID == 4
            disp('-------------------------------------------------------------');
            disp('---------- Start Transmitter Beam Training Program ----------');
            disp('-------------------------------------------------------------');
        else %prmPAControl.Program_ID == 6
            disp('-------------------------------------------------------------');
            disp('---------- Start Overrall Beam Training Program -------------');
            disp('-------------------------------------------------------------');           
        end
        %Open serial
        coder.ceval('OpenSerialPortRx', prmPAControl.SerialPortRx);
        coder.ceval('Operation_Interval_Time',Holding_Time_Open_Serial);
        %Open receiver antenna
        if prmPAControl.Program_ID ~= 6
            coder.ceval('WriteCommandRx', prmPAControl.rxon_single(1,:), prmPAControl.SerialPortRx);
            coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
            coder.ceval('ReadCommandRxPrint');
        else
            coder.ceval('WriteCommandRx', prmPAControl.rxonall, prmPAControl.SerialPortRx);
            coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
            coder.ceval('ReadCommandRxPrint');    
            %Apply the rx beam to train
            coder.ceval('WriteCommandRx', prmPAControl.CodeBook_Rx(prmPAControl.nthRxBeam,:), prmPAControl.SerialPortRx);
            coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
            coder.ceval('ReadCommandRxPrint'); 
        end
        %Set rxvga
        coder.ceval('WriteCommandRx',prmPAControl.rxvga(prmPAControl.rx_rxvga_to_use+1,:),prmPAControl.SerialPortRx);
        coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
        coder.ceval('ReadCommandRxPrint');
        %Set rfvga
        coder.ceval('WriteCommandRx',prmPAControl.rfvga(prmPAControl.rx_rfvga_to_use+1,:),prmPAControl.SerialPortRx);
        coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
        coder.ceval('ReadCommandRxPrint');
        coder.ceval('CloseSerialPortRx');  
        state = 1;
        if prmPAControl.Program_ID == 2
            N_State = prmPAControl.Maximum_Number_State_To_Test;
        else
            N_State = prmPAControl.Number_CodeBook_Tx_to_train;
        end
        disp('********** Start Receiving Signal **********');
        for i = 1:1:N_State
            NthFrame = 1;
            len=uint32(0);    
            while NthFrame <= prmPAControl.Number_of_Frames_Per_Test_Rx_Buffered
                %Keep accessing the SDRu System object output until it is valid       
                while len <= 0
                    [corruptSignal, len] = step(radio);
                    Raw_corruptSignal(:, NthFrame, state) = corruptSignal;
                    %Raw_corruptSignal(:, NthFrame, state) = Raw_corruptSignal(:, NthFrame, state);
                    %len = uint32(1);
                end
                len=uint32(0);
                NthFrame = NthFrame + 1;
            end
            state = state + 1;
        end
        disp('Number of state recorded:');
        disp(state-1);
        disp('********** Stop Receiving Signal **********');
        if prmPAControl.Program_ID == 2
            disp('------------------------------------------------------------------------');
            disp('---------- Stop Transmitter Phase Shifter Calibration Program ----------');
            disp('------------------------------------------------------------------------');
        elseif prmPAControl.Program_ID == 4
            disp('------------------------------------------------------------');
            disp('---------- Stop Transmitter Beam Training Program ----------');
            disp('------------------------------------------------------------');
        else %prmPAControl.Program_ID == 6
            disp('------------------------------------------------------------');
            disp('---------- Stop Overrall Beam Training Program -------------');
            disp('------------------------------------------------------------');               
        end

        
    %% Program_ID = 3: Receiver Phase Calibration
    elseif prmPAControl.Program_ID == 3
        disp('-----------------------------------------------------------------------');
        disp('---------- Start Receiver Phase Shifter Calibration Program -----------');
        disp('-----------------------------------------------------------------------');
        %Open serial
        coder.ceval('OpenSerialPortRx',prmPAControl.SerialPortRx);
        coder.ceval('Operation_Interval_Time',Holding_Time_Open_Serial);
        %Turn on 1-st antenna
        %coder.ceval('WriteCommandRx',prmPAControl.rxoffall,prmPAControl.SerialPortRx);
        coder.ceval('WriteCommandRx',prmPAControl.rxon(1,:),prmPAControl.SerialPortRx);
        coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
        coder.ceval('ReadCommandRxPrint');
        %Set rxvga
        coder.ceval('WriteCommandRx',prmPAControl.rxvga(prmPAControl.rx_rxvga_to_use+1,:),prmPAControl.SerialPortRx);
        coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
        coder.ceval('ReadCommandRxPrint');
        %Set rfvga
        coder.ceval('WriteCommandRx',prmPAControl.rfvga(prmPAControl.rx_rfvga_to_use+1,:),prmPAControl.SerialPortRx);
        coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
        coder.ceval('ReadCommandRxPrint');
        None_Zero_List_Index_antenna_to_calibrate = zeros(1,prmPAControl.Number_antenna_to_calibrate);
        ind2 = 1;
        for ind1 = 1:1:length(prmPAControl.List_Index_antenna_to_calibrate)
            if prmPAControl.List_Index_antenna_to_calibrate(ind1) == 1
                None_Zero_List_Index_antenna_to_calibrate(ind2) = ind1; 
                ind2 = ind2 + 1;
            end
        end
        disp('********** Start Receiving Signal **********');
        state = 1;
        flag = 0;
        for i = 1:1:length(prmPAControl.List_Index_antenna_to_calibrate) 
            if i > prmPAControl.Number_antenna_to_calibrate
                break;
            end
            Current_Index_antenna = None_Zero_List_Index_antenna_to_calibrate(i);
            coder.ceval('WriteCommandRx', prmPAControl.rxon(Current_Index_antenna,:), prmPAControl.SerialPortRx);
            %coder.ceval('WriteCommandRx', prmPAControl.rxon_single(Current_Index_antenna,:), prmPAControl.SerialPortRx);
            coder.ceval('Operation_Interval_Time',Holding_Time_Change_Phase_Shifter_State);
            %coder.ceval('ReadCommandRxPrint');
            for j=1:1:prmPAControl.Number_Phase_Shifter_State
                coder.ceval('WriteCommandRx', prmPAControl.phvec((Current_Index_antenna-2)*4+j,:), prmPAControl.SerialPortRx);
                coder.ceval('Operation_Interval_Time',Holding_Time_Change_Phase_Shifter_State);
                %coder.ceval('ReadCommandRxPrint');
                NthFrame = 1;
                len=uint32(0);    
                while NthFrame <= prmPAControl.Number_of_Frames_Per_Test_Rx_Buffered
                    %Keep accessing the SDRu System object output until it is valid       
                    while len <= 0
                        [corruptSignal, len] = step(radio);
                        Raw_corruptSignal(:, NthFrame, state) = corruptSignal;
                        %Raw_corruptSignal(:, NthFrame, state) = Raw_corruptSignal(:, NthFrame, state);
                        %len = uint32(1);
                    end
                    len=uint32(0);
                    NthFrame = NthFrame + 1;
                end
                state = state + 1;
                if state > prmPAControl.Maximum_Number_State_To_Test
                    flag = 1;
                    break;                
                end
            end
            if flag == 1
                break;
            end
        end
        coder.ceval('CloseSerialPortRx');
        disp('Number of state recorded:');
        disp(state-1);
        disp('********** Stop Receiving Signal **********');
        disp('----------------------------------------------------------------------');
        disp('---------- Stop Receiver Phase Shifter Calibration Program -----------');    
        disp('----------------------------------------------------------------------');

        
    %% Program_ID = 5: Receiver Beam Training
    elseif prmPAControl.Program_ID == 5 || prmPAControl.Program_ID == 8 || prmPAControl.Program_ID == 8.1
        if prmPAControl.Program_ID == 5
            disp('-------------------------------------------------------------');
            disp('---------- Start Receiver Beam Training Program -------------');
            disp('-------------------------------------------------------------');
        elseif prmPAControl.Program_ID == 8
            disp('-------------------------------------------------------------');
            disp('---------- Start Channel Estimation Program -----------------');
            disp('-------------------------------------------------------------'); 
        else %
            disp('-------------------------------------------------------------');
            disp('---------- Start Overall Beam Training Program (wired)-------');
            disp('-------------------------------------------------------------');            
        end        
        %Open Serial
        coder.ceval('OpenSerialPortRx', prmPAControl.SerialPortRx);
        coder.ceval('Operation_Interval_Time',Holding_Time_Open_Serial);
        %Turn on all antennas
        coder.ceval('WriteCommandRx', prmPAControl.rxonall, prmPAControl.SerialPortRx);
        coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
        coder.ceval('ReadCommandRxPrint');
        %Reset phase vector to "0x0000000"
        coder.ceval('WriteCommandRx', 'phvec 0', prmPAControl.SerialPortRx);
        coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
        coder.ceval('ReadCommandRxPrint');
        %Set rxvga
        coder.ceval('WriteCommandRx',prmPAControl.rxvga(prmPAControl.rx_rxvga_to_use+1,:),prmPAControl.SerialPortRx);
        coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
        coder.ceval('ReadCommandRxPrint');
        %Set rfvga
        coder.ceval('WriteCommandRx',prmPAControl.rfvga(prmPAControl.rx_rfvga_to_use+1,:),prmPAControl.SerialPortRx);
        coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
        coder.ceval('ReadCommandRxPrint');
        
        None_Zero_List_Index_CodeBook_Rx_to_train = zeros(1,prmPAControl.Number_CodeBook_Rx_to_train);
        ind2 = 1;
        for ind1 = 1:1:length(prmPAControl.List_Index_CodeBook_Rx_to_train)
            if prmPAControl.List_Index_CodeBook_Rx_to_train(ind1) == 1
                None_Zero_List_Index_CodeBook_Rx_to_train(ind2) = ind1; 
                ind2 = ind2 + 1;
            end
        end
        disp('********** Start Receiving Signal **********');
        state = 1; 
        for i = 1:1:prmPAControl.Number_CodeBook_Rx  
            if i > prmPAControl.Number_CodeBook_Rx_to_train
                break;
            end
            Current_Index_CodeBook_Rx = None_Zero_List_Index_CodeBook_Rx_to_train(i);
            coder.ceval('WriteCommandRx', prmPAControl.CodeBook_Rx(Current_Index_CodeBook_Rx,:), prmPAControl.SerialPortRx);
            coder.ceval('Operation_Interval_Time',Holding_Time_Change_Phase_Shifter_State);
            coder.ceval('ReadCommandRxPrint');
            NthFrame = 1;
            len=uint32(0);    
            while NthFrame <= prmPAControl.Number_of_Frames_Per_Test_Rx_Buffered
                %Keep accessing the SDRu System object output until it is valid       
                while len <= 0
                    [corruptSignal, len] = step(radio);
                    Raw_corruptSignal(:, NthFrame, state) = corruptSignal;
                    %Raw_corruptSignal(:, NthFrame, state) = Raw_corruptSignal(:, NthFrame, state);
                    %len = uint32(1);
                end
                len=uint32(0);
                NthFrame = NthFrame + 1;
            end
            state = state + 1;
            if state > prmPAControl.Maximum_Number_State_To_Test
                break;                
            end
        end
        coder.ceval('CloseSerialPortRx');
        disp('Number of state recorded:');
        disp(state-1);
        disp('********** Stop Receiving Signal **********');
        if prmPAControl.Program_ID == 5
            disp('----------------------------------------------------------');
            disp('---------- Stop Receiver Beam Training Program -----------');
            disp('----------------------------------------------------------');
        elseif prmPAControl.Program_ID == 8
            disp('----------------------------------------------------------');
            disp('---------- Stop Channel Estimation Program ---------------');
            disp('----------------------------------------------------------');
        else
            disp('----------------------------------------------------------');
            disp('---------- Stop Overall Beam Training Program (wired) ----');
            disp('----------------------------------------------------------'); 
        end
        
    %% Program_ID = 2.1: Transmitter Phase Calibration (wired cable version)
    elseif prmPAControl.Program_ID == 2.1
        disp('----------------------------------------------------------------------------------');
        disp('---------- Start Transmitter Phase Shifter Calibration Program (Wired) -----------');
        disp('----------------------------------------------------------------------------------');
        %Open serial
        coder.ceval('OpenSerialPortTx',prmPAControl.SerialPortTx);
        coder.ceval('Operation_Interval_Time',Holding_Time_Open_Serial);
        %Turn on 1-st transmitter antenna
        coder.ceval('WriteCommandTx',prmPAControl.txon(1,:),prmPAControl.SerialPortTx);
        coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
        coder.ceval('ReadCommandTxPrint');
        %Reset phase vector to "0x0000000"
        coder.ceval('WriteCommandTx', 'phvec 0', prmPAControl.SerialPortTx);
        coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
        coder.ceval('ReadCommandTxPrint');       
        %%
        None_Zero_List_Index_antenna_to_calibrate = zeros(1,prmPAControl.Number_antenna_to_calibrate);
        ind2 = 1;
        for ind1 = 1:1:length(prmPAControl.List_Index_antenna_to_calibrate)
            if prmPAControl.List_Index_antenna_to_calibrate(ind1) == 1
                None_Zero_List_Index_antenna_to_calibrate(ind2) = ind1; 
                ind2 = ind2 + 1;
            end
        end
        disp('********** Start Receiving Signal **********');
        state = 1;
        flag = 0;
        for i = 1:1:length(prmPAControl.List_Index_antenna_to_calibrate) 
            if i > prmPAControl.Number_antenna_to_calibrate
                break;
            end
            Current_Index_antenna = None_Zero_List_Index_antenna_to_calibrate(i);
            coder.ceval('WriteCommandTx', prmPAControl.txon(Current_Index_antenna,:), prmPAControl.SerialPortTx);
            coder.ceval('Operation_Interval_Time',Holding_Time_Change_Phase_Shifter_State);
            %coder.ceval('ReadCommandTxPrint');
            for j=1:1:prmPAControl.Number_Phase_Shifter_State
                coder.ceval('WriteCommandTx', prmPAControl.phvec((Current_Index_antenna-2)*4+j,:), prmPAControl.SerialPortTx);
                coder.ceval('Operation_Interval_Time',Holding_Time_Change_Phase_Shifter_State);
                %coder.ceval('ReadCommandTxPrint');
                NthFrame = 1;
                len=uint32(0);    
                while NthFrame <= prmPAControl.Number_of_Frames_Per_Test_Rx_Buffered
                    %Keep accessing the SDRu System object output until it is valid       
                    while len <= 0
                        [corruptSignal, len] = step(radio);
                        Raw_corruptSignal(:, NthFrame, state) = corruptSignal;
                        %Raw_corruptSignal(:, NthFrame, state) = Raw_corruptSignal(:, NthFrame, state);
                        %len = uint32(1);
                    end
                    len=uint32(0);
                    NthFrame = NthFrame + 1;
                end
                state = state + 1;
                if state > prmPAControl.Maximum_Number_State_To_Test
                    flag = 1;
                    break;                
                end
            end
            if flag == 1
                break;
            end
        end
        coder.ceval('CloseSerialPortTx');
        disp('Number of state recorded:');
        disp(state-1);
        disp('********** Stop Receiving Signal **********');
        disp('------------------------------------------------------------------------------');
        disp('---------- Stop Transmitter Phase Shifter Calibration Program (Wired) --------');    
        disp('------------------------------------------------------------------------------');

        
    %% Program_ID = 4.1: Transmitter Beam Training (wired cable version)
    elseif prmPAControl.Program_ID == 4.1
        disp('---------------------------------------------------------------------');
        disp('---------- Start Transmitter Beam Training Program (Wired) -------------');
        disp('---------------------------------------------------------------------');
        %Open serial
        coder.ceval('OpenSerialPortTx',prmPAControl.SerialPortTx);
        coder.ceval('Operation_Interval_Time',Holding_Time_Open_Serial);
        %Turn on 1-st transmitter antenna
        coder.ceval('WriteCommandTx',prmPAControl.txonall,prmPAControl.SerialPortTx);
        coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
        coder.ceval('ReadCommandTxPrint');
        %Reset phase vector to "0x0000000"
        coder.ceval('WriteCommandTx', 'phvec 0', prmPAControl.SerialPortTx);
        coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
        coder.ceval('ReadCommandTxPrint');
        %%
        None_Zero_List_Index_CodeBook_Tx_to_train = zeros(1,prmPAControl.Number_CodeBook_Tx_to_train);
        ind2 = 1;
        for ind1 = 1:1:length(prmPAControl.List_Index_CodeBook_Tx_to_train)
            if prmPAControl.List_Index_CodeBook_Tx_to_train(ind1) == 1
                None_Zero_List_Index_CodeBook_Tx_to_train(ind2) = ind1; 
                ind2 = ind2 + 1;
            end
        end
        disp('********** Start Receiving Signal **********');
        state = 1; 
        for i = 1:1:prmPAControl.Number_CodeBook_Tx  
            if i > prmPAControl.Number_CodeBook_Tx_to_train
                break;
            end
            Current_Index_CodeBook_Tx = None_Zero_List_Index_CodeBook_Tx_to_train(i);
            coder.ceval('WriteCommandTx', prmPAControl.CodeBook_Tx(Current_Index_CodeBook_Tx,:), prmPAControl.SerialPortTx);
            coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
            coder.ceval('ReadCommandTxPrint');
            NthFrame = 1;
            len=uint32(0);    
            while NthFrame <= prmPAControl.Number_of_Frames_Per_Test_Rx_Buffered
                %Keep accessing the SDRu System object output until it is valid       
                while len <= 0
                    [corruptSignal, len] = step(radio);
                    Raw_corruptSignal(:, NthFrame, state) = corruptSignal;
                    %Raw_corruptSignal(:, NthFrame, state) = Raw_corruptSignal(:, NthFrame, state);
                    %len = uint32(1);
                end
                len=uint32(0);
                NthFrame = NthFrame + 1;
            end
            state = state + 1;
            if state > prmPAControl.Maximum_Number_State_To_Test
                break;                
            end
        end
        coder.ceval('CloseSerialPortTx');
        disp('Number of state recorded:');
        disp(state-1);
        disp('********** Stop Receiving Signal **********');
        disp('--------------------------------------------------------------------');
        disp('---------- Stop Transmitter Beam Training Program (Wired)-----------');
        disp('--------------------------------------------------------------------');
   
    %% Program_ID = 7
    elseif prmPAControl.Program_ID == 7    
        disp('---------------------------------------------------------------------');
        disp('---------- Start Showing Beam Patter Program -------------');
        disp('---------------------------------------------------------------------');  
        coder.ceval('WriteCommandRx', prmPAControl.CodeBook_Rx(1,:), prmPAControl.SerialPortRx);
        coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
        coder.ceval('ReadCommandRxPrint');
        %Set rxvga
        coder.ceval('WriteCommandRx',prmPAControl.rxvga(prmPAControl.rx_rxvga_to_use+1,:),prmPAControl.SerialPortRx);
        coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
        coder.ceval('ReadCommandRxPrint');
        %Set rfvga
        coder.ceval('WriteCommandRx',prmPAControl.rfvga(prmPAControl.rx_rfvga_to_use+1,:),prmPAControl.SerialPortRx);
        coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
        coder.ceval('ReadCommandRxPrint');
             
%         scope = dsp.TimeScope('NumInputPorts',1,...
%                               'TimeSpan',prmQPSKReceiver.FrameTime,...
%                               'AxesScaling','Auto',...
%                               'YLimits', [1 15],...
%                               'TimeSpanOverrunAction', 'Scroll',...
%                               'ReduceUpdates',true);
        len=uint32(0);        
        ite = 1;
        corruptSignal = 0;
        while 1 >= 1
            while len <= 0 
                [corruptSignal, len] = step(radio);
            end
            ite = ite + 1;
            if ite == 300
                disp(rms(corruptSignal));
                ite = 1;
            end             
            len=uint32(0);
        end
            
    %% Program_ID = 1: Power Level test
    else %prmPAControl.Program_ID == 0
        disp('---------------------------------------------------------------------------------------');
        disp('---------- No Operation Is Required for Receiver in Power Calibration Program----------');
        disp('---------------------------------------------------------------------------------------');
    end
end


%% release the radio
%release(radio);




