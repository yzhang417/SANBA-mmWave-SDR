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
% This function runs the Tx USRP and feed the digital data stream into the
% Tx USRP. For different objectives, different operations of phased array
% configuration will be performed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% prmQPSKTransmitter: a structure array that groups the related parameters
% on the QPSK transmitter including frame format, modulation scheme and
% etc. With prmQPSKTransmitter, the baseband signal processing is set.
% prmPAControl: a structure array that groups the related parameters
% on phased arrays including available codebooks, receiving power gain,
% antenna switching on/off settings and etc. With prmPAControl, the RF 
% configuration is set.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function run_SDRuQPSKTransmitter(prmQPSKTransmitter, prmPAControl)
%#codegen
global Raw_originalSignal

%% Create and configure the SDRu system object for N210
persistent radio
if isempty(radio)   
     radio = comm.SDRuTransmitter(...
         'Platform',                prmQPSKTransmitter.Platform, ...
         'IPAddress',               prmQPSKTransmitter.Address, ...
         'CenterFrequency',         prmQPSKTransmitter.USRPCenterFrequency, ...
         'Gain',                    prmQPSKTransmitter.USRPGain, ...
         'InterpolationFactor',     prmQPSKTransmitter.USRPInterpolationFactor, ...
         'LocalOscillatorOffset',   prmQPSKTransmitter.USRPLocalOscillatorOffset, ...
         'TransportDataType',       prmQPSKTransmitter.USRPTransportDataType, ...
         'EnableBurstMode',         prmQPSKTransmitter.USRPEnableBurstMode, ...
         'NumFramesInBurst',        prmQPSKTransmitter.USRPNumFramesInBurst, ...
         'ClockSource',             prmQPSKTransmitter.USRPClockSource, ...
         'PPSSource',               prmQPSKTransmitter.USRPPPSSource); 
end    


%% Get information about the used USRP if the program is not run with MATLAB Codegen
if coder.target('MATLAB')
    radioinfo = info(radio);
    disp(radioinfo);
end


%% Open serial port and start transmission process
% The setting of holding time is to make sure the configuration of the
% phased array is successfully complete. The following deterministic values
% is going to be used by phased array interface functions detailed in
% PART.c
Holding_Time_Open_Serial = prmPAControl.Holding_Time_Open_Serial;
Holding_Time_For_Writing_Buffer = prmPAControl.Holding_Time_For_Writing_Buffer;
Holding_Time_Turn_On_Antenna = prmPAControl.Holding_Time_Turn_On_Antenna;
Holding_Time_Change_Phase_Shifter_State = prmPAControl.Holding_Time_Change_Phase_Shifter_State;


%% The following parts works only with MATLAB codegen
if ~coder.target('MATLAB')
    switch prmPAControl.Program_ID
        %% To measure the beam pattern (unfinished)
        case 7
            disp('-----------------------------------------------------');
            disp('---------- Start Show Beam Pattern Program ----------');
            disp('-----------------------------------------------------');
            %Open Tx serial
            coder.ceval('OpenSerialPortTx',prmPAControl.SerialPortTx);
            coder.ceval('Operation_Interval_Time',Holding_Time_Open_Serial);
            %Turn on 1-st transmitter antenna
            coder.ceval('WriteCommandTx',prmPAControl.txon(1,:),prmPAControl.SerialPortTx);
            coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
            coder.ceval('ReadCommandTxPrint');
            %Close serial
            coder.ceval('CloseSerialPortTx');
            currentTime = 0;
            % Display the start of the transmission
            disp('********** Start Transmitting Signal **********');        
            while currentTime < prmPAControl.StopTransmissionTime
                step(radio, Raw_originalSignal(:,1)); % Always send the first frame
                currentTime = currentTime+prmQPSKTransmitter.FrameTime;
            end    
            disp('********** Stop Transmitting Signal **********');
            disp('----------------------------------------------------');
            disp('---------- Stop Show Beam Pattern Program ----------');
            disp('----------------------------------------------------');
            
        %% Power calibration
        case 0
            disp('-----------------------------------------------------');
            disp('---------- Start Power Calibration Program ----------');
            disp('-----------------------------------------------------');
            currentTime = 0;
            % Display the start of the transmission
            disp('********** Start Transmitting Signal **********');        
            while currentTime < prmPAControl.StopTransmissionTime
                step(radio, Raw_originalSignal(:,1)); %Always send the first frame
                currentTime = currentTime+prmQPSKTransmitter.FrameTime;
            end    
            disp('********** Stop Transmitting Signal **********');
            disp('----------------------------------------------------');
            disp('---------- Stop Power Calibration Program ----------');
            disp('----------------------------------------------------');
            
            
        %% Transmitter phased array calibration
        case 2
            disp('-------------------------------------------------------------------------');
            disp('---------- Start Transmitter Phase Shifter Calibration Program ----------');
            disp('-------------------------------------------------------------------------');
            %Open serial
            coder.ceval('OpenSerialPortTx',prmPAControl.SerialPortTx);
            coder.ceval('Operation_Interval_Time',Holding_Time_Open_Serial);
            %Turn on 1-st antenna
            coder.ceval('WriteCommandTx',prmPAControl.txon(1,:),prmPAControl.SerialPortTx);
            coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer*2);
            coder.ceval('ReadCommandTxPrint');
            %Reset phase shifters' state
            coder.ceval('WriteCommandTx','phvec 0',prmPAControl.SerialPortTx);
            coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer*2);
            coder.ceval('ReadCommandTxPrint');
            % Display the start of the transmission
            disp('********** Start Transmitting Signal **********');
            repeat = prmPAControl.Total_Frame_To_Transmit/prmPAControl.Number_of_Frames_Per_Test_Tx/...
                    (prmPAControl.Number_antenna_to_calibrate*prmPAControl.Number_Phase_Shifter_State);
            for r=1:1:repeat
                for state=1:1:prmPAControl.Number_antenna_to_calibrate*prmPAControl.Number_Phase_Shifter_State                 
                    for nf = 1:1:prmPAControl.Number_of_Frames_Per_Test_Tx
                        step(radio, Raw_originalSignal(:,state));
                    end
                    coder.ceval('WriteCommandTx',prmPAControl.txon(ceil(state/4)+1,:),prmPAControl.SerialPortTx);
                    coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
                    coder.ceval('WriteCommandTx',prmPAControl.phvec(state,:),prmPAControl.SerialPortTx);
                    %--------------------------------------------------------------------------------------------------
                    %---------------------If other settings are changed, time should be remeasured---------------------
                    %--------------------------------------------------------------------------------------------------
                    coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer*2); 
                    %--------------------------------------------------------------------------------------------------
                end
            end
            coder.ceval('CloseSerialPortTx');
            disp('********** Stop Transmitting Signal **********');
            disp('------------------------------------------------------------------------');
            disp('---------- Stop Transmitter Phase Shifter Calibration Program ----------');
            disp('------------------------------------------------------------------------');
            
            
        %% Transmitter or overall beam training
        case {4,6}
            if prmPAControl.Program_ID == 4
                disp('-------------------------------------------------------------');
                disp('---------- Start Transmitter Beam Training Program ----------');    
                disp('-------------------------------------------------------------');
            else
                disp('-------------------------------------------------------------');
                disp('---------- Start Overall Beam Training Program --------------');    
                disp('-------------------------------------------------------------');
            end
            %Open serial
            coder.ceval('OpenSerialPortTx',prmPAControl.SerialPortTx);
            coder.ceval('Operation_Interval_Time',Holding_Time_Open_Serial);
            %Turn on all antennas
            coder.ceval('WriteCommandTx',prmPAControl.txonall,prmPAControl.SerialPortTx);
            coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
            coder.ceval('ReadCommandTxPrint'); 
            %Reset phase shifters' state
            coder.ceval('WriteCommandTx','phvec 0',prmPAControl.SerialPortTx);
            coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
            coder.ceval('ReadCommandTxPrint');
            % Display the start of the transmission
            disp('********** Start Transmitting Signal **********');   
            for r=1:1:(prmPAControl.Total_Frame_To_Transmit/prmPAControl.Number_of_Frames_Per_Test_Tx/11)
                for state=1:1:prmPAControl.Number_CodeBook_Tx 
                    for nf = 1:1:prmPAControl.Number_of_Frames_Per_Test_Tx
                        step(radio, Raw_originalSignal(:,state));
                    end
                    coder.ceval('WriteCommandTx',prmPAControl.CodeBook_Tx(state,:),prmPAControl.SerialPortTx);
                    %--------------------------------------------------------------------------------------------------
                    %---------------------If other settings are changed, time should be remeasured---------------------
                    %--------------------------------------------------------------------------------------------------
                    coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer*5.5); 
                    %--------------------------------------------------------------------------------------------------
                end
            end
            disp('********** Stop Transmitting Signal **********');
            coder.ceval('CloseSerialPortTx');
            if prmPAControl.Program_ID == 4
                disp('-------------------------------------------------------------');
                disp('---------- Stop Transmitter Beam Training Program -----------');    
                disp('-------------------------------------------------------------');
            else %prmPAControl.Program_ID == 6
                disp('-------------------------------------------------------------');
                disp('---------- Stop Overall Beam Training Program ---------------');    
                disp('-------------------------------------------------------------');
            end
            
            
        %% Others
        case {1, 3, 5, 8, 2.1, 4.1, 8.1}
            if prmPAControl.Program_ID == 1
                disp('-----------------------------------------------------');
                disp('---------- Start Power Level Test Program -----------');   
                disp('-----------------------------------------------------');
            elseif prmPAControl.Program_ID == 3
                disp('-----------------------------------------------------------------------');
                disp('---------- Start Receiver Phase Shifter Calibration Program -----------');
                disp('-----------------------------------------------------------------------');
            elseif prmPAControl.Program_ID == 5
                disp('-----------------------------------------------------------');
                disp('---------- Start Receiver Beam Training Program -----------');
                disp('-----------------------------------------------------------');
            elseif prmPAControl.Program_ID == 8
                disp('-----------------------------------------------------------');
                disp('---------- Start Chanel Estimation Program ---------------');
                disp('-----------------------------------------------------------');
            elseif prmPAControl.Program_ID == 2.1
                disp('---------------------------------------------------------------------------------');
                disp('---------- Start Transmitter Phase Shifter Calibration Program (Wired)-----------');
                disp('---------------------------------------------------------------------------------');
            elseif prmPAControl.Program_ID == 4.1
                disp('----------------------------------------------------------------------');
                disp('---------- Start Transmitter Beam Training Program (Wired) -----------');
                disp('----------------------------------------------------------------------');
            else %prmPAControl.Program_ID == 8.1
                disp('----------------------------------------------------------------------');
                disp('---------- Start Overall Beam Training Program (Wired) ---------------');
                disp('----------------------------------------------------------------------');
            end
                            
            if prmPAControl.Program_ID == 1 || prmPAControl.Program_ID == 3 || prmPAControl.Program_ID == 5 ...
                    || prmPAControl.Program_ID == 8
                %Open Tx serial
                coder.ceval('OpenSerialPortTx',prmPAControl.SerialPortTx);
                coder.ceval('Operation_Interval_Time',Holding_Time_Open_Serial);
                %Turn on 1-st transmitter antenna
                ind_on = 1;  
                coder.ceval('WriteCommandTx',prmPAControl.txon_single(ind_on,:),prmPAControl.SerialPortTx);
                %coder.ceval('WriteCommandTx','txon 0x001',prmPAControl.SerialPortTx);
                coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
                coder.ceval('ReadCommandTxPrint');
                %Reset phase vector to "0x0000000"
                coder.ceval('WriteCommandTx', 'phvec 0', prmPAControl.SerialPortTx);
                coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
                %phase_on = 3;
                %coder.ceval('WriteCommandTx',prmPAControl.phvec(4*(ind_on-2)+phase_on,:),prmPAControl.SerialPortTx);
                coder.ceval('Operation_Interval_Time',Holding_Time_For_Writing_Buffer);
                coder.ceval('ReadCommandTxPrint');
                %Close serial
                coder.ceval('CloseSerialPortTx');
            elseif prmPAControl.Program_ID == 8.1
                ;
            else          
                %Open Rx serial
                coder.ceval('OpenSerialPortRx',prmPAControl.SerialPortRx);
                coder.ceval('Operation_Interval_Time',Holding_Time_Open_Serial);
                %Turn on 1-st antenna
                coder.ceval('WriteCommandRx',prmPAControl.rxon(1,:),prmPAControl.SerialPortRx);
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
                %Close Rx serial
                coder.ceval('CloseSerialPortRx');
            end
            
            % Start of the transmission
            currentTime = 0.00;
            Total_Frame_Transmitted = 0;
            NthFrame = 0;
            Count_Transmitted_Frame_Normal = 1;
            disp('********** Start Transmitting Signal **********')
            while currentTime < prmPAControl.StopTransmissionTime
                if Count_Transmitted_Frame_Normal > prmPAControl.Number_of_Frames_Per_Test_Tx  
                    NthFrame = NthFrame + 1;
                    Count_Transmitted_Frame_Normal = 1;
                end
                Index_of_raw_signal = mod(NthFrame,prmPAControl.Number_State_To_Test)+1;
                underrun = step(radio, Raw_originalSignal(:, Index_of_raw_signal));
                if underrun && ~prmQPSKTransmitter.USRPEnableBurstMode
                    disp('Underrun detected');
                end
                % Update simulation time
                currentTime = currentTime+prmQPSKTransmitter.FrameTime;
                Count_Transmitted_Frame_Normal = Count_Transmitted_Frame_Normal + 1;
                Total_Frame_Transmitted = Total_Frame_Transmitted + 1;
            end
            disp('********** Stop Transmitting Signal **********');
            if prmPAControl.Program_ID == 1
                disp('----------------------------------------------------');
                disp('---------- Stop Power Level Test Program -----------');   
                disp('----------------------------------------------------');
            elseif prmPAControl.Program_ID == 3
                disp('----------------------------------------------------------------------');
                disp('---------- Stop Receiver Phase Shifter Calibration Program -----------');
                disp('----------------------------------------------------------------------');
            elseif prmPAControl.Program_ID == 5
                disp('----------------------------------------------------------');
                disp('---------- Stop Receiver Beam Training Program -----------');
                disp('----------------------------------------------------------');
            elseif prmPAControl.Program_ID == 8
                disp('----------------------------------------------------------');
                disp('---------- Stop Channel Estimation Program -----------');
                disp('----------------------------------------------------------');
            elseif prmPAControl.Program_ID == 2.1
                disp('--------------------------------------------------------------------------------');
                disp('---------- Stop Transmitter Phase Shifter Calibration Program (Wired)-----------');
                disp('--------------------------------------------------------------------------------');
            elseif  prmPAControl.Program_ID == 4.1
                disp('----------------------------------------------------------------------');
                disp('---------- Stop Transmitter Beam Training Program (Wired) -----------');
                disp('----------------------------------------------------------------------'); 
            else  % case prmPAControl.Program_ID == 8.1
                disp('----------------------------------------------------------------------');
                disp('---------- Stop Overrall Beam Training Program (Wired) -----------');
                disp('----------------------------------------------------------------------');
            end  
    end        
end
%% Release system objects
release(radio);


