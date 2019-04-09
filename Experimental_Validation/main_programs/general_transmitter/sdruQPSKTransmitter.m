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
% This main script runs the transmitter which includes the process of data
% generation, Tx USRP configuration, antenna configuration and data 
% transmission.
%--------Program_ID = 0: Power calibration for the google phased array.
%--------Program_ID = 1: Power level test. This program demonstrates a 
%                        representative part of the available receiving 
%                        power setting. 
%--------Program_ID = 2: Transmitter phased array calibration. This program 
%                        requires the receiver to recognize that the 
%                        received signal is transmitted via which beam 
%                        pattern by decoding the information of the frame 
%                        transmitted.
%--------Program_ID = 2.1: Transmitter phased array calibration with phased 
%                          array control at Rx side. This program does not 
%                          require the decoding part at the Rx side as the 
%                          Tx phased array configuration is controlled at
%                          the Rx host computer for convenience.
%--------Program_ID = 3: Receiver phased array calibration. This program 
%                        does not require the decoding part at the Rx side 
%                        as the Tx phased array is fixed during the
%                        receiver phased array calibration.
%--------Program_ID = 4: Transmitter beam training. This program 
%                        requires the receiver to recognize that the 
%                        received signal is transmitted via which beam 
%                        pattern by decoding the information of the frame 
%                        transmitted.
%--------Program_ID = 4.1: Transmitter Beam Training with phased 
%                          array control at Rx side. This program does not 
%                          require the decoding part at the Rx side as the 
%                          Tx phased array configuration is controlled at
%                          the Rx host computer for convenience.
%--------Program_ID = 5: Receiver beam training. This program 
%                        does not require the decoding part at the Rx side 
%                        as the Tx phased array is fixed during the
%                        receiver phased array calibration.
%--------Program_ID = 6: Overall beam training. This program 
%                        requires the receiver to recognize that the 
%                        received signal is transmitted via which beam 
%                        pattern by decoding the information of the frame 
%                        transmitted. It performs double-side beam
%                        training.
%--------Program_ID = 7: Show beam pattern. This program allows measuring
%                        beam pattern.
%--------Program_ID = 8: Channel estimation. This program allows performing
%                        channel estimation in a non-coherent way.
%--------Program_ID = 8.1: Overall beam training with phased array control
%                          at Rx side. This program does not require the 
%                          decoding part at the Rx side as both Tx and Rx 
%                          phased arrays are controlled at the Rx side.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Add paths
folder_name = ...
[
"/Experimental_Validation/";
];
for i=1:length(folder_name)
    Lib_path = char(strcat(pwd,folder_name(i)));
    addpath(genpath(Lib_path));
end


%% Connect to the USRP and phased array
clear all
address = '192.168.10.2';
%address = '192.168.20.2';
platform = 'N200/N210/USRP2';
SerialPortTx = ['/dev/tty.usbmodem3441']; % USB serial port for Tx phased array, which need to be checked before experiments.
SerialPortRx = ['/dev/tty.usbmodem221'];  % USB serial port for Rx phased array, which need to be checked before experiments.
Program_ID = 8.1;

%% Initialization
% Transmitter parameter structure
prmQPSKTransmitter = sdruqpsktransmitter_init(platform);
prmQPSKTransmitter.Platform = platform;
prmQPSKTransmitter.Address = address;

% Phase array control parameter structure
[prmPAControl, type_prmPAControl] = phase_array_control_init(prmQPSKTransmitter, SerialPortTx, SerialPortRx, Program_ID);

%For default setting, the burst mode is false
prmQPSKTransmitter.USRPEnableBurstMode = false;
prmQPSKTransmitter.USRPNumFramesInBurst = prmPAControl.Number_of_Frames_Per_Test_Tx;

% Code generation usage
reGeneration = true;
compileIt = true; % true if code is to be compiled for accelerated execution 
useCodegen = true; % true to run the latest generated mex file


%% Raw transmitted data storage and generate offline data
global Raw_originalSignal
global Raw_originalSignal_Tx_CE
global Raw_originalSignal_Normal
global Raw_originalSignal_Tx_Calibration
global Raw_originalSignal_Tx_Beamtraining
global Baseband_y_original
SamplesPerFrame_Initial = prmQPSKTransmitter.FrameSize*prmQPSKTransmitter.Upsampling;
Raw_originalSignal_Tx_CE = complex(zeros(SamplesPerFrame_Initial,int32(44)));
Raw_originalSignal_Normal = complex(zeros(SamplesPerFrame_Initial,int32(44)));
Raw_originalSignal_Tx_Calibration = complex(zeros(SamplesPerFrame_Initial,int32(44)));
Raw_originalSignal_Tx_Beamtraining = complex(zeros(SamplesPerFrame_Initial,int32(11)));
Raw_originalSignal = complex(zeros(SamplesPerFrame_Initial,int32(44)));
Baseband_y_original = complex(zeros(prmQPSKTransmitter.FrameSize,int32(44)));


%% Compliation
%Complie and regenerate transmit data set
if reGeneration
    codegen('run_Offline_Data_Generator', '-args', {coder.Constant(prmQPSKTransmitter), type_prmPAControl}); %#ok<UNRCH>
    clear run_Offline_Data_Generator_mex %#ok<UNRCH>
    offline_transmit_data_path = char(strcat(pwd,'/Experimental_Validation/src/offline_transmit_data'));
    save(char(strcat(offline_transmit_data_path,'/Raw_originalSignal_Normal.mat')),'Raw_originalSignal_Normal');
    save(char(strcat(offline_transmit_data_path,'/Raw_originalSignal_Tx_Calibration.mat')),'Raw_originalSignal_Tx_Calibration');
    save(char(strcat(offline_transmit_data_path,'/Raw_originalSignal_Tx_Beamtraining.mat')),'Raw_originalSignal_Tx_Beamtraining');
    save(char(strcat(offline_transmit_data_path,'/Raw_originalSignal_Tx_CE.mat')),'Raw_originalSignal_Tx_CE');
    save(char(strcat(offline_transmit_data_path,'/Baseband_y_original.mat')),'Baseband_y_original');
end
switch prmPAControl.Program_ID
    case {0,1,3,5,2.1,4.1,7}
        load('Raw_originalSignal_Normal.mat');
        Raw_originalSignal(:,1:prmPAControl.Number_kinds_of_Offline_Frames_Normal) = Raw_originalSignal_Normal;
    case 2
        load('Raw_originalSignal_Tx_Calibration.mat');
        Raw_originalSignal(:,1:prmPAControl.Number_kinds_of_Offline_Frames_Calibration) = Raw_originalSignal_Tx_Calibration;
    case {4, 6}
        load('Raw_originalSignal_Tx_Beamtraining.mat');
        Raw_originalSignal(:,1:prmPAControl.Number_kinds_of_Offline_Frames_Beamtraining) = Raw_originalSignal_Tx_Beamtraining;
    case {8, 8.1}
        load('Raw_originalSignal_Tx_CE.mat');
        Raw_originalSignal(:,1:prmPAControl.Number_kinds_of_Offline_Frames_Normal) = Raw_originalSignal_Tx_CE;
end
%Complie the transmit
if compileIt
    codegen('run_SDRuQPSKTransmitter','PATR.c', 'PATR.h', '-args', {coder.Constant(prmQPSKTransmitter), type_prmPAControl},...
                                                          '-I', char(strcat(pwd,'/Experimental_Validation/src/array_control')));%#ok<UNRCH>
end


%% Execution of the transmitter
if useCodegen
    clear run_SDRuQPSKTransmitter_mex %#ok<UNRCH>
    clear run_SDRuQPSKReceiver_mex %#ok<UNRCH>
    clear Control_Phase_Array_External_mex %#ok<UNRCH>
    clear Control_Phase_Array_External_Phase_Error_Calibration_mex %#ok<UNRCH>
    run_SDRuQPSKTransmitter_mex(prmQPSKTransmitter, prmPAControl);  
else
    run_Offline_Data_Generator(prmQPSKTransmitter, prmPAControl);
    offline_transmit_data_path = char(strcat(pwd,'/Transmitter/offline_transmit_data'));
    save(char(strcat(offline_transmit_data_path,'/Raw_originalSignal_Normal.mat')),'Raw_originalSignal_Normal');
    save(char(strcat(offline_transmit_data_path,'/Raw_originalSignal_Tx_Calibration.mat')),'Raw_originalSignal_Tx_Calibration');
    save(char(strcat(offline_transmit_data_path,'/Raw_originalSignal_Tx_Beamtraining.mat')),'Raw_originalSignal_Tx_Beamtraining');
    save(char(strcat(offline_transmit_data_path,'/Raw_originalSignal_Tx_CE.mat')),'Raw_originalSignal_Tx_CE');
    save(char(strcat(offline_transmit_data_path,'/Baseband_y_original.mat')),'Baseband_y_original');
    switch prmPAControl.Program_ID
        case {0,1,3,5,2.1,4.1,7}
            load('Raw_originalSignal_Normal.mat');
            Raw_originalSignal(:,1:prmPAControl.Number_kinds_of_Offline_Frames_Normal) = Raw_originalSignal_Normal;
        case 2
            load('Raw_originalSignal_Tx_Calibration.mat');
            Raw_originalSignal(:,1:prmPAControl.Number_kinds_of_Offline_Frames_Calibration) = Raw_originalSignal_Tx_Calibration;
        case {4, 6}
            load('Raw_originalSignal_Tx_Beamtraining.mat');
            Raw_originalSignal(:,1:prmPAControl.Number_kinds_of_Offline_Frames_Beamtraining) = Raw_originalSignal_Tx_Beamtraining;
        case {8, 8.1}
            load('Raw_originalSignal_Tx_CE.mat');
            Raw_originalSignal(:,1:prmPAControl.Number_kinds_of_Offline_Frames_Normal) = Raw_originalSignal_Tx_CE;
    end
    run_SDRuQPSKTransmitter(prmQPSKTransmitter, prmPAControl);
end


