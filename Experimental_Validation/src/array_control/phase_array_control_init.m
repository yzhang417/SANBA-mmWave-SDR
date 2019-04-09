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
% This function Outputs a structure array prmPAControl that groups the 
% related parameters on the phased array configuration parameters. This
% structure array is widely used in the testbed as it serves as an
% instruction for most programs involving the phased array configuration.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% prmQPSK: it could be a structure array that either groups the related 
% parameters on the QPSK transmitter or receiver, i.e., variable 
% prmQPSKReceiver in [Experimental_Validation\main_programs\
% general_receiver_decoder\sdruQPSKReceiver.m] or prmQPSKtransmitter in 
% [Experimental_Validation\main_programs\general_transmitter\
% sdruQPSKTransmitter.m]
% SerialPortTx: USB serial port for Tx phased array, which needs to be 
% checked before experiments.
% SerialPortRx: USB serial port for Rx phased array, which needs to be 
% checked before experiments.
% Program_ID:
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
%--------Program_ID = 9: Antenna element gain calibration. This program
%                        collects the data for calculating the antenna gain
%                        for each single element.
%--------Program_ID = 9.1: Fine phase error calibration for TX. This
%                          program collects the data for our phase match 
%                          based calibration (PMBC). 
%--------Program_ID = 9.2: Fine phase error calibration for RX. This
%                          program collects the data for our phase match 
%                          based calibration (PMBC). 
% varargin: used to revise the number of codebook to train, whose default 
% value is 11. For the experiment, the number of codebook to train is 6 and
% 8. % Format for varargin: 'property_name_1", value_1, 
% "property_name_2", value_2, ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output arguments:
% prmPAControl: a structure array prmPAControl that groups the 
% related parameters on the phased array configuration parameters. This
% structure array is widely used in this testbed as it serves as an
% instructions for most programs involving the phased array configurations.
% type_prmPAControl: the type of prmPAControl, which will be used when
% using Codegen to compile some functions for speeding up the program.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [prmPAControl, type_prmPAControl] = phase_array_control_init(prmQPSK, SerialPortTx, SerialPortRx, Program_ID, varargin)
%% Serial port to Tx ans Rx phased array
%
prmPAControl.SerialPortTx = SerialPortTx;
type_SerialPortTx = coder.newtype('char',[1 inf]);
%
prmPAControl.SerialPortRx = SerialPortRx;
type_SerialPortRx = coder.newtype('char',[1 inf]);                  

%% Holding time for phased array operations
prmPAControl.Holding_Time_Open_Serial = 1000000;
[r, c] = size(prmPAControl.Holding_Time_Open_Serial);
type_Holding_Time_Open_Serial = coder.newtype('double',[r c]);

prmPAControl.Holding_Time_For_Writing_Buffer = 3000;
[r, c] = size(prmPAControl.Holding_Time_For_Writing_Buffer);
type_Holding_Time_For_Writing_Buffer = coder.newtype('double',[r c]);

prmPAControl.Holding_Time_Turn_On_Antenna = 2000;
[r, c] = size(prmPAControl.Holding_Time_Turn_On_Antenna);
type_Holding_Time_Turn_On_Antenna = coder.newtype('double',[r c]);

prmPAControl.Holding_Time_Change_Phase_Shifter_State = 1500;
[r, c] = size(prmPAControl.Holding_Time_Change_Phase_Shifter_State);
type_Holding_Time_Change_Phase_Shifter_State = coder.newtype('double',[r c]);

prmPAControl.Holding_Time_Self_Defined = 1000;
[r, c] = size(prmPAControl.Holding_Time_Self_Defined);
type_Holding_Time_Self_Defined = coder.newtype('double',[r c]);

%% Antenna and phase shifter number
%
prmPAControl.Number_Phase_Shifter_State = 4;
[r, c] = size(prmPAControl.Number_Phase_Shifter_State);
type_Number_Phase_Shifter_State = coder.newtype('double',[r c]);

%
prmPAControl.Number_antenna = 12;
[r, c] = size(prmPAControl.Number_antenna);
type_Number_antenna = coder.newtype('double',[r c]);

%% Turning on/off antenna command (UART commands)
%
prmPAControl.txoffall = ['txon 0'];
[r, c] = size(prmPAControl.txoffall);
type_txoffall = coder.newtype('char',[r c]);
%
prmPAControl.txonall = ['tx12'];
[r, c] = size(prmPAControl.txonall);
type_txonall = coder.newtype('char',[r c]);
%
prmPAControl.rxoffall = ['rxon 0'];
[r, c] = size(prmPAControl.rxoffall);
type_rxoffall = coder.newtype('char',[r c]);
%
prmPAControl.rxonall = ['rx12'];
[r, c] = size(prmPAControl.rxonall);
type_rxonall = coder.newtype('char',[r c]);

%% Turning on/off single antennas (UART commands)
%
prmPAControl.txon_single = ['txon 0x001';'txon 0x002';'txon 0x004';'txon 0x008';...
                            'txon 0x010';'txon 0x020';'txon 0x040';'txon 0x080';...
                            'txon 0x100';'txon 0x200';'txon 0x400';'txon 0x800'];  
[r, c] = size(prmPAControl.txon_single);
type_txon_single = coder.newtype('char',[r c]);

%                
prmPAControl.rxon_single = ['rxon 0x001';'rxon 0x002';'rxon 0x004';'rxon 0x008';...
                            'rxon 0x010';'rxon 0x020';'rxon 0x040';'rxon 0x080';...
                            'rxon 0x100';'rxon 0x200';'rxon 0x400';'rxon 0x800'];
[r, c] = size(prmPAControl.rxon_single);
type_rxon_single = coder.newtype('char',[r c]);

%% Turning on only two antennas for calibrations (UART commands)
%
prmPAControl.txon = ['txon 0x001';'txon 0x003';'txon 0x005';'txon 0x009';...
                     'txon 0x011';'txon 0x021';'txon 0x041';'txon 0x081';...
                     'txon 0x101';'txon 0x201';'txon 0x401';'txon 0x801'];  
[r, c] = size(prmPAControl.txon);
type_txon = coder.newtype('char',[r c]);

%                
prmPAControl.rxon = ['rxon 0x001';'rxon 0x003';'rxon 0x005';'rxon 0x009';...
                     'rxon 0x011';'rxon 0x021';'rxon 0x041';'rxon 0x081';...
                     'rxon 0x101';'rxon 0x201';'rxon 0x401';'rxon 0x801'];
[r, c] = size(prmPAControl.rxon);
type_rxon = coder.newtype('char',[r c]);

%
prmPAControl.phvec_Antenna1 = ['phvec 0x000000';'phvec 0x000001';'phvec 0x001000';'phvec 0x001001'];
[r, c] = size(prmPAControl.phvec_Antenna1);
type_phvec_Antenna1 = coder.newtype('char',[r c]);

prmPAControl.phvec = ['phvec 0x000000';'phvec 0x000002';'phvec 0x002000';'phvec 0x002002';...
                      'phvec 0x000000';'phvec 0x000004';'phvec 0x004000';'phvec 0x004004';...
                      'phvec 0x000000';'phvec 0x000008';'phvec 0x008000';'phvec 0x008008';...
                      'phvec 0x000000';'phvec 0x000010';'phvec 0x010000';'phvec 0x010010';...
                      'phvec 0x000000';'phvec 0x000020';'phvec 0x020000';'phvec 0x020020';...
                      'phvec 0x000000';'phvec 0x000040';'phvec 0x040000';'phvec 0x040040';...
                      'phvec 0x000000';'phvec 0x000080';'phvec 0x080000';'phvec 0x080080';...
                      'phvec 0x000000';'phvec 0x000100';'phvec 0x100000';'phvec 0x100100';...
                      'phvec 0x000000';'phvec 0x000200';'phvec 0x200000';'phvec 0x200200';...
                      'phvec 0x000000';'phvec 0x000400';'phvec 0x400000';'phvec 0x400400';...
                      'phvec 0x000000';'phvec 0x000800';'phvec 0x800000';'phvec 0x800800'];                                    
[r, c] = size(prmPAControl.phvec);
type_phvec = coder.newtype('char',[r c]);

%
prmPAControl.Number_phvec = length(prmPAControl.phvec(:,1));
[r, c] = size(prmPAControl.Number_phvec);
type_Number_phvec = coder.newtype('double',[r c]);

%
List_antenna_to_calibrate = [2:12]; %Start with 2 to 12
prmPAControl.List_Index_antenna_to_calibrate = zeros(1,12);
prmPAControl.List_Index_antenna_to_calibrate(List_antenna_to_calibrate) = 1;
[r, c] = size(prmPAControl.List_Index_antenna_to_calibrate);
type_List_Index_antenna_to_calibrate = coder.newtype('double',[r c]);

%
prmPAControl.Number_antenna_to_calibrate = length(List_antenna_to_calibrate);
[r, c] = size(prmPAControl.Number_antenna_to_calibrate);
type_Number_antenna_to_calibrate = coder.newtype('double',[r c]);

%% Revise the number of codebook to train, whose default value is 11
num_codebook_to_train = 11;
n_var_in = length(varargin);
if mod(n_var_in,2) ~= 0
    error('Number of input arguments should be even.')
    % Format for varargin: 'property_name_1", value_1, "property_name_2", value_2, ...
else
    for k = 1 : 2 : n_var_in
        field_name = varargin{k};
        new_value = varargin{k+1};
        if strcmp(field_name,'num_codebook_to_train')
            num_codebook_to_train = new_value;
        end
    end
end
if num_codebook_to_train == 6
    list_active_beam_label = linspace(-10,10,num_codebook_to_train);
else
    list_active_beam_label = [linspace(-10,10,num_codebook_to_train-1) 0];
end
prmPAControl.beam_label = zeros(1,101);
prmPAControl.beam_label(round(list_active_beam_label+50)) = 1;
[r, c] = size(prmPAControl.beam_label);
type_beam_label = coder.newtype('double',[r c]);

%% Transmitter codebook   
% Phvec_command_no_calibration =['phvec 0x70e936'];
% Phvec_command_coarse =        ['phvec 0x464480'];
% Phvec_command_fine_m1 =       ['phvec 0x460014'];
% Phvec_command_fine_m2 =       ['phvec 0x1e0016'];
%
Phvec_command =['phvec xxxxxxxx'];
CodeBook_Size = 12; % not greater than 44
prmPAControl.CodeBook_Tx = repmat(Phvec_command,CodeBook_Size,1);
[r, c] = size(prmPAControl.CodeBook_Tx);
type_CodeBook_Tx = coder.newtype('char',[r c]);

%
prmPAControl.Number_CodeBook_Tx = length(prmPAControl.CodeBook_Tx(:,1));
[r, c] = size(prmPAControl.Number_CodeBook_Tx);
type_Number_CodeBook_Tx = coder.newtype('double',[r c]);

%
List_CodeBook_Tx_to_train = 1:num_codebook_to_train; %Between 1 to 12
prmPAControl.List_Index_CodeBook_Tx_to_train = zeros(1,CodeBook_Size);
prmPAControl.List_Index_CodeBook_Tx_to_train(List_CodeBook_Tx_to_train) = 1;
[r, c] = size(prmPAControl.List_Index_CodeBook_Tx_to_train);
type_List_Index_CodeBook_Tx_to_train = coder.newtype('double',[r c]);

%
prmPAControl.Number_CodeBook_Tx_to_train = length(List_CodeBook_Tx_to_train);
[r, c] = size(prmPAControl.Number_CodeBook_Tx_to_train);
type_Number_CodeBook_Tx_to_train = coder.newtype('double',[r c]);

%% Receiver Codebook
% Phvec_command_no_calibration =['phvec 0x000000'];
% Phvec_command_coarse =        ['phvec 0xa5cdb6'];
% Phvec_command_fine_m1 =       ['phvec 0xa5cdb6'];
% Phvec_command_fine_m2 =       ['phvec 0xffe180'];
%
Phvec_command =['phvec xxxxxxxx'];
prmPAControl.CodeBook_Rx = repmat(Phvec_command,CodeBook_Size,1);   
[r, c] = size(prmPAControl.CodeBook_Rx);
type_CodeBook_Rx = coder.newtype('char',[r c]);

%
prmPAControl.Number_CodeBook_Rx = length(prmPAControl.CodeBook_Rx(:,1));
[r, c] = size(prmPAControl.Number_CodeBook_Rx);
type_Number_CodeBook_Rx = coder.newtype('double',[r c]);

%
List_CodeBook_Rx_to_train = 1:num_codebook_to_train; %Between 1 to 12
prmPAControl.List_Index_CodeBook_Rx_to_train = zeros(1,CodeBook_Size);
prmPAControl.List_Index_CodeBook_Rx_to_train(List_CodeBook_Rx_to_train) = 1;
[r, c] = size(prmPAControl.List_Index_CodeBook_Rx_to_train);
type_List_Index_CodeBook_Rx_to_train = coder.newtype('double',[r c]);

%
prmPAControl.Number_CodeBook_Rx_to_train = length(List_CodeBook_Rx_to_train);
[r, c] = size(prmPAControl.Number_CodeBook_Rx_to_train);
type_Number_CodeBook_Rx_to_train = coder.newtype('double',[r c]);


%% Receiving power configuration (UART commands)
%
prmPAControl.rfvga = ['rfvga 0';'rfvga 1';'rfvga 2';'rfvga 3';'rfvga 4'];
[r, c] = size(prmPAControl.rfvga);
type_rfvga = coder.newtype('char',[r c]);

%
prmPAControl.rxvga = ['rxvga  0';'rxvga  1';'rxvga  2';'rxvga  3';'rxvga  4';...
                      'rxvga  5';'rxvga  6';'rxvga  7';'rxvga  8';'rxvga  9';...
                      'rxvga 10';'rxvga 11';'rxvga 12';'rxvga 13';'rxvga 14';...
                      'rxvga 15';'rxvga 16';'rxvga 17';'rxvga 18';'rxvga 19';...
                      'rxvga 20';'rxvga 21';'rxvga 22';'rxvga 23';'rxvga 24';...
                      'rxvga 25';'rxvga 26';'rxvga 27';'rxvga 28';'rxvga 29';... 
                      'rxvga 30';'rxvga 31';'rxvga 32';'rxvga 33';'rxvga 34';...
                      'rxvga 35';'rxvga 36';'rxvga 37';'rxvga 38';'rxvga 39';...
                      'rxvga 40';'rxvga 41';'rxvga 42';'rxvga 43';'rxvga 44';...
                      'rxvga 45';'rxvga 46';'rxvga 47'];   
[r, c] = size(prmPAControl.rxvga);
type_rxvga = coder.newtype('char',[r c]);

%
List_rfvga_to_test = [1:4]+1;
prmPAControl.List_Index_rfvga_to_test = zeros(1,5);
prmPAControl.List_Index_rfvga_to_test(List_rfvga_to_test) = 1;
[r, c] = size(prmPAControl.List_Index_rfvga_to_test);
type_List_Index_rfvga_to_test = coder.newtype('double',[r c]);

%
prmPAControl.Number_rfvga_to_test = length(List_rfvga_to_test);
[r, c] = size(prmPAControl.Number_rfvga_to_test);
type_Number_rfvga_to_test = coder.newtype('double',[r c]);

%
List_rxvga_to_test = [0 5 10 15 20 25 30 35 40 45 47]+1;
prmPAControl.List_Index_rxvga_to_test = zeros(1,48);
prmPAControl.List_Index_rxvga_to_test(List_rxvga_to_test) = 1;
[r, c] = size(prmPAControl.List_Index_rxvga_to_test);
type_List_Index_rxvga_to_test = coder.newtype('double',[r c]);

%
prmPAControl.Number_rxvga_to_test = length(List_rxvga_to_test);
[r, c] = size(prmPAControl.Number_rxvga_to_test);
type_Number_rxvga_to_test = coder.newtype('double',[r c]);

%% Program ID
prmPAControl.Program_ID = Program_ID;
[r, c] = size(prmPAControl.Program_ID);
type_Program_ID = coder.newtype('double',[r c]);

%% Number of state to test
switch prmPAControl.Program_ID
    case 0
        prmPAControl.Number_State_To_Test = 1;
    case 1
        prmPAControl.Number_State_To_Test = prmPAControl.Number_rfvga_to_test*prmPAControl.Number_rxvga_to_test;
    case {2, 3, 2.1, 9, 9.1, 9.2}
        prmPAControl.Number_State_To_Test = prmPAControl.Number_antenna_to_calibrate*prmPAControl.Number_Phase_Shifter_State;
    case {4, 5, 4.1, 8}
        prmPAControl.Number_State_To_Test = prmPAControl.Number_CodeBook_Tx_to_train;
    case {6, 8.1}
        prmPAControl.Number_State_To_Test = prmPAControl.Number_CodeBook_Tx_to_train;
    case 7
        prmPAControl.Number_State_To_Test = prmPAControl.Number_CodeBook_Tx_to_train;
end
[r, c] = size(prmPAControl.Number_State_To_Test);
type_Number_State_To_Test = coder.newtype('double',[r c]);

%% Once the below variables are changed, recompliation is required
% Receiver Side
prmPAControl.Maximum_Number_State_To_Test = 44;
[r, c] = size(prmPAControl.Maximum_Number_State_To_Test);
type_Maximum_Number_State_To_Test = coder.newtype('double',[r c]);

prmPAControl.Number_of_Frames_Per_Test_Tx = prmQPSK.RxBufferedFrames; % Key parameters which impacts the collected data organizations
%prmPAControl.Number_of_Frames_Per_Test_Tx = 800;
[r, c] = size(prmPAControl.Number_of_Frames_Per_Test_Tx);
type_Number_of_Frames_Per_Test_Tx = coder.newtype('double',[r c]);

prmPAControl.Number_of_Frames_Per_Test_Rx_Buffered = prmPAControl.Number_of_Frames_Per_Test_Tx/prmQPSK.RxBufferedFrames;
%prmPAControl.Number_of_Frames_Per_Test_Rx_Buffered = prmPAControl.Number_of_Frames_Per_Test_Tx*prmPAControl.Number_CodeBook_Rx_to_train;
[r, c] = size(prmPAControl.Number_of_Frames_Per_Test_Rx_Buffered);
type_Number_of_Frames_Per_Test_Rx_Buffered = coder.newtype('double',[r c]);

prmPAControl.Remove_Percent_Beginning_Per_Test = 0.4;
[r, c] = size(prmPAControl.Remove_Percent_Beginning_Per_Test);
type_Remove_Percent_Beginning_Per_Test = coder.newtype('double',[r c]);

prmPAControl.Remove_Percent_Ending_Per_Test = prmPAControl.Remove_Percent_Beginning_Per_Test;
[r, c] = size(prmPAControl.Remove_Percent_Ending_Per_Test);
type_Remove_Percent_Ending_Per_Test = coder.newtype('double',[r c]);

%% Configuration of transmitter 
Repeat = 10000;
switch prmPAControl.Program_ID
    case {2, 2.1}
        prmPAControl.Total_Frame_To_Transmit = ...
            prmPAControl.Number_antenna_to_calibrate*prmPAControl.Number_Phase_Shifter_State*prmPAControl.Number_of_Frames_Per_Test_Tx*Repeat;
    case {4, 4.1, 6}
        prmPAControl.Total_Frame_To_Transmit = ...
            prmPAControl.Number_CodeBook_Tx*prmPAControl.Number_of_Frames_Per_Test_Tx*Repeat;
    case {0, 1, 3, 5, 7, 8, 8.1, 9, 9.1, 9.2}
        prmPAControl.Total_Frame_To_Transmit = 50000000;
        %prmPAControl.Total_Frame_To_Transmit = 500;
end
[r, c] = size(prmPAControl.Total_Frame_To_Transmit);
type_Total_Frame_To_Transmit = coder.newtype('double',[r c]);

prmPAControl.StopTransmissionTime = prmPAControl.Total_Frame_To_Transmit * prmQPSK.FrameTime;
[r, c] = size(prmPAControl.StopTransmissionTime);
type_StopTransmissionTime = coder.newtype('double',[r c]);

%% Once the below variable is changed, recompliation is required
% Number of frame for different types
prmPAControl.Number_kinds_of_Offline_Frames_Normal = prmPAControl.Maximum_Number_State_To_Test;
[r, c] = size(prmPAControl.Number_kinds_of_Offline_Frames_Normal);
type_Number_kinds_of_Offline_Frames_Normal = coder.newtype('double',[r c]);

prmPAControl.Number_kinds_of_Offline_Frames_Calibration = prmPAControl.Number_Phase_Shifter_State*prmPAControl.Number_antenna_to_calibrate;
[r, c] = size(prmPAControl.Number_kinds_of_Offline_Frames_Calibration);
type_Number_kinds_of_Offline_Frames_Calibration = coder.newtype('double',[r c]);

prmPAControl.Number_kinds_of_Offline_Frames_Beamtraining = prmPAControl.Number_CodeBook_Tx_to_train;
[r, c] = size(prmPAControl.Number_kinds_of_Offline_Frames_Beamtraining);
type_Number_kinds_of_Offline_Frames_Beamtraining = coder.newtype('double',[r c]);

%% Receiving power configurations
%
prmPAControl.rx_rxvga_to_use = 5;
[r, c] = size(prmPAControl.rx_rxvga_to_use);
type_rx_rxvga_to_use = coder.newtype('double',[r c]);
%
prmPAControl.rx_rfvga_to_use = 2;
[r, c] = size(prmPAControl.rx_rfvga_to_use);
type_rx_rfvga_to_use = coder.newtype('double',[r c]);

%% Indicaor of testing receiver beams
prmPAControl.nthRxBeam = 1;
[r, c] = size(prmPAControl.nthRxBeam);
type_nthRxBeam = coder.newtype('double',[r c]);

%% Indicaor of testing transmitter beams
prmPAControl.nthTxBeam = 1;
[r, c] = size(prmPAControl.nthTxBeam);
type_nthTxBeam = coder.newtype('double',[r c]);

%% Antenna activation setting for fine phased error calibration
prmPAControl.ControlTx = -1;
[r, c] = size(prmPAControl.ControlTx);
type_ControlTx = coder.newtype('double',[r c]);
prmPAControl.Activate_Tx_ID = -1;
[r, c] = size(prmPAControl.Activate_Tx_ID);
type_Activate_Tx_ID = coder.newtype('double',[r c]);
prmPAControl.Activate_Tx_Phase_ID = -1;
[r, c] = size(prmPAControl.Activate_Tx_Phase_ID);
type_Activate_Tx_Phase_ID = coder.newtype('double',[r c]);
prmPAControl.Activate_Rx_ID = -1;
[r, c] = size(prmPAControl.Activate_Rx_ID);
type_Activate_Rx_ID = coder.newtype('double',[r c]);
prmPAControl.Activate_Rx_Phase_ID = -1;
[r, c] = size(prmPAControl.Activate_Rx_Phase_ID);
type_Activate_Rx_Phase_ID = coder.newtype('double',[r c]);

%% Structure of prmPAControl and type of prmPAControl
type_prmPAControl = coder.newtype('struct',struct(...
                                 'SerialPortTx', type_SerialPortTx,...
                                 'SerialPortRx', type_SerialPortRx,...
                                 'Holding_Time_Open_Serial', type_Holding_Time_Open_Serial,...
                                 'Holding_Time_For_Writing_Buffer', type_Holding_Time_For_Writing_Buffer,...
                                 'Holding_Time_Turn_On_Antenna', type_Holding_Time_Turn_On_Antenna,...
                                 'Holding_Time_Change_Phase_Shifter_State', type_Holding_Time_Change_Phase_Shifter_State,...
                                 'Holding_Time_Self_Defined', type_Holding_Time_Self_Defined,...
                                 'Number_Phase_Shifter_State', type_Number_Phase_Shifter_State,...
                                 'Number_antenna', type_Number_antenna,...
                                 'txoffall', type_txoffall,...
                                 'txonall', type_txonall,...
                                 'rxoffall', type_rxoffall,...
                                 'rxonall', type_rxonall,...
                                 'txon_single', type_txon_single,...
                                 'rxon_single', type_rxon_single,...
                                 'txon', type_txon,...
                                 'rxon', type_rxon,...
                                 'phvec_Antenna1', type_phvec_Antenna1, ...
                                 'phvec', type_phvec,...
                                 'Number_phvec', type_Number_phvec,...
                                 'List_Index_antenna_to_calibrate', type_List_Index_antenna_to_calibrate,...
                                 'Number_antenna_to_calibrate', type_Number_antenna_to_calibrate,...
                                 'beam_label', type_beam_label, ...
                                 'CodeBook_Tx', type_CodeBook_Tx,...
                                 'Number_CodeBook_Tx', type_Number_CodeBook_Tx,...
                                 'List_Index_CodeBook_Tx_to_train', type_List_Index_CodeBook_Tx_to_train,...
                                 'Number_CodeBook_Tx_to_train', type_Number_CodeBook_Tx_to_train,...
                                 'CodeBook_Rx', type_CodeBook_Rx,...
                                 'Number_CodeBook_Rx', type_Number_CodeBook_Rx,...
                                 'List_Index_CodeBook_Rx_to_train', type_List_Index_CodeBook_Rx_to_train,...
                                 'Number_CodeBook_Rx_to_train', type_Number_CodeBook_Rx_to_train,...
                                 'rfvga', type_rfvga,...
                                 'rxvga', type_rxvga,...
                                 'List_Index_rfvga_to_test', type_List_Index_rfvga_to_test,...
                                 'Number_rfvga_to_test', type_Number_rfvga_to_test,...
                                 'List_Index_rxvga_to_test', type_List_Index_rxvga_to_test,...
                                 'Number_rxvga_to_test', type_Number_rxvga_to_test,...
                                 'Program_ID', type_Program_ID,...
                                 'Number_State_To_Test', type_Number_State_To_Test,...
                                 'Maximum_Number_State_To_Test', type_Maximum_Number_State_To_Test,...
                                 'Number_of_Frames_Per_Test_Tx', type_Number_of_Frames_Per_Test_Tx,...
                                 'Number_of_Frames_Per_Test_Rx_Buffered', type_Number_of_Frames_Per_Test_Rx_Buffered,...
                                 'Remove_Percent_Beginning_Per_Test', type_Remove_Percent_Beginning_Per_Test,...
                                 'Remove_Percent_Ending_Per_Test', type_Remove_Percent_Ending_Per_Test,...
                                 'Total_Frame_To_Transmit', type_Total_Frame_To_Transmit,...
                                 'StopTransmissionTime', type_StopTransmissionTime,...
                                 'Number_kinds_of_Offline_Frames_Normal', type_Number_kinds_of_Offline_Frames_Normal,...
                                 'Number_kinds_of_Offline_Frames_Calibration', type_Number_kinds_of_Offline_Frames_Calibration,...
                                 'Number_kinds_of_Offline_Frames_Beamtraining', type_Number_kinds_of_Offline_Frames_Beamtraining,...
                                 'rx_rxvga_to_use', type_rx_rxvga_to_use,...
                                 'rx_rfvga_to_use', type_rx_rfvga_to_use,...
                                 'nthRxBeam', type_nthRxBeam,...
                                 'nthTxBeam', type_nthTxBeam,...
                                 'ControlTx', type_ControlTx,...
                                 'Activate_Tx_ID', type_Activate_Tx_ID,...
                                 'Activate_Tx_Phase_ID', type_Activate_Tx_Phase_ID,...
                                 'Activate_Rx_ID', type_Activate_Rx_ID,...
                                 'Activate_Rx_Phase_ID', type_Activate_Rx_Phase_ID));                      
