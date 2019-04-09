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
% This main script runs the receiver which includes the process of Rx USRP 
% configuration, Rx phased array configuration, data reception and 
% data arrangement.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Add paths
folder_name = ...
[
"/Experimental_Validation";
];
for i=1:length(folder_name)
    Lib_path = char(strcat(pwd,folder_name(i)));
    addpath(genpath(Lib_path));
end


%% Connect to the USRP
clear all
address = '192.168.20.2';
platform = 'N200/N210/USRP2';
SerialPortTx = ['/dev/tty.usbmodem2441']; % USB serial port for Tx phased array, which need to be checked before experiments.
SerialPortRx = ['/dev/tty.usbmodem321'];  % USB serial port for Rx phased array, which need to be checked before experiments.
Program_ID = 6;
Detect_Noise_Floor = 0; % This is set to 1 to estimate the noise floor for Program_ID equals 9/9.1/9.2

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


%% Initialization
% Receiver parameter structure
prmQPSKReceiver = sdruqpskreceiver_init(platform);
prmQPSKReceiver.Platform = platform;
prmQPSKReceiver.Address = address;

% Phase array control parameter structure
[prmPAControl, type_prmPAControl] = phase_array_control_init(prmQPSKReceiver, SerialPortTx, SerialPortRx, Program_ID);
type_portion_of_boundary_frame_to_remove = coder.newtype('double',[1 1]);

% Calibration control parameter
[Calibration, type_Calibration] = Calibration_control_init();

% Code generation usage
compileIt  = true; % true if code is to be compiled for accelerated execution
useCodegen = true; % true to run the latest generated code (mex file) instead of MATLAB code                             


%% Global storage for corrupted signal and calculated BERMER
global Raw_corruptSignal
SamplesPerFrame = prmQPSKReceiver.FrameSize * prmQPSKReceiver.Upsampling * prmQPSKReceiver.RxBufferedFrames;
Raw_corruptSignal = complex(zeros(...
    SamplesPerFrame,...
    int32(prmPAControl.Number_of_Frames_Per_Test_Rx_Buffered),...
    prmPAControl.Maximum_Number_State_To_Test));

global BERMER_Non_Count
BERMER_Non_Count = zeros(7, ...
    prmQPSKReceiver.RxBufferedFrames, ...
    prmPAControl.Maximum_Number_State_To_Test*prmPAControl.Number_of_Frames_Per_Test_Rx_Buffered);

global Baseband_y
Baseband_y = complex(zeros(prmQPSKReceiver.FrameSize, ...
    prmQPSKReceiver.RxBufferedFrames, ...
    prmPAControl.Maximum_Number_State_To_Test*prmPAControl.Number_of_Frames_Per_Test_Rx_Buffered));

global FreqOffset
FreqOffset = complex(ones(1,prmPAControl.Maximum_Number_State_To_Test*prmPAControl.Number_of_Frames_Per_Test_Rx_Buffered));


%% Compilation
if compileIt
    codegen('run_SDRuQPSKReceiver', 'PATR.c', 'PATR.h', '-args', {coder.Constant(prmQPSKReceiver), type_prmPAControl}, ...
                                                        '-I', char(strcat(pwd,'/Experimental_Validation/src/array_control')));%#ok<UNRCH>
    codegen('Control_Phase_Array_External', 'PATR.c', 'PATR.h', '-args', {type_prmPAControl}, ...
                                                        '-I', char(strcat(pwd,'/Experimental_Validation/src/array_control')));%#ok<UNRCH>
    codegen('Control_Phase_Array_External_Phase_Error_Calibration', 'PATR.c', 'PATR.h', '-args', {type_prmPAControl,type_Calibration},...
                                                        '-I', char(strcat(pwd,'/Experimental_Validation/src/array_control')));%#ok<UNRCH>
    codegen('run_Offline_SDRuQPSKDecoder', '-args', {coder.Constant(prmQPSKReceiver), type_prmPAControl, type_portion_of_boundary_frame_to_remove});%#ok<UNRCH>
end


%% Clearance
if useCodegen
    clear run_SDRuQPSKTransmitter_mex %#ok<UNRCH>
    clear run_SDRuQPSKReceiver_mex %#ok<UNRCH>
    clear Control_Phase_Array_External_mex %#ok<UNRCH>
    clear Control_Phase_Array_External_Phase_Error_Calibration_mex %#ok<UNRCH>
end


%% Program_ID 1-8
if prmPAControl.Program_ID < 9 && prmPAControl.Program_ID > 0
    Number_Test = 1;
    for nthTest = 1 : Number_Test
        disp('-------------------------------------------------------------')
        fprintf('-----------%d-th collection of received raw signal------------\n',nthTest);
        disp('-------------------------------------------------------------')
        if prmPAControl.Program_ID ~=6 && prmPAControl.Program_ID ~=8.1
            if useCodegen
               clear run_SDRuQPSKReceiver_mex %#ok<UNRCH>
               run_SDRuQPSKReceiver_mex(prmQPSKReceiver, prmPAControl);
            else
               run_SDRuQPSKReceiver(prmQPSKReceiver, prmPAControl); 
            end
            fprintf('\n-----------Saving data for %d test-----------\n', nthTest);
            save(['Raw_corruptSignal' num2str(nthTest) '.mat'],'Raw_corruptSignal');
            fprintf('\n-----------Saved data for %d test-----------\n', nthTest);
        elseif prmPAControl.Program_ID == 6
            for nthRxBeam = 1:1:11
                prmPAControl.nthRxBeam = nthRxBeam;
                if useCodegen
                   clear run_SDRuQPSKReceiver_mex %#ok<UNRCH>
                   run_SDRuQPSKReceiver_mex(prmQPSKReceiver, prmPAControl);
                else
                   run_SDRuQPSKReceiver(prmQPSKReceiver, prmPAControl); 
                end
                fprintf('\n-----------Saving data (%d receiver beaming pattern applied)-----------  \n', nthRxBeam);
                save(['Raw_corruptSignal' num2str(nthTest) '_' num2str(nthRxBeam) '.mat'],'Raw_corruptSignal');
                fprintf('\n-----------Saved data (%d receiver beaming pattern applied)----------- %d \n', nthRxBeam)
            end
        else 
            for nthTxBeam = 1:1:11
                prmPAControl.nthTxBeam = nthTxBeam;
                if useCodegen
                   clear Control_Phase_Array_External_mex %#ok<UNRCH>
                   Control_Phase_Array_External_mex(prmPAControl);
                   clear run_SDRuQPSKReceiver_mex %#ok<UNRCH>
                   run_SDRuQPSKReceiver_mex(prmQPSKReceiver, prmPAControl);
                else
                   run_SDRuQPSKReceiver(prmQPSKReceiver, prmPAControl); 
                end
                fprintf('\n-----------Saving data (%d-th transmitter beaming pattern applied)-----------  \n', nthTxBeam);
                save(['Raw_corruptSignal' num2str(nthTest) '_' num2str(nthTxBeam) '.mat'],'Raw_corruptSignal');
                fprintf('\n-----------Saved data (%d-th transmitter beaming pattern applied)----------- %d \n', nthTxBeam)
            end
        end
    end
end


%% Element gain estimation and PMBC data collection
Fixed_Side_Antenna_ID_Range = [5];
Fixed_Side_Antenna_Phase_Range = [1];

% Element gain measurement with Program_ID equaling 9.
if prmPAControl.Program_ID == 9
   ind = 0;
   for Activate_Tx_ID = 1:1:12
       for Activate_Tx_Phase_ID = 1:1:4
           for Activate_Rx_ID = 1:1:12
               for Activate_Rx_Phase_ID = 1:1:4
                   if (ismember(Activate_Tx_ID,Fixed_Side_Antenna_ID_Range) && ismember(Activate_Tx_Phase_ID,Fixed_Side_Antenna_Phase_Range)) ...
                           || (ismember(Activate_Rx_ID,Fixed_Side_Antenna_ID_Range) && ismember(Activate_Rx_Phase_ID,Fixed_Side_Antenna_Phase_Range) )
                       ind = ind + 1;
                       if ind >= 1
                           % Set phased array  
                           Configure_Phase_Array_Setting_Standard(prmPAControl,Activate_Tx_ID,Activate_Tx_Phase_ID,...
                                                        Activate_Rx_ID,Activate_Rx_Phase_ID);

                           % Run receiver
                           clear run_SDRuQPSKReceiver_mex %#ok<UNRCH>
                           run_SDRuQPSKReceiver_mex(prmQPSKReceiver, prmPAControl);

                           % Save data
                           fprintf('\n-----------Saving data (%d )-----------  \n', ind);
                           if Detect_Noise_Floor == 1
                               save(['Noise_Raw_corruptSignal_TX' num2str(Activate_Tx_ID) '_'...
                                                            num2str(Activate_Tx_Phase_ID) '_Rx' ...
                                                            num2str(Activate_Rx_ID) '_' ...
                                                            num2str(Activate_Rx_Phase_ID) '.mat'],'Raw_corruptSignal');                            
                           else
                               save(['Raw_corruptSignal_TX' num2str(Activate_Tx_ID) '_'...
                                                            num2str(Activate_Tx_Phase_ID) '_Rx' ...
                                                            num2str(Activate_Rx_ID) '_' ...
                                                            num2str(Activate_Rx_Phase_ID) '.mat'],'Raw_corruptSignal');
                           end
                           fprintf('\n-----------Saved data (%d )----------- %d \n', ind);
                       end
                       % Restart tx program
                       if mod(ind,40) == 0 && Detect_Noise_Floor ~= 1 
                           fprintf('\n---------------------------------------------\n');
                           fprintf('-----Please restart transmitting program-----');
                           fprintf('\n---------------------------------------------\n');
                           pause;
                       end
                   end
               end
           end
       end

   end
end
   
% PMBC data collection with Program_ID equaling 9.1 and 9.2.
Ant_Ref_Ind_Range = [1];
Ant_Ref_Phase_Range = [1 2 3 4];
if prmPAControl.Program_ID == 9.1 || prmPAControl.Program_ID == 9.2
   ind = 0;
   for Fixed_Antenna_ID = Fixed_Side_Antenna_ID_Range
       for Fixed_Antenna_Phase = Fixed_Side_Antenna_Phase_Range
           for Ant_Ref_Ind = Ant_Ref_Ind_Range
               for Ant_Ref_Ind_Phase = Ant_Ref_Phase_Range
                   for Ant_To_Be_Calibrated_Ind = 2:12
                       fprintf('\nCalibration of %d-th antenna\n', Ant_To_Be_Calibrated_Ind)                  
                       for Ant_To_Be_Calibrated_Phase = 1:1:4
                           ind = ind + 1;
                           if ind >= 1  
                               % Set phased array            
                                Configure_Phase_Array_Setting(prmPAControl,...
                                                        Fixed_Antenna_ID,...
                                                        Fixed_Antenna_Phase,...
                                                        Ant_Ref_Ind,...
                                                        Ant_Ref_Ind_Phase,...
                                                        Ant_To_Be_Calibrated_Ind,...
                                                        Ant_To_Be_Calibrated_Phase)             
                               % Run receiver
                               clear run_SDRuQPSKReceiver_mex %#ok<UNRCH>
                               run_SDRuQPSKReceiver_mex(prmQPSKReceiver, prmPAControl);       
                               % Save data
                               fprintf('\n-----------Saving data (%d )-----------  \n', ind);
                               if Detect_Noise_Floor == 1 
                                   save(['Noise_Raw_corruptSignal_To_Cal_' num2str(Ant_To_Be_Calibrated_Ind) '_'...
                                                                num2str(Ant_To_Be_Calibrated_Phase) '_Ref' ...
                                                                num2str(Ant_Ref_Ind) '_' ...
                                                                num2str(Ant_Ref_Ind_Phase) '_Fix' ...
                                                                num2str(Fixed_Antenna_ID) '_' ...
                                                                num2str(Fixed_Antenna_Phase) ...
                                                                num2str(prmPAControl.Program_ID) '.mat'],'Raw_corruptSignal');                          
                               else
                                   save(['Raw_corruptSignal_To_Cal_' num2str(Ant_To_Be_Calibrated_Ind) '_'...
                                                                num2str(Ant_To_Be_Calibrated_Phase) '_Ref' ...
                                                                num2str(Ant_Ref_Ind) '_' ...
                                                                num2str(Ant_Ref_Ind_Phase) '_Fix' ...
                                                                num2str(Fixed_Antenna_ID) '_' ...
                                                                num2str(Fixed_Antenna_Phase) ...
                                                                num2str(prmPAControl.Program_ID) '.mat'],'Raw_corruptSignal'); 
                               end
                               fprintf('\n-----------Saved data (%d )----------- %d \n', ind);       
                           end
                           % Restart tx program
                           if mod(ind,60) == 0 && Detect_Noise_Floor ~= 1 
                               fprintf('\n---------------------------------------------\n');
                               fprintf('-----Please restart transmitting program-----');
                               fprintf('\n---------------------------------------------\n');
                               pause;
                           end
                       end
                   end
               end
           end
       end
   end
end


%% Save configuration parameters
save('prmQPSKReceiver.mat','prmQPSKReceiver');
save('prmPAControl.mat','prmPAControl');

