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
% This function decodes the recived signal in an offline mode.
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
% portion_of_boundary_frame_to_remove: proportion of the frames that are
% removed to ensure that the counted frame is correctly transmitted by the
% intended array configuration.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output and modified global variable:
% BERMER: a structure array that groups the related evaluation metrics
% which includes 'BER' (1st-row), 'MER in dB' (4-th row), '95% MER in dB' 
% (5-th row), 'Received signal strength indication (RSSI)' (7-th row)
% BERMER_Non_Count: a structure array that groups the related evaluation
% metrics as BERMER array. However, BERMER_Non_Count does not perform the
% average operation, which allows BER and MER analysis for each received
% frame separately.
% Raw_corruptSignal would be updated.
% Baseband_y would be updated as the demodulated corrupted signal.
% FreqOffset would be updated as the CFO estimated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function BERMER = run_Offline_SDRuQPSKDecoder(prmQPSKReceiver, prmPAControl, portion_of_boundary_frame_to_remove)
%#codegen
global Raw_corruptSignal
global BERMER_Non_Count
global Baseband_y
global FreqOffset
persistent hRx
BERMER_Non_Count(:,:,:) = 0;
Baseband_y(:,:,:) = 0;

%% Initialize variables
BERMER = zeros(7,prmPAControl.Number_State_To_Test);
BERMER_Temp = zeros(6,1);
BERMER_Non_Count_Temp = zeros(7,prmQPSKReceiver.RxBufferedFrames);
Baseband_y_Temp = complex(zeros(prmQPSKReceiver.FrameSize,prmQPSKReceiver.RxBufferedFrames));
FreqOffset_Temp = complex(zeros(1,1));
coder.extrinsic('disp');

%% Create decoder object
if isempty(hRx)
    hRx = sdruQPSKRx( ...
        'DesiredAmplitude',               1/sqrt(prmQPSKReceiver.Upsampling), ...
        'ModulationOrder',                prmQPSKReceiver.M, ...
        'DownsamplingFactor',             prmQPSKReceiver.Downsampling, ...
        'CoarseCompFrequencyResolution',  prmQPSKReceiver.CoarseCompFrequencyResolution, ...
        'PhaseRecoveryLoopBandwidth',     prmQPSKReceiver.PhaseRecoveryLoopBandwidth, ...
        'PhaseRecoveryDampingFactor',     prmQPSKReceiver.PhaseRecoveryDampingFactor, ...
        'TimingRecoveryLoopBandwidth',    prmQPSKReceiver.TimingRecoveryLoopBandwidth, ...
        'TimingRecoveryDampingFactor',    prmQPSKReceiver.PhaseRecoveryDampingFactor, ...
        'PostFilterOversampling',         prmQPSKReceiver.Upsampling/prmQPSKReceiver.Downsampling, ...
        'PhaseErrorDetectorGain',         prmQPSKReceiver.PhaseErrorDetectorGain, ...
        'PhaseRecoveryGain',              prmQPSKReceiver.PhaseRecoveryGain, ...
        'TimingErrorDetectorGain',        prmQPSKReceiver.TimingErrorDetectorGain, ...
        'TimingRecoveryGain',             prmQPSKReceiver.TimingRecoveryGain, ...
        'FrameSize',                      prmQPSKReceiver.FrameSize, ...
        'BarkerLength',                   prmQPSKReceiver.BarkerLength, ...
        'MessageLength',                  prmQPSKReceiver.MessageLength, ...
        'SampleRate',                     prmQPSKReceiver.Fs, ...
        'DataLength',                     prmQPSKReceiver.DataLength, ...
        'ReceiverFilterCoefficients',     prmQPSKReceiver.ReceiverFilterCoefficients, ...
        'DescramblerBase',                prmQPSKReceiver.ScramblerBase, ...
        'DescramblerPolynomial',          prmQPSKReceiver.ScramblerPolynomial, ...
        'DescramblerInitialConditions',   prmQPSKReceiver.ScramblerInitialConditions,...
        'PrintOption',                    true);
end    

%% Decoding
disp('**********Start Decoding Message**********');
nth_received_frame = 1;
for i=1:1:prmPAControl.Number_State_To_Test
    %disp(i);
    Beginning_Frame = 1 + prmPAControl.Number_of_Frames_Per_Test_Rx_Buffered*portion_of_boundary_frame_to_remove;
    Ending_Frame = prmPAControl.Number_of_Frames_Per_Test_Rx_Buffered - ...
        prmPAControl.Number_of_Frames_Per_Test_Rx_Buffered*portion_of_boundary_frame_to_remove;
    NthFrame = Beginning_Frame;
    nth_received_frame = Beginning_Frame+(i-1)*prmPAControl.Number_of_Frames_Per_Test_Rx_Buffered;
    % Calculate the bit error rate (BER) and modulation error rate (MER)
    while NthFrame <= Ending_Frame
        % When the SDRu system object output is valid, decode the received message
        corruptSignal = Raw_corruptSignal(:,NthFrame,i);
        % Calculate the BER and MER
        [BERMER_Temp, BERMER_Non_Count_Temp, Baseband_y_Temp, FreqOffset_Temp] = step(hRx, corruptSignal, prmPAControl.Program_ID);
        % Update simulation time
        NthFrame = NthFrame+1;
        % Store BERMER_Non_Count
        BERMER_Non_Count(:,:,nth_received_frame) = BERMER_Non_Count_Temp;
        Baseband_y(:,:,nth_received_frame) = Baseband_y_Temp;
        FreqOffset(1,nth_received_frame) = FreqOffset_Temp;
        nth_received_frame = nth_received_frame + 1;
    end
    BERMER(1:6,i) = BERMER_Temp;
    reset(hRx);
    % Calculate the power of signal for counted BER setting
    rms2d = dsp.RMS('Dimension','All');    
    BERMER(7,i) = rms2d(Raw_corruptSignal(:,Beginning_Frame:Ending_Frame,i));
    release(rms2d);
end
disp('**********Stop Decoding Message**********');
release(hRx);

