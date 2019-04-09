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
% This script is specifically for checking details of the performance of 
% the proposed calibration methods. It is a 
% simplified version of [Experimental_Validation/src/general_decoder
% /run_Offline_SDRuQPSKDecoder.m].
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


function BERMER = run_offline_decoder(prmQPSKReceiver, prmPAControl, portion_of_boundary_frame_to_remove)
%#codegen
global Raw_corruptSignal

%% Initialize variables
BERMER = zeros(7,prmPAControl.Number_State_To_Test);
BERMER_Temp = zeros(6,1);
coder.extrinsic('disp');

%% Decoding
disp('**********Start Decoding Message**********');
for i=1:1:prmPAControl.Number_State_To_Test
    %disp(i);
    Beginning_Frame = 1 + prmPAControl.Number_of_Frames_Per_Test_Rx_Buffered*portion_of_boundary_frame_to_remove;
    Ending_Frame = prmPAControl.Number_of_Frames_Per_Test_Rx_Buffered - ...
        prmPAControl.Number_of_Frames_Per_Test_Rx_Buffered*portion_of_boundary_frame_to_remove;
    % Calculate the Bit Error Rate (BER) and Modulation Error Rate (MER)
    % XXXXX
    BERMER(1:6,i) = BERMER_Temp;
    % Calculate the Power of Signal for Counted BER Setting
    rms2d = dsp.RMS('Dimension','All');    
    BERMER(7,i) = rms2d(Raw_corruptSignal(:,Beginning_Frame:Ending_Frame,i));
    release(rms2d);
end
disp('**********Stop Decoding Message**********');

