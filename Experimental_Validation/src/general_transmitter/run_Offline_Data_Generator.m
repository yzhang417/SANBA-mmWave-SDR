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
% This function generates the transmitted data in an offline way. The
% function can be run by using MATLAB Codegen. The generated data are
% stored in the global variables Raw_originalSignal_Normal, 
% Raw_originalSignal_Tx_Calibration (used for calibration), 
% Raw_originalSignal_Tx_Beamtraining (used for beam training),
% Raw_originalSignal_Tx_CE (used for channel estimation) and 
% Baseband_y_original (baseband signal associated with 
% Raw_originalSignal_Tx_CE).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function run_Offline_Data_Generator(prmQPSKTransmitter, prmPAControl)
%#codegen
global Raw_originalSignal_Normal
global Raw_originalSignal_Tx_Calibration
global Raw_originalSignal_Tx_Beamtraining
global Raw_originalSignal_Tx_CE
global Baseband_y_original

%%
persistent hTx
if isempty(hTx)
    %% Initialize the components
     % Create and configure the transmitter system object
     hTx = QPSKTransmitter(...
         'UpsamplingFactor', prmQPSKTransmitter.Upsampling, ...
         'MessageLength', prmQPSKTransmitter.MessageLength, ...
         'TransmitterFilterCoefficients',prmQPSKTransmitter.TransmitterFilterCoefficients, ...
         'DataLength', prmQPSKTransmitter.DataLength, ...
         'ScramblerBase', prmQPSKTransmitter.ScramblerBase, ...
         'ScramblerPolynomial', prmQPSKTransmitter.ScramblerPolynomial, ...
         'ScramblerInitialConditions', prmQPSKTransmitter.ScramblerInitialConditions);
end    

%% Offline data generation process
coder.extrinsic('disp')
disp('--------------------------------------------------------');
disp('********** Start Generating Transmit Data Set **********')
disp('--------------------------------------------------------');
NthFrame = 1;
while NthFrame <= prmPAControl.Maximum_Number_State_To_Test
    % Bit generation, modulation and transmission filtering
    [transmittedSignal_Normal, ...
     transmittedSignal_Tx_Calibration, ...
     transmittedSignal_Tx_Beamtraining, ...
     transmittedSignal_Tx_CE, ...
     modulatedData_CE] = step(hTx);
    % Store the data
    if NthFrame <= 11
        Raw_originalSignal_Tx_Beamtraining(:, NthFrame) = transmittedSignal_Tx_Beamtraining;
    end
    Raw_originalSignal_Normal(:, NthFrame) = transmittedSignal_Normal;
    Raw_originalSignal_Tx_Calibration(:, NthFrame) = transmittedSignal_Tx_Calibration;
    Raw_originalSignal_Tx_CE(:, NthFrame) = transmittedSignal_Tx_CE;
    Baseband_y_original(:, NthFrame) = modulatedData_CE;
    % Generate next data frame
    NthFrame = NthFrame + 1;
end   
disp('-------------------------------------------------------');
disp('********** Stop Generating Transmit Data Set **********') 
disp('-------------------------------------------------------');

%% Release system objects
release(hTx);
