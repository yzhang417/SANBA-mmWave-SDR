%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2019 Yi Zhang and The University of Texas at Austin 
% Copyright 2012-2016 The MathWorks, Inc.      
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
% Input arguments:
% prmQPSKReceiver: a structure array that groups the related 
% parameters on the QPSK receiver including frame format, modulation 
% scheme and etc. With the generated structure array, the baseband post 
% signal processing is set.
% prmPAControl: a structure array that groups the related parameters
% on phased array including available codebooks, receiving power gain,
% antenna switching on/off and etc. With prmPAControl, the RF configuration
% is set.
% portion_of_boundary_frame_to_remove: proportion of the frames that are
% removed to ensure that the counted frame is correctly transmitted by the
% intended array configuration.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Class description:
% This class generates a matlab system object which serves as a decoder of
% the receiver.
% Output of this class are given as below:
% BERMER: a structure array that groups the related evaluation metrics
% which includes 'BER' (1st-row), 'MER in dB' (4-th row), '95% MER in dB' 
% (5-th row), 'Received signal strength indication (RSSI)' (7-th row)
% BERMER_Non_Count: a structure array that groups the related evaluation
% metrics as BERMER array. However, BERMER_Non_Count does not perform the
% average operation, which allows BER and MER analysis for each received
% frame separately.
% Baseband_y: the demodulated raw corrupted signal.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef sdruQPSKDataDecoder < matlab.System
%#codegen    
    properties (Nontunable)
        FrameSize
        BarkerLength
        ModulationOrder
        DataLength
        MessageLength
        DescramblerBase
        DescramblerPolynomial
        DescramblerInitialConditions
        PrintOption
    end
    
    properties (Access=private)
        pCount %useless
        pDelay
        pPhase
        pBuffer
        pModulator
        pModulatedHeader
        pCorrelator
        pQPSKDemodulator
        pDescrambler
        pBitGenerator
        pBitGeneratorSync
        pSyncFlag  %useless
        pSyncIndex %useless
        pPreviousEstimatedFrameIndex
        pBER
        pMER
        pBER_Non_Count
        pMER_Non_Count
    end
     
    methods
        function obj = sdruQPSKDataDecoder(varargin)
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods (Access=protected)
        function setupImpl(obj, ~, ~)
            [obj.pCount, obj.pDelay, obj.pPhase] = deal(0);
            obj.pPreviousEstimatedFrameIndex = 100;
            obj.pSyncIndex=0;
            obj.pSyncFlag = true;
            obj.pBuffer=dsp.Buffer(obj.FrameSize*2, obj.FrameSize);
            bbc = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1]; % Bipolar Barker Code
            ubc = ((bbc + 1) / 2)'; % Unipolar Barker Code
            header = (repmat(ubc,1,2))';
            header = header(:);
            obj.pModulator = comm.QPSKModulator('BitInput', true, 'PhaseOffset', pi/4);
            obj.pModulatedHeader = step(obj.pModulator, header); % Modulate the header
            obj.pCorrelator = dsp.Crosscorrelator;
            obj.pQPSKDemodulator = comm.QPSKDemodulator('PhaseOffset',pi/4, 'BitOutput', true);
            obj.pDescrambler = comm.Descrambler(obj.DescramblerBase, obj.DescramblerPolynomial, obj.DescramblerInitialConditions);
            % BER and MER calculator
            obj.pBER = comm.ErrorRate;
            obj.pBER_Non_Count = comm.ErrorRate('ResetInputPort',true);
            RefConst = [   0.7071 + 0.7071i;  -0.7071 + 0.7071i;  -0.7071 - 0.7071i;   0.7071 - 0.7071i];
            obj.pMER = comm.MER('ReferenceSignalSource', 'Estimated from reference constellation', ...
                                                     'ReferenceConstellation',  RefConst, ...
                                                     'XPercentileMEROutputPort', true, ...
                                                     'XPercentileValue', 95,...
                                                     'MinimumMEROutputPort', true, ...
                                                     'SymbolCountOutputPort', true);  
            obj.pMER_Non_Count = comm.MER('ReferenceSignalSource', 'Estimated from reference constellation', ...
                                                     'ReferenceConstellation',  RefConst, ...
                                                     'XPercentileMEROutputPort', true, ...
                                                     'XPercentileValue', 95,...
                                                     'MinimumMEROutputPort', true, ...
                                                     'SymbolCountOutputPort', false);                                                 
        end
        
        function [BERMER, BERMER_Non_Count, Baseband_y] = stepImpl(obj, DataIn, Program_ID)
            % Buffer one frame in case that contiguous data scatter across two adjacent frames
            rxData = step(obj.pBuffer,DataIn);
            
            % Get a frame of data aligned on the frame boundary
            Data = rxData(obj.pDelay+1:obj.pDelay+length(rxData)/2);
            
            % Recovered Baseband signal
            Baseband_y = Data;
            
            % Phase estimation
            y = mean(conj(obj.pModulatedHeader) .* Data(1:obj.BarkerLength));
            
            % Compensating for the phase offset
            if Data(1)~=0
                phShiftedData = Data .* exp(-1j*obj.pPhase);
            else
                phShiftedData = complex(zeros(size(Data)));
            end
           
            % Calculate the MER (Modulation Error Ratio)       
            [MERdB, MinMER, PercentileMER, NumSym] = step(obj.pMER, phShiftedData);
            MER = [MERdB; MinMER; PercentileMER];
            [MERdB, MinMER, PercentileMER] = step(obj.pMER_Non_Count, phShiftedData);
            MER_Non_Count = [MERdB; MinMER; PercentileMER];
            
            % Demodulate the phase recovered data
            demodOut = step(obj.pQPSKDemodulator, phShiftedData);
            
            % Perform descrambling
            if Program_ID == 8 || Program_ID == 8.1
                deScrData = demodOut(obj.BarkerLength*log2(obj.ModulationOrder)+1 : obj.FrameSize*log2(obj.ModulationOrder));
            else    
                deScrData = step(obj.pDescrambler, demodOut(obj.BarkerLength*log2(obj.ModulationOrder)+1 : obj.FrameSize*log2(obj.ModulationOrder)));
            end
            % Recovering the message from the data
            Received = deScrData(1:obj.MessageLength);
            
            % Finding the delay to achieve frame synchronization
            z=abs(step(obj.pCorrelator,obj.pModulatedHeader,DataIn));
            [~, ind] = max(z);
            obj.pDelay = mod(length(DataIn)-ind,(length(DataIn)-1));
            
            % Phase ambiguity correction
            obj.pPhase = round(angle(y)*2/pi)/2*pi;     
            
            % Print received frame and estimate the received frame index
            [estimatedFrameIndex,syncIndex]=bits2ASCII(obj,Received);
            obj.pSyncIndex = syncIndex;
                        
            % With the estimated frame index, estimate the transmitted message
            if Program_ID == 8 || Program_ID == 8.1
                transmittedMessage = ones(obj.MessageLength,1);
            else
                transmittedMessage = messEstimator(estimatedFrameIndex, Program_ID, obj);
            end
            % Update the PreviousEstimatedFrameIndex
            obj.pPreviousEstimatedFrameIndex = estimatedFrameIndex;
            
            % Calculate the BER
            BER = step(obj.pBER, transmittedMessage, Received);
            BER_Non_Count = step(obj.pBER_Non_Count, transmittedMessage, Received, 1);
                        
            %Combine result for output
            BERMER = [BER; MER];
            BERMER_Non_Count = [BER_Non_Count; MER_Non_Count; estimatedFrameIndex];
           
        end
        
        function resetImpl(obj)
            reset(obj.pBuffer);
            reset(obj.pBER);
            reset(obj.pMER);
            reset(obj.pBER_Non_Count);
            reset(obj.pMER_Non_Count);
        end
        
        function releaseImpl(obj)
            release(obj.pBuffer);
        end
        
        function N = getNumInputsImpl(~)
            N = 2;
        end
        
    end
    
    methods (Access=private)
        function [estimatedFrameIndex,syncIndex]=bits2ASCII(obj,u)
            coder.extrinsic('disp')
            
            % Convert binary-valued column vector to 7-bit decimal values.
            w = [64 32 16 8 4 2 1]; % binary digit weighting
            Nbits = numel(u);
            Ny = Nbits/7;
            y = zeros(1,Ny);
            % Obtain ASCII values of received frame
            for i = 0:Ny-1
                y(i+1) = w*u(7*i+(1:7));
            end
            
            % Display ASCII message to command window   
%             if(obj.PrintOption)
%                 disp(char(y));
%             end

            % Retrieve last 2 ASCII values
            decodedNumber=y(Ny-1:end);
            % Create lookup table of ASCII values and corresponding integer numbers 
            look_tab=zeros(2,10);
            look_tab(1,:)=0:9;
            look_tab(2,:)=48:57;
            % Initialize variables
            estimatedFrameIndex=100;
            syncIndex=0;
            onesPlace=0;
            tensPlace=0;
            dec_found=false;
            unity_found=false;
            
            % Index lookup table with decoded ASCII values
            % There are more efficient ways to perform vector indexing
            % using MATLAB functions like find(). However, to meet codegen
            % requirements, the usage of the four loop was necessary.
            
            for ii=1:10
                % Find the ones place in the lookup table
                if ( decodedNumber(1) == look_tab(2,ii) )
                    onesPlace=10*look_tab(1,ii);
                    dec_found=true;
                end
                % Find the tens place in the lookup table
                if ( decodedNumber(2) == look_tab(2,ii) )
                    tensPlace=look_tab(1,ii);
                    unity_found=true;
                end
            end
                        
            % Estimate the frame index
            if(dec_found && unity_found)
                estimatedFrameIndex=onesPlace+tensPlace;
                syncIndex=obj.pSyncIndex+1;
            end
        end
        
        function msg = messEstimator(ind, Program_ID, obj)     
            MsgStrSet_Normal = [...
              'Hello world 001';'Hello world 002';'Hello world 003';'Hello world 004';...
              'Hello world 005';'Hello world 006';'Hello world 007';'Hello world 008';...
              'Hello world 009';'Hello world 010';'Hello world 011';'Hello world 012';...
              'Hello world 013';'Hello world 014';'Hello world 015';'Hello world 016';...
              'Hello world 017';'Hello world 018';'Hello world 019';'Hello world 020';...
              'Hello world 021';'Hello world 022';'Hello world 023';'Hello world 024';...
              'Hello world 025';'Hello world 026';'Hello world 027';'Hello world 028';...
              'Hello world 029';'Hello world 030';'Hello world 031';'Hello world 032';...
              'Hello world 033';'Hello world 034';'Hello world 035';'Hello world 036';...
              'Hello world 037';'Hello world 038';'Hello world 039';'Hello world 040';...
              'Hello world 041';'Hello world 042';'Hello world 043';'Hello world 044']; 
          
            MsgStrSet_Tx_Calibration = [...
              'Phvec calib 001';'Phvec calib 002';'Phvec calib 003';'Phvec calib 004';...
              'Phvec calib 005';'Phvec calib 006';'Phvec calib 007';'Phvec calib 008';...
              'Phvec calib 009';'Phvec calib 010';'Phvec calib 011';'Phvec calib 012';...
              'Phvec calib 013';'Phvec calib 014';'Phvec calib 015';'Phvec calib 016';...
              'Phvec calib 017';'Phvec calib 018';'Phvec calib 019';'Phvec calib 020';...
              'Phvec calib 021';'Phvec calib 022';'Phvec calib 023';'Phvec calib 024';...
              'Phvec calib 025';'Phvec calib 026';'Phvec calib 027';'Phvec calib 028';...
              'Phvec calib 029';'Phvec calib 030';'Phvec calib 031';'Phvec calib 032';...
              'Phvec calib 033';'Phvec calib 034';'Phvec calib 035';'Phvec calib 036';...
              'Phvec calib 037';'Phvec calib 038';'Phvec calib 039';'Phvec calib 040';...
              'Phvec calib 041';'Phvec calib 042';'Phvec calib 043';'Phvec calib 044'];
          
            MsgStrSet_Tx_Beamtraining = [...
              'Phvec state 001';'Phvec state 002';'Phvec state 003';'Phvec state 004';...
              'Phvec state 005';'Phvec state 006';'Phvec state 007';'Phvec state 008';...
              'Phvec state 009';'Phvec state 010';'Phvec state 011';]; 
                        
            if Program_ID == 2
                cycle = mod(ind,44);
                if cycle == 0
                    cycle = 44;
                end
                msgStr = MsgStrSet_Tx_Calibration(cycle,:); 
            elseif Program_ID == 4 || Program_ID == 6
                cycle = mod(ind,11);
                if cycle == 0
                    cycle = 11;
                end
                msgStr = MsgStrSet_Tx_Beamtraining(cycle,:);                
            else
                cycle = mod(ind,44);
                if cycle == 0
                    cycle = 44;
                end
                msgStr = MsgStrSet_Normal(cycle,:);                   
            end
            
            msgBin = de2bi(int8(msgStr),7,'left-msb');
            msg = reshape(double(msgBin).',obj.MessageLength,1);
        end
    end
end

