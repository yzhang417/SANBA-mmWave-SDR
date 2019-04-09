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
% Class description:
% This class is a MATLAB system object which serves as a receiver.
% The output of the key function of the class are given as below:
% BERMER: a structure array that groups the related evaluation metrics
% which includes 'BER' (1st-row), 'MER in dB' (4-th row), '95% MER in dB' 
% (5-th row), 'Received signal strength indication (RSSI)' (7-th row)
% BERMER_Non_Count: a structure array that groups the related evaluation
% metrics as BERMER array. However, BERMER_Non_Count does not perform the
% average operation, which allows BER and MER analysis for each received
% frame separately.
% Baseband_y: the demodulated raw corrupted signal.
% FreqOffset: the CFO estimated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
classdef sdruQPSKRx < matlab.System
%#codegen

    properties (Nontunable)
        DesiredAGCAmplitude
        DesiredAmplitude
        ModulationOrder
        DownsamplingFactor
        CoarseCompFrequencyResolution
        PhaseRecoveryLoopBandwidth
        PhaseRecoveryDampingFactor
        TimingRecoveryLoopBandwidth
        TimingRecoveryDampingFactor
        PostFilterOversampling
        PhaseErrorDetectorGain
        PhaseRecoveryGain
        TimingErrorDetectorGain
        TimingRecoveryGain
        FrameSize
        BarkerLength
        MessageLength
        SampleRate
        DataLength
        ReceiverFilterCoefficients
        DescramblerBase
        DescramblerPolynomial
        DescramblerInitialConditions
        PrintOption
    end
    
    properties (Access=private)
        pAGC
        pRxFilter
        pCoarseFreqCompensator
        pFineFreqCompensator
        pTimingRec
        pDataDecod
        pOldOutput % Stores the previous output of fine frequency compensation which is used by the same System object for phase error detection
    end
    
    methods
        function obj = sdruQPSKRx(varargin)
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods (Access=protected)
        function setupImpl(obj, ~, ~)
            obj.pAGC = comm.AGC;
            
            obj.pRxFilter = dsp.FIRDecimator(obj.DownsamplingFactor,obj.ReceiverFilterCoefficients);
            
            obj.pCoarseFreqCompensator = QPSKCoarseFrequencyCompensator(...
                'ModulationOrder', obj.ModulationOrder, ...
                'CoarseCompFrequencyResolution', obj.CoarseCompFrequencyResolution, ...
                'SampleRate', obj.SampleRate, ...
                'DownsamplingFactor', obj.DownsamplingFactor);
            
            % Refer C.57 to C.61 in Michael Rice's "Digital Communications - A Discrete-Time Approach" for K1 and K2
            theta = obj.PhaseRecoveryLoopBandwidth/(obj.PhaseRecoveryDampingFactor + 0.25/obj.PhaseRecoveryDampingFactor)/obj.PostFilterOversampling;
            d = 1 + 2*obj.PhaseRecoveryDampingFactor*theta + theta*theta;
            K1 = (4*obj.PhaseRecoveryDampingFactor*theta/d)/(obj.PhaseErrorDetectorGain*obj.PhaseRecoveryGain);
            K2 = (4*theta*theta/d)/(obj.PhaseErrorDetectorGain*obj.PhaseRecoveryGain);
            obj.pOldOutput = complex(0); % used to store past value
            obj.pFineFreqCompensator = QPSKFineFrequencyCompensator( ...
                'ProportionalGain', K1, ...
                'IntegratorGain', K2, ...
                'DigitalSynthesizerGain', -1*obj.PhaseRecoveryGain);
            
            % Refer C.57 to C.61 in Michael Rice's "Digital Communications - A Discrete-Time Approach" for K1 and K2
            theta = obj.TimingRecoveryLoopBandwidth/(obj.TimingRecoveryDampingFactor + 0.25/obj.TimingRecoveryDampingFactor)/obj.PostFilterOversampling;
            d = 1 + 2*obj.TimingRecoveryDampingFactor*theta + theta*theta;
            K1 = (4*obj.TimingRecoveryDampingFactor*theta/d)/(obj.TimingErrorDetectorGain*obj.TimingRecoveryGain);
            K2 = (4*theta*theta/d)/(obj.TimingErrorDetectorGain*obj.TimingRecoveryGain);           
            obj.pTimingRec = QPSKTimingRecovery('ProportionalGain', K1,...
                'IntegratorGain', K2, ...
                'PostFilterOversampling', obj.PostFilterOversampling, ...
                'BufferSize', obj.FrameSize);
                
            obj.pDataDecod = sdruQPSKDataDecoder('FrameSize', obj.FrameSize, ...
                'BarkerLength', obj.BarkerLength, ...
                'ModulationOrder', obj.ModulationOrder, ...
                'DataLength', obj.DataLength, ...
                'MessageLength', obj.MessageLength, ...
                'DescramblerBase', obj.DescramblerBase, ...
                'DescramblerPolynomial', obj.DescramblerPolynomial, ...
                'DescramblerInitialConditions', obj.DescramblerInitialConditions, ...
                'PrintOption', obj.PrintOption);   
        end
        
        function [BERMER, BERMER_Non_Count, Baseband_y, FreqOffset] =  stepImpl(obj, bufferSignal, Program_ID)
            
            % Apply automatic gain control to the signal
            % AGCSignal = obj.DesiredAmplitude*step(obj.pAGC, bufferSignal);
            AGCSignal = step(obj.pAGC, bufferSignal);
            %AGCSignal = bufferSignal;
            
            % Pass the signal through square root raised cosine received filter
            RCRxSignal = step(obj.pRxFilter,AGCSignal);
            
            % Coarsely compensate for the frequency offset
            [coarseCompSignal, FreqOffset]= step(obj.pCoarseFreqCompensator, RCRxSignal);
            
            % Buffers to store values required for plotting
            coarseCompBuffer = coder.nullcopy(complex(zeros(size(coarseCompSignal))));
            timingRecBuffer = coder.nullcopy(zeros(size(coarseCompSignal)));
            
            % Result vectors
            BERMER = zeros(6,1);
            BERMER_Non_Count = zeros(7,int32(length(coarseCompSignal)/200));
            Baseband_y = complex(zeros(100,int32(length(coarseCompSignal)/200)));
            j_Non_Count = 0;
            
            % Scalar processing for fine frequency compensation and timing recovery 
            for i=1:length(coarseCompSignal)             
                % Fine frequency compensation
                fineCompSignal = step(obj.pFineFreqCompensator, [obj.pOldOutput coarseCompSignal(i)]);
                %fineCompSignal = coarseCompSignal(i);
                coarseCompBuffer(i) = fineCompSignal;
                obj.pOldOutput = fineCompSignal;
                
                % Timing recovery of the received data
                [dataOut, isDataValid, timingRecBuffer(i)] = step(obj.pTimingRec, fineCompSignal);
                       
                if isDataValid
                    j_Non_Count = j_Non_Count + 1;
                    % Decoding the received data
                    if j_Non_Count > int32(length(coarseCompSignal)/200)
                        disp('Inconsistency of number of buffer frames in receiver occurs, it is ignorable.');
                        break;
                    end
                    [BERMER, BERMER_Non_Count(:,j_Non_Count), Baseband_y(:,j_Non_Count)] = step(obj.pDataDecod, dataOut, Program_ID); 
                end
            end   
            if j_Non_Count < int32(length(coarseCompSignal)/200)
                disp('Buffer frames less than expected.');
            end
        end
        
        function resetImpl(obj)
            reset(obj.pDataDecod);
        end
        
        function N = getNumInputsImpl(~)
            N = 2;
        end        
    end
end

