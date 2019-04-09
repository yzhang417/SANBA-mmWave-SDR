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
% This class generates the transmitted baseband data stream which is going
% to be fed into the USRP hardware. They are transmittedSignal_Normal, 
% transmittedSignal_Tx_Calibration, transmittedSignal_Tx_Beamtraining, 
% and transmittedSignal_Tx_CE. Besides, modulatedData_CE is the baseband
% signal corresponds transmittedSignal_Tx_CE before pulse filtering.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef QPSKTransmitter < matlab.System  
%#codegen
    % Generates the QPSK signal to be transmitted    
    properties (Nontunable)
        UpsamplingFactor = 4;
        MessageLength = 105;
        DataLength = 174;
        TransmitterFilterCoefficients = 1;
        ScramblerBase = 2;
        ScramblerPolynomial = [1 1 1 0 1];
        ScramblerInitialConditions = [0 0 0 0];
    end
    
    properties (Access=private)
        pBitGenerator
        pQPSKModulator
        pQPSKModulator2
        pTransmitterFilter
        pTransmitterFilter2
    end
    
    methods
        function obj = QPSKTransmitter(varargin)
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods (Access=protected)
        function setupImpl(obj)
            obj.pBitGenerator = QPSKBitsGenerator(...
                'MessageLength', obj.MessageLength, ...
                'BernoulliLength', obj.DataLength-obj.MessageLength, ...
                'ScramblerBase', obj.ScramblerBase, ...
                'ScramblerPolynomial', obj.ScramblerPolynomial, ...
                'ScramblerInitialConditions', obj.ScramblerInitialConditions);
            obj.pQPSKModulator  = comm.QPSKModulator('BitInput', true, 'PhaseOffset', pi/4);
            obj.pQPSKModulator2  = comm.QPSKModulator('BitInput', true, 'PhaseOffset', pi/4);
            obj.pTransmitterFilter = dsp.FIRInterpolator(obj.UpsamplingFactor, obj.TransmitterFilterCoefficients);
            obj.pTransmitterFilter2 = dsp.FIRInterpolator(obj.UpsamplingFactor, obj.TransmitterFilterCoefficients);
        end
        
        function [transmittedSignal_Normal, transmittedSignal_Tx_Calibration, transmittedSignal_Tx_Beamtraining, transmittedSignal_Tx_CE, modulatedData_CE] = stepImpl(obj)    
            % Generates the data to be transmitted 
            [transmittedData_Normal, transmittedData_Tx_Calibration, transmittedData_Tx_Beamtraining, transmittedData_CE] = obj.pBitGenerator(); 
            % Modulates the bits into QPSK symbols 
            % Square root Raised Cosine Transmit Filter
            
            modulatedData_Normal = obj.pQPSKModulator(transmittedData_Normal); 
            transmittedSignal_Normal = obj.pTransmitterFilter(modulatedData_Normal);  
            
            modulatedData_Calibration = obj.pQPSKModulator(transmittedData_Tx_Calibration); 
            transmittedSignal_Tx_Calibration = obj.pTransmitterFilter(modulatedData_Calibration); 
            
            modulatedData_Beamtraining = obj.pQPSKModulator(transmittedData_Tx_Beamtraining);
            transmittedSignal_Tx_Beamtraining = obj.pTransmitterFilter(modulatedData_Beamtraining);
            
            modulatedData_CE = obj.pQPSKModulator2(transmittedData_CE);
            transmittedSignal_Tx_CE = obj.pTransmitterFilter2(modulatedData_CE);
        end
        
        function resetImpl(obj)
            reset(obj.pBitGenerator);
            reset(obj.pQPSKModulator);
            reset(obj.pTransmitterFilter);
            reset(obj.pTransmitterFilter2);
        end
        
        function releaseImpl(obj)
            release(obj.pBitGenerator);
            release(obj.pQPSKModulator );
            release(obj.pTransmitterFilter);
            release(obj.pTransmitterFilter2);
        end
        
        function N = getNumInputsImpl(~)
            N = 0;
        end
    end
end

