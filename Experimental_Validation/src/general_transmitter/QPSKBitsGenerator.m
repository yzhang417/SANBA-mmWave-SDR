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
% This class generates bits stream for the transmitted data stream which 
% are y_Normal, y_Tx_Calibration, y_Tx_Beamtraining and 
% y_channel_estimation. They are the bit streams before modulation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef QPSKBitsGenerator < matlab.System
%#codegen
    % Generates the bits for each frame
    properties (Nontunable)
        MessageLength = 105;
        BernoulliLength = 69;
        ScramblerBase = 2;
        ScramblerPolynomial = [1 1 1 0 1];
        ScramblerInitialConditions = [0 0 0 0];
    end
    
    properties (Access=private)
        pHeader
        pScrambler
        pMsgStrSet_Normal
        pMsgStrSet_Tx_Calibration
        pMsgStrSet_Tx_Beamtraining
        pCount
    end
    
    methods
        function obj = QPSKBitsGenerator(varargin)
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods (Access=protected)
        function setupImpl(obj, ~)
            bbc = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1]; % Bipolar Barker Code
            ubc = ((bbc + 1) / 2)'; % Unipolar Barker Code
            temp = (repmat(ubc,1,2))';
            obj.pHeader = temp(:);
            obj.pCount = 0;
            obj.pScrambler = comm.Scrambler(obj.ScramblerBase, ...
            obj.ScramblerPolynomial, obj.ScramblerInitialConditions);
            %The following set could be passed by parameters.
            obj.pMsgStrSet_Normal = [...
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
          
            obj.pMsgStrSet_Tx_Calibration = [...
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
         
            obj.pMsgStrSet_Tx_Beamtraining = [...
              'Phvec state 001';'Phvec state 002';'Phvec state 003';'Phvec state 004';...
              'Phvec state 005';'Phvec state 006';'Phvec state 007';'Phvec state 008';...
              'Phvec state 009';'Phvec state 010';'Phvec state 011']; 
        end
        
        function [y_Normal, y_Tx_Calibration, y_Tx_Beamtraining, y_channel_estimation] = stepImpl(obj)
            % Converts the message string to bit format
            cycle = mod(obj.pCount,length(obj.pMsgStrSet_Normal(:,1)));
            msgStr = obj.pMsgStrSet_Normal(cycle+1,:);
            msgBin = de2bi(int8(msgStr),7,'left-msb');
            msg = reshape(double(msgBin).',obj.MessageLength,1);
            data = [msg ; randi([0 1], obj.BernoulliLength, 1)];
            
            % Scramble the data
            scrambledData = obj.pScrambler(data);
            
            % Append the scrambled bit sequence to the header
            y_Normal = [obj.pHeader ; scrambledData];

            % Converts the message string to bit format
            cycle = mod(obj.pCount,length(obj.pMsgStrSet_Tx_Calibration(:,1)));
            msgStr = obj.pMsgStrSet_Tx_Calibration(cycle+1,:);
            msgBin = de2bi(int8(msgStr),7,'left-msb');
            msg = reshape(double(msgBin).',obj.MessageLength,1);
            data = [msg ; randi([0 1], obj.BernoulliLength, 1)];
            
            % Scramble the data
            scrambledData = obj.pScrambler(data);
            
            % Append the scrambled bit sequence to the header
            y_Tx_Calibration = [obj.pHeader ; scrambledData];
            
            % Converts the message string to bit format
            cycle = mod(obj.pCount,length(obj.pMsgStrSet_Tx_Beamtraining(:,1)));
            msgStr = obj.pMsgStrSet_Tx_Beamtraining(cycle+1,:);
            msgBin = de2bi(int8(msgStr),7,'left-msb');
            msg = reshape(double(msgBin).',obj.MessageLength,1);
            data = [msg ; randi([0 1], obj.BernoulliLength, 1)];
            
            % Scramble the data
            scrambledData = obj.pScrambler(data);
            
            % Append the scrambled bit sequence to the header
            y_Tx_Beamtraining = [obj.pHeader ; scrambledData];
            
            % Message for channel estimation
            y_channel_estimation = [obj.pHeader ; ones(length(scrambledData),1)];
            
            obj.pCount = obj.pCount+1;            
        end
        
        function resetImpl(obj)
            obj.pCount = 0;
            reset(obj.pScrambler);
        end
        
        function releaseImpl(obj)
            release(obj.pScrambler);
        end
        
        function N = getNumInputsImpl(~)
            N = 0; 
        end
        
        function N = getNumOutputsImpl(~)
            N = 4;
        end
    end
end

