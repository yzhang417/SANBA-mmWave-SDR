classdef QPSKCoarseFrequencyCompensator < matlab.System
%#codegen
%   This object is used only in supporting packages.
% 
%   Copyright 2012-2016 The MathWorks, Inc.

    properties (Nontunable)
        ModulationOrder = 4;
        CoarseCompFrequencyResolution = 50;        
        SampleRate = 200000;
        DownsamplingFactor = 2;
    end
    
    properties (Access=private)
        pPhaseFreqOffset
        pCoarseFreqEst
    end
    
    methods
        function obj = QPSKCoarseFrequencyCompensator(varargin)
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods (Access=protected)
        function setupImpl(obj, ~)
            currentSampleRate = obj.SampleRate/obj.DownsamplingFactor;
            obj.pPhaseFreqOffset = comm.PhaseFrequencyOffset(...
                'PhaseOffset', 0, ...
                'FrequencyOffsetSource', 'Input port' , ...
                'SampleRate', currentSampleRate);
            obj.pCoarseFreqEst = comm.PSKCoarseFrequencyEstimator( ...
                'ModulationOrder', obj.ModulationOrder, ...
                'FrequencyResolution', obj.CoarseCompFrequencyResolution, ...
                'SampleRate', currentSampleRate);
        end
        
        function [compensatedSignal, FreqOffset]= stepImpl(obj, filteredSignal)
            
            % Find the frequency used for correction (the negative of the
            % actual offset)
            FreqOffset = -obj.pCoarseFreqEst(filteredSignal);
          
            % Remove the frequency offset
            compensatedSignal = ...
                obj.pPhaseFreqOffset(filteredSignal,FreqOffset);
            
        end
        
        function resetImpl(obj)
            reset(obj.pPhaseFreqOffset);
            reset(obj.pCoarseFreqEst);
        end
        
        function releaseImpl(obj)
            release(obj.pPhaseFreqOffset);
            release(obj.pCoarseFreqEst);            
        end
    end
end

