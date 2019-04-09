classdef QPSKFineFrequencyCompensator < matlab.System
%#codegen
%   This object is used only in supporting packages.
% 
%   Copyright 2012-2016 The MathWorks, Inc.

    properties (Nontunable)
        ProportionalGain = 0.008;
        IntegratorGain = 3e-5;
        DigitalSynthesizerGain = -1;
    end
    
    properties (Access=private)
        pPhase
        pLoopFilter
        pIntegrator
    end
    
    methods
        function obj = QPSKFineFrequencyCompensator(varargin)
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods (Access=protected)
        function setupImpl(obj, ~)
            obj.pPhase = 0;
            obj.pLoopFilter = dsp.IIRFilter( ...
                'Structure', 'Direct form II transposed', ...
                'Numerator', [1 0], 'Denominator', [1 -1]);
            obj.pIntegrator = dsp.IIRFilter(...
                'Structure', 'Direct form II transposed', ...
                'Numerator', [0 1], 'Denominator', [1 -1]);
        end
        
        function y = stepImpl(obj, u)
            
            y = u(2) * exp(1j*obj.pPhase);        % Complex phase shift          
            phErr = sign(real(u(1)))*imag(u(1)) -...  % Find phase error
                sign(imag(u(1)))*real(u(1));           
            loopFiltOut = obj.pLoopFilter(...     % Loop Filter
                phErr*obj.IntegratorGain);            
            DDSOut = ...                              % Direct Digital Synthesizer
                obj.pIntegrator(phErr*obj.ProportionalGain + loopFiltOut);
            obj.pPhase =  obj.DigitalSynthesizerGain * DDSOut;

        end
        
        function resetImpl(obj)
            reset(obj.pLoopFilter);
            reset(obj.pIntegrator);
        end
        
        function releaseImpl(obj)
            release(obj.pLoopFilter);
            release(obj.pIntegrator);            
        end
        
    end
end

