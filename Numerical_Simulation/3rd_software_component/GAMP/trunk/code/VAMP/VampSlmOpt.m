% Options class for VampSlmEst.m
%
% USAGE
% -----
% "opt = VampSlmOpt" uses default values
%
% "opt = VampSlmOpt('property1',value1,'property2',value2, ...)" uses default
%       values, except for the specified property/value pairs
%
% 

classdef VampSlmOpt

    properties
        nitMax = 50;    % maximum number of VAMP iterations
        tol = 1e-4;     % stopping tolerance
        gamMin = 1e-11; % minimum allowed precision 
        gamMax = 1e11;  % maximum allowed precision 
        damp = 0.97;    % damping parameter in (0,1] on r2 
        dampGam = [];   % damping parameter in (0,1] on gam2 
                        % default =damp if not specified

        Ah = [];        % fxn handle for A', needed only if A is a fxn handle
        N = [];         % =size(A,2), needed only if A is a fxn handle

        U = [];         % matrix of eigenvectors of A*A', or fxn handle
                        %   [U,D]=eig(A*A'); d=diag(D);
        d = [];         % vector of eigenvalues of A*A'
        Uh = [];        % fxn handle for U', needed only if U is a fxn handle

        UhA = [];       % fxn handle for U'*A
        AhU = [];       % fxn handle for A'*U, needed only if UhA is a fxn handle

        r1init = [];    % initial value of vector r1 (estimate of x)
        gam1init = [];  % initial value of scalar gam1 (precision on r1)
        r2init = [];    % initial value of vector r2 (estimate of x) ... if damping
        gam2init = [];  % initial value of scalar gam2 (precision on r2) ... if damping

        learnNoisePrec = true; % learn the noise precision?
        learnNoisePrecAlg = 'EM'; % learning type in {'EM','Newton'}
        learnNoisePrecNit = 100; % iterations used for noise learning
        learnNoisePrecTol = 1e-2; % tolerance used for noise learning 
        NoisePrecMax = 1e9; % max allowed value of noise precision 
        NoisePrecInit = [];% noise precision initialization (can leave empty)

        verbose = false;% verbose switch 
        silent = false; % disable all warnings 
        fxnErr = [];    % handle to a fxn of x2 for error reporting, e.g.,
                        %   fxnErr = @(x2) 10*log10( ...
                        %                    sum(abs(x2-xTrue).^2,1) ...
                        %                    ./sum(abs(xTrue).^2,1) );
        fxnGen = [];    % handle to a general fxn of the form
                        %   fxnGen = @(r1old,r1,gam1,x1,eta1,...
                        %              r2old,r2,gam2,x2,eta2)
        fxnStop = [];   % handle to a stopping fxn of the form
                        %   fxnStop = @(it,err,...
                        %               r1old,r1,gam1,x1,eta1,...
                        %               r2old,r2,gam2,x2,eta2) ...
                        %             all(err<-SNRdB);
        histIntvl = 1;  % can save memory by decimating the saved history

        divChange = 1e-3; % amount to perturb input for denoiser Monte-Carlo
                        % divergence estimate; passed to FxnhandleEstimIn 
                        % when the denoiser reports a single output
    end

    methods
        
        % Constructor with default options
        function opt = VampSlmOpt(varargin)
            if nargin == 0
                % No custom parameters values, thus create default object
                return
            elseif mod(nargin, 2) == 0
                % User is providing property/value pairs
                names = fieldnames(opt);    % Get names of class properties

                % Iterate through every property/value pair, assigning
                % user-specified values.  Note that the matching is NOT
                % case-sensitive
                for i = 1:2:nargin-1
                    if any(strcmpi(varargin{i}, names))
                        % User has specified a valid property
                        propName = names(strcmpi(varargin{i}, names));
                        opt.(propName{1}) = varargin{i+1};
                    else
                        % Unrecognized property name
                        error('VampSlmOpt: %s is an unrecognized option', ...
                            num2str(varargin{i}));
                    end
                end
                return
            else
                error(['The VampSlmOpt constructor requires arguments ' ...
                    'be provided in pairs, e.g., VampSlmOpt(''verbose'',' ...
                    ' false, ''nitMax'', 50)'])
            end
        end % VampSlmOpt

    end % methods

end
