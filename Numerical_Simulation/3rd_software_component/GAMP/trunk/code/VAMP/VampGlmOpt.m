classdef VampGlmOpt
    % Options class for VampGlmEstim.m

    properties
        nitMax = 50;    % maximum number of VAMP iterations
        tol = 1e-4;     % stopping tolerance

        gamMin = 1e-8;  % minimum allowed precision [1e-8]
        gamMax = 1e11;  % maximum allowed precision [1e11]

        damp = 0.9;     % damping parameter in (0,1] 
        dampGam = [];   % damping parameter in (0,1] on precisions
        % dampConfig is a binary vector that enables/disables damping of...
        %    [z1,gam2z,p2,x1,gam2x,r2, z2,gam1z,p1,x2,gam1x,r1]
       %dampConfig = [0,1,0,1,0,0, 1,0,0,0,1,0]; % original, from VampGlmEst.m
       %dampConfig = [0,1,1,1,0,0, 0,0,0,0,1,1]; % best from dampTest
        dampConfig = [0,1,1,1,0,0, 0,0,0,0,0,1]; % best from dampTest
       %dampConfig = [0,1,0,0,1,1, 0,1,1,0,0,0]; % symmetric choice

        Ah = [];        % fxn handle for A', needed only if A is a fxn handle
        N = [];         % =size(A,2), needed only if A is a fxn handle

        U = [];         % matrix of eigenvectors of A*A', or fxn handle,
                        % used (if present) when M<=N
                        %   [U,D]=eig(A*A'); d=diag(D);
        V = [];         % matrix of eigenvectors of A'*A, or fxn handle,
                        % used (if present) when M>N
                        %   [V,D]=eig(A'*A); d=diag(D);
        d = [];         % vector of eigenvalues of A*A' or A'*A
        Uh = [];        % fxn handle for U', needed only if U is a fxn handle
        Vh = [];        % fxn handle for V', needed only if V is a fxn handle

        r1init = [];    % initial value of vector r1 (estimate of x)
        gam1xinit = 1e-8;  % initial value of scalar gam1x (precision on r1)
        p1init = [];    % initial value of vector p1 (estimate of z)
        gam1zinit = 1e-8;  % initial value of scalar gam1z (precision on p1)

        altUpdate = false; % alternately update x and z, or simultaneous?

        debug = true;   % debug switch: set =false for mild speedup 
        verbose = false;% verbose switch 
        silent = false; % disable all warnings 
        fxnErr1 = [];   % handle to a fxn of x1,z1 for error reporting, e.g.,
                        %   fxnErr1 = @(x1,z1) 10*log10( ...
                        %                    sum(abs(x1-xTrue).^2,1) ...
                        %                    ./sum(abs(xTrue).^2,1) );
        fxnErr2 = [];   % handle to another fxn of x1,z1 for error reporting
        fxnStop = [];   % handle to a stopping fxn of the form
                        %   fxnStop(i,err1(:,i),err2(:,i),...
                        %           r1old,r1,gam1x,x1,eta1x,...
                        %           p1old,p1,gam1z,z1,eta1z,...
                        %           r2old,r2,gam2x,x2,eta2x,...
                        %           p2old,p2,gam2z,z2,eta2z);
        histIntvl = 1;  % can save memory by decimating the saved history

        divChange = 1e-3; % amount to perturb input for denoiser Monte-Carlo
                        % divergence estimate; passed to FxnhandleEstimIn 
                        % when the denoiser reports a single output 

        r1old = [];     % used with warm-start and damping
        r2old = [];     % used with warm-start and damping
        p1old = [];     % used with warm-start and damping
        p2old = [];     % used with warm-start and damping
        x1old = [];     % used with warm-start and damping
        x2old = [];     % used with warm-start and damping
        z1old = [];     % used with warm-start and damping
        z2old = [];     % used with warm-start and damping
        gam1xold = [];  % used with warm-start and damping
        gam1zold = [];  % used with warm-start and damping
        gam2xold = [];  % used with warm-start and damping
        gam2zold = [];  % used with warm-start and damping
    end

    methods
        
        % Constructor with default options
        function opt = VampGlmOpt(varargin)
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
                        error('VampGlmOpt: %s is an unrecognized option', ...
                            num2str(varargin{i}));
                    end
                end
                return
            else
                error(['The VampGlmOpt constructor requires arguments ' ...
                    'be provided in pairs, e.g., VampGlmOpt(''verbose'',' ...
                    ' false, ''nitMax'', 50)'])
            end
        end % VampGlmOpt

        % Warm-start copy function
        function new = warmStart(old,estFin,varargin)
        %
        %   opt2 = opt1.warmStart(estFin[,'key1',value1,...])
        %
        % returns a copy of a VampOpt object with specified changes:
        % 1) warm-start values from the provided estFin struct
        % 2) any additional {key,value} pairs, e.g., 
        %       opt1 = VampOpt('verbose = false');
        %       estFin1 = VampGlmEst2(...,opt1);
        %       opt2 = opt1.warmStart(estFin1,'tol',1e-5,'nit',20);
        %       estFin2 = VampGlmEst2(...,opt2);

            new = old; % copy all options
            new.debug = false; % dont' debug during warm start

            % populate the warm start fields used by VampEstGlm2
            if ~isempty(estFin)
                new.r1init = estFin.r1; 
                new.p1init = estFin.p1; 
                new.gam1xinit = estFin.gam1x; 
                new.gam1zinit = estFin.gam1z; 
                new.r1old = estFin.r1old; 
                new.r2old = estFin.r2old; 
                new.p1old = estFin.p1old; 
                new.p2old = estFin.p2old; 
                new.x1old = estFin.x1old; 
                new.x2old = estFin.x2old; 
                new.z1old = estFin.z1old; 
                new.z2old = estFin.z2old; 
                new.gam1xold = estFin.gam1xold; 
                new.gam1zold = estFin.gam1zold; 
                new.gam2xold = estFin.gam2xold; 
                new.gam2zold = estFin.gam2zold; 
            end

            % override any additional fields
            for i = 1:2:length(varargin)
                new.(varargin{i}) = varargin{i+1};
            end

        end % warmStart

    end % methods

end
