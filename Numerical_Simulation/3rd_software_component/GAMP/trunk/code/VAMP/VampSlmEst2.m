function [x2,estFin,estHist] = VampSlmEst(denoiser,y,A,opt)
% VampSlmEst:  VAMP for the Standard Linear Model 
%
% DESCRIPTION
% -----------
% This implements VAMP for the 
% standard linear model, i.e., estimating a vector x from measurements 
%    y = A*x + w ,
% where w is additive white Gaussian noise and A is a known linear 
% operator.
%
% VAMP-SLM is an iterative algorithm where each iteration has two stages.
% The first stage performs MMSE (linear) estimation of x under the prior 
% p(x)=N(x;r1,1/gam1) and the likelihood p(y|x)=N(y;A*x,1/gamwHat), 
% where gamwHat is (an estimate of) the noise precision.  This stage 
% yields an estimate x1 with precision eta1.  The pair (x1,eta1) is then 
% converted into the pair (r2,gam2), where r2 is treated as an 
% AWGN-corrupted version of the true x with noise precision gam2.  
% The second stage of the iteration denoises r2 to produce the estimate 
% x2 with precision eta2.  Finally, the pair (x2,eta2) is converted to 
% the pair (r1,gam1), concluding the iteration.
% These iterations are repeated until r1 converges, and the final
% value of x2 reported as the estimate of x.
%
%
% USAGE
% -----
% [x2,estFin,estHist] = VampSlmEst(denoiser,y,A,[opt])
%
% denoiser: implements denoising of r2 to produce x2.  Several options:
%   1) handle to fxn of form [x2,xvar] = denoiser(r2,rvar), where r2 is 
%      the noisy signal and rvar is an estimate of the noise variance
%   2) handle to fxn of form x2 = denoiser(r2,rvar)
%   3) EstimIn object from GAMPmatlab
%
% y: an Mx1 vector of measurements (or MxL matrix in MMV case) 
%
% A: linear operator.  Several options:
%   1) matrix of size MxN
%   2) fxn handle, in which case opt.N must specify the # of columns in A
%   3) LinTrans object from GAMPmatlab
%
% opt: options structure, left empty or initialized via opt=VampSlmOpt.
%   See VampSlmOpt.m for details...
%
%
% FAST OPERATORS
% --------------
% VampSlmOpt calls "[U,D]=eig(A*A')" unless U and d=diag(D) are specified.
% Several options:
% 1) opt.U is an MxM matrix 
% 2) opt.U is a fxn handle for U and opt.Uh is a fxn handle for U' 
% In all cases, opt.d must be an Mx1 vector


% Process output format
if (nargout > 3),
    error('too many output arguments')
end
saveHist = (nargout >= 3);

% Get options
if (nargin < 3) || isempty(opt)
    opt = VampSlmOpt;   % default options
end
nitMax = opt.nitMax;  % maximum number of iterations 
tol = opt.tol;  % stopping tolerance
gamMin = opt.gamMin;  % minimum allowed precision
gamMax = opt.gamMax;  % maximum allowed precision
damp = opt.damp;  % damping parameter
dampGam = opt.dampGam;  % damping parameter on precisions
if isempty(dampGam), dampGam = damp; end;
learnNoisePrec = opt.learnNoisePrec; % learn noise precision? 
learnNoisePrecAlg = opt.learnNoisePrecAlg; % learning algorithm
learnNoisePrecNit = opt.learnNoisePrecNit; % nit for noise learning
learnNoisePrecTol = opt.learnNoisePrecTol; % tol for noise learning 
gamwHatMax = opt.NoisePrecMax; % max allowed noise precision 
verbose = opt.verbose; % verbose output?
silent = opt.silent; % disable all warnings?

% Process measurements y
if isa(y,'numeric')
    [M,L] = size(y);
    if (~silent) && all(y(:)==zeros(M*L,1))
      warning('y is zero valued!  May cause problems...') 
    end
else
    error('2nd input (y) must be a vector (or matrix) of measurements')
end

% Process linear transform A, length N, eigenvectors U, eigenvalues d
time_eig = nan;
if isa(A,'numeric')

    % treat the case where A is an explicit matrix
    [Mtest,N] = size(A);

    if Mtest~=M, 
        error('Number of rows in 2nd & 3rd inputs (y and A) do not match'); 
    end;

    % if needed, compute eigendecomposition of A*A'
    if isempty(opt.U)||isempty(opt.d)
        if verbose, 
            fprintf('Computing eigendecomposition of A*Ah...\n'); 
            pause(eps);
        end;
        tstart = tic;
          AAh = A*A';
          [U,D] = eig(0.5*(AAh+AAh'));
          d = diag(D);
          clear AAh D;
        time_eig = toc(tstart);
    else
        U = opt.U;
        d = opt.d;
        if length(d)~=M, error('Need length(d)==M'); end
    end

    % create function handles
    fxnA = @(x) A*x;
    fxnAh = @(z) A'*z;
    clear A;
    if isa(U,'function_handle')
      fxnU = U;
      if (~isempty(opt.Uh))&&isa(opt.Uh,'function_handle')
        fxnUh = opt.Uh;
      else
        error('Since opt.U is a fxn handle, opt.Uh must also be one')
      end
    else
      fxnU = @(x) U*x;
      fxnUh = @(x) U'*x;
      clear U;
    end

elseif isa(A,'function_handle')||any(strcmp(superclasses(A),'LinTrans'))

    % treat the case where A is a function handle or a LinTrans object
    if isa(A,'function_handle')

        fxnA = A;
        if isa(opt.Ah,'function_handle'), 
            fxnAh = opt.Ah;
        else
            error('opt.Ah must be a fxn handle to a linear operator');  
        end;
        if ~isempty(opt.N)
            if floor(opt.N)==opt.N, 
                N = opt.N;
            else
                error('opt.N must be an integer');  
            end
        else
            error('Since 3rd argument is a fxn handle, must specify opt.N');  
        end

    else % A is a LinTrans object

        fxnA = @(x) A.mult(x);
        fxnAh = @(x) A.multTr(x);
        [~,N] = size(A);

    end

    % test fxnA for size
    try
        z = fxnA(zeros(N,L));
        if any(size(z)~=[M,L]), 
            error('output size of 3rd input (A) doesnt match 2nd input (y)');
        end
        clear z;
    catch
        error('3rd input (A) doesnt accept inputs of correct size')
    end

    % test fxnAh for size
    try
        x = fxnAh(zeros(M,L));
        if any(size(x)~=[N,L]), 
            error('output size of opt.Ah doesnt match 2nd input (y)');
        end
        clear x;
    catch
        error('opt.Ah doesnt accept inputs of correct size')
    end

    % if needed, compute eigendecomposition of A*A'
    if isempty(opt.U)||isempty(opt.d)
        if verbose, 
            fprintf('Computing eigendecomposition of A*Ah...\n'); 
            pause(eps);
        end;
        tstart = tic;
          try
            AAh = fxnA(fxnAh(speye(M)));
          catch
            if verbose
                fprintf('fxnA doesnt support matrix argument!  Will slow down eigendecomposition...\n')
            end
            AAh = zeros(M);
            I = eye(M);
            for m=1:M, AAh(:,m) = fxnA(fxnAh(I(:,m))); end
          end
          [U,D] = eig(0.5*(AAh+AAh'));
          d = diag(D);
          clear AAh D;
        time_eig = toc(tstart);
    else
        U = opt.U;
        d = opt.d;
        if length(d)~=M, error('Need length(d)==M'); end
    end

    % create function handles for U and U'
    if isa(U,'function_handle')
      fxnU = U;
      if (~isempty(opt.Uh))&&isa(opt.Uh,'function_handle')
        fxnUh = opt.Uh;
      else
        error('Since opt.U is a fxn handle, opt.Uh must also be one')
      end
    else
      fxnU = @(x) U*x;
      fxnUh = @(x) U'*x;
    end
    clear U;

else
    error('3rd input must be a matrix or fxn handle to a linear operator')
end

% Load or create U'*A and A'*U handles
if isa(opt.UhA,'function_handle') 
  fxnUhA = opt.UhA;
else
  fxnUhA = @(x) fxnUh(fxnA(x));
end
if isa(opt.AhU,'function_handle') 
  fxnAhU = opt.AhU;
else
  fxnAhU = @(z) fxnAh(fxnU(z));
end

% Process denoiser
if isa(denoiser,'function_handle')
    if nargin(denoiser)~=2, error('need nargin(denoiser)==2'); end;
    % check if this function handle returns [xhat,xvar] or just [xhat]
    try % to produce two outputs
        [xhat,xvar] = denoiser(zeros(N,L),ones(1,L));
        fxnDenoise = @(rhat,rvar) denoiser(rhat,rvar); 
    catch % else turn into an EstimIn object
        divAvg = 1;
        denoiser1 = FxnhandleEstimIn(denoiser,...
                                     'changeFactor',opt.divChange,...
                                     'avg',divAvg); 
        fxnDenoise = @(rhat,rvar) denoiser1.estim(rhat,rvar); 
    end 
elseif any(strcmp(superclasses(denoiser),'EstimIn')) 
    % turns EstimIn object into a function handle
    fxnDenoise = @(rhat,rvar) denoiser.estim(rhat,rvar); 
else
    error('First input (denoiser) must be either an EstimIn object, a fxn handle that accepts [rhat,rvar] and produces [xhat,xvar], or a fxn handle that accepts [rhat,rvar] and produces only xhat.')
end 

% Process error-reporting function
if isa(opt.fxnErr,'function_handle')
    fxnErr = opt.fxnErr;
else
    fxnErr = [];
end

% Process general reporting function
if isa(opt.fxnGen,'function_handle')
    fxnGen = opt.fxnGen;
else
    fxnGen = [];
end

% Process stop function
if isa(opt.fxnStop,'function_handle')
    fxnStop = opt.fxnStop;
else
    fxnStop = [];
end

% Initialize noise precision & learning
gamwHatMin = M*L/norm(y,'fro')^2;
if ~isempty(opt.NoisePrecInit)
    gamwHat = opt.NoisePrecInit;
    if gamwHat>gamwHatMax
        warning('opt.NoisePrecInit too large')
        gamwHat = gamwHatMax;
    elseif gamwHat<gamwHatMin
        warning('opt.NoisePrecInit too small')
        gamwHat = gamwHatMin;
    end
else
    gamwHat = gamwHatMin;
end

% Compute some constants
d = max(d,eps); % avoid zero-valued eigenvalues
s = sqrt(d); % singular values of A
Uhy = fxnUh(y); % U'*y where [U,D]=eig(A*A')

% Prepare for saving history 
if saveHist
    histIntvl = opt.histIntvl;
    nitSave = floor(nitMax/histIntvl);
    estHist.r1 = nan(N,L,nitSave);
    estHist.gam1 = nan(L,nitSave);
    estHist.x1 = nan(N,L,nitSave);
    estHist.gamwHat = nan(1,nitSave);
    estHist.eta1 = nan(L,nitSave);
    estHist.r2 = nan(N,L,nitSave);
    estHist.gam2 = nan(L,nitSave);
    estHist.x2 = nan(N,L,nitSave); 
    estHist.xvar2 = nan(N,L,nitSave);
    estHist.eta2 = nan(L,nitSave);
    estHist.err = nan(L,nitSave);
    estHist.gen = nan(L,nitSave);
end

% Initialize VAMP
if ~isempty(opt.r1init)
    r1 = opt.r1init;    
    if size(r1,2)==1, r1 = r1*ones(1,L); end;
else
    r1 = zeros(N,L); % default initialization
end
if ~isempty(opt.gam1init)
    gam1 = opt.gam1init;    
else
%   gam1 = 0.1*sum(d.^2)*L/norm(fxnAh(y),'fro').^2;  % default, ~= 0.1/xvar0
    gam1 = 0.1*sum(d.^2)*L/norm(bsxfun(@times,s,Uhy),'fro').^2;  % default, ~= 0.1/xvar0
end
if size(gam1,2)==1, gam1 = gam1*ones(1,L); end;
if ~isempty(opt.r2init)
    r2 = opt.r2init;    
    if size(r2,2)==1, r2 = r2*ones(1,L); end;
else
    r2 = nan(N,L); % default initialization
end
if ~isempty(opt.gam2init)
    gam2 = opt.gam2init;    
else
    gam2 = nan(1,L); % default initialization
end
if size(gam2,2)==1, gam2 = gam2*ones(1,L); end;
r1old = nan(N,L); % used by fxnGen
r2old = nan(N,L); % used by fxnGen
i = 0;
stop = false;

% Run VAMP
err = nan(L,nitMax);
gen = nan(L,nitMax);
gamwHat_ = nan(1,nitMax);
nitNoise = nan(1,nitMax);
time = nan(1,nitMax);
tstart = tic;
while ~stop

  %-----update counter
  i = i + 1;

  %-----damping 
  dampG1 = 1;
  dampR1 = 1;
  if i>1 
    gam1 = dampG1*gam1 + (1-dampG1)*gam1old;
    r1 = dampR1*r1 + (1-dampR1)*r1old;
  end;

  %-----first half of iteration
  UyAr1 = Uhy-fxnUhA(r1);
  gam1overD = (1./d)*gam1;

  if learnNoisePrec
     gam1overD_UyAr1_Sq = (gam1overD.^2).*abs(UyAr1).^2;
     switch learnNoisePrecAlg
       case 'EM'
         for ii=1:learnNoisePrecNit 
           gamwHat_old = gamwHat;
           % note that resNormSq = sum(abs(y-A*x1).^2,1);
           resNormSq = sum(gam1overD_UyAr1_Sq./((gamwHat+gam1overD).^2)); 
           gamwHat = M/mean(resNormSq + sum(1./(gamwHat+gam1overD),1) );
           gamwHat = min(gamwHat,gamwHatMax);
           gamwHat = max(gamwHat,gamwHatMin);
           if abs(gamwHat_old-gamwHat)/gamwHat < learnNoisePrecTol, break; end;
         end
         nitNoise(i)=ii;
       case 'Newton'
         hessReg = 1e-15; % Hessian regularization
         gamwHatMin = 1e-1; % minimum allowed value of gamwHat
         for ii=1:learnNoisePrecNit 
           grad = -M/gamwHat +mean(sum(...
                  gam1overD_UyAr1_Sq./((gamwHat+gam1overD).^2) ...
                  + 1./(gamwHat+gam1overD) ));
           hess = M/gamwHat^2 -mean(sum(...
                  2*gam1overD_UyAr1_Sq./((gamwHat+gam1overD).^3) ...
                  + 1./((gamwHat+gam1overD).^2) ));
           gamwHat_old = gamwHat;
           gamwHat = gamwHat - grad/(abs(hess)+hessReg);
           gamwHat = min(max(gamwHat,gamwHatMin),gamwHatMax);
           if abs(gamwHat_old-gamwHat)/gamwHat < learnNoisePrecTol, break; end;
         end
         nitNoise(i)=ii;
       otherwise
         error('unknown type of learnNoisePrecAlg')
     end
     if (~silent)&&(gamwHat<gamMin), 
         warning('gamwHat=%g too small, iter=%i',gamwHat,i); 
     end 
  end % learnNoisePrec
  gamwHat_(i) = gamwHat;

  oneMinusAlf1 = (gamwHat/N)*sum(1./(gamwHat+gam1overD),1); 
  x1 = r1 + fxnAhU(bsxfun(@rdivide,UyAr1,bsxfun(@plus,d,(1/gamwHat)*gam1)));

  gam2old = gam2;
  gam2 = oneMinusAlf1./(1-oneMinusAlf1).*gam1;
  eta1 = gam1 + gam2; % reported but not used
  if any(gam2<gamMin), 
    bad = find(gam2<gamMin);
    if ~silent
      warning('gam2=%g too small at column %d, iter=%d\n',[gam2(bad);bad;i*ones(size(bad))]); 
    end
    gam2(bad) = gam2old(bad);
  end 
  if any(gam2>gamMax), 
    bad = find(gam2>gamMax);
    if ~silent
      warning('gam2=%g too large at column %d, iter=%d\n',[gam2(bad);bad;i*ones(size(bad))]); 
    end
    gam2(bad) = gam2old(bad);
  end 
  gam2 = min(max(gam2,gamMin),gamMax);

  r2old = r2;
  r2 = r1 + bsxfun(@times,x1-r1,1./oneMinusAlf1);

  %-----damping 
  if damp>0
    dampR2 = damp;
    dampG2 = dampGam;
  else % specifying damp<=0 triggers adaptive damping
    dampR2 = 2./(1+max([gam2./gam1;gam1./gam2],[],1));
    dampG2 = min(dampR2,dampGam);
  end
  if i>1 % apply damping 
    gam2 = dampG2.*gam2 + (1-dampG2).*gam2old;
    r2 = bsxfun(@times,dampR2,r2) + bsxfun(@times,(1-dampR2),r2old);
  end;

  %-----second half of iteration
  [x2,xvar2] = fxnDenoise(r2,ones(N,1)*(1./gam2));
  eta2 = 1./mean(xvar2+eps,1); % 1xL row vector in MMV case

  %-----record progress
  if ~isempty(fxnErr)
    err(:,i) = fxnErr(x2).'; % evaluate error function
  end
  if ~isempty(fxnGen)
    gen(:,i) = fxnGen(r1old,r1,gam1,x1,eta1,...
                      r2old,r2,gam2,x2,eta2); % evaluate general function
  end

    %-----save history
    if saveHist && rem(i,histIntvl)==0
      iHist = i/histIntvl;
      estHist.r1(:,:,iHist) = r1;
      estHist.gam1(:,iHist) = gam1;
      estHist.x1(:,:,iHist) = x1;
      estHist.gamwHat(iHist) = gamwHat;
      estHist.eta1(:,iHist) = eta1;
      estHist.r2(:,:,iHist) = r2;
      estHist.gam2(:,iHist) = gam2;
      estHist.x2(:,:,iHist) = x2;
      estHist.xvar2(:,:,iHist) = xvar2;
      estHist.eta2(:,iHist) = eta2.';
      estHist.err(:,iHist) = err(:,i).';
      estHist.gen(:,iHist) = gen(:,i).';
    end

  %-----second half of iteration (continued)
  gam1old = gam1;
  gam1 = eta2 - gam2; % compute L copies 
  if any(gam1<gamMin), 
    bad = find(gam1<gamMin);
    if ~silent
      warning('gam1=%g too small at column %d, iter=%d\n',[gam1(bad);bad;i*ones(size(bad))]); 
    end
    gam1(bad) = gam1old(bad);
    %fprintf('corrected small gam1: %g\n',gam1(bad))
  end 
  if any(gam1>gamMax), 
    bad = find(gam1>gamMax);
    if ~silent
      warning('gam1=%g too large at column %d, iter=%d\n',[gam1(bad);bad;i*ones(size(bad))]); 
    end
    gam1(bad) = gam1old(bad);
    %fprintf('corrected large gam1: %g\n',gam1(bad))
  end 
  gam1 = min(max(gam1,gamMin),gamMax);

  r1old = r1;
  r1 = bsxfun(@rdivide, bsxfun(@times,x2,eta2)-bsxfun(@times,r2,gam2), gam1);

  %-----check for NaN
  if (~silent)
    if any(isnan(r1(:))), 
      bad = find(sum(isnan(r1),1));
      warning('r1=NaN at column %d, iter=%i\n',[bad;i*ones(size(bad))]);
    end
    if any(isnan(x1(:))), 
      bad = find(sum(isnan(x1),1));
      warning('x1=NaN at column %d, iter=%i\n',[bad;i*ones(size(bad))]);
    end;
    if any(isnan(r2(:))), 
      bad = find(sum(isnan(r2),1));
      warning('r2=NaN at column %d, iter=%i\n',[bad;i*ones(size(bad))]);
    end;
    if any(isnan(x2(:))), 
      bad = find(sum(isnan(x2),1));
      warning('x2=NaN at column %d, iter=%i\n',[bad;i*ones(size(bad))]);
    end;
  end;

    %-----report progress
    if verbose
      if ~isempty(fxnErr)
          fprintf('i=%3i: eta2/eta1=%8.5f, |dr|/|r|=%9.3e, err=%8.5g\n',...
                  [i*ones(L,1),...
                   (eta2./eta1).',...
                   (sqrt(sum(abs(r1-r1old).^2,1)./sum(abs(r1).^2,1))).',...
                   err(:,i)].')
      else
          fprintf('i=%3i: eta2/eta1=%8.5f, |dr|/|r|=%9.3e\n',...
                  [i*ones(L,1),...
                   (eta2./eta1).',...
                   (sqrt(sum(abs(r1-r1old).^2)./sum(abs(r1).^2))).'].')
      end
    end

  %-----stopping rule
  if ~isempty(fxnStop)
    stop = fxnStop(i,err(:,i),r1old,r1,gam1,x1,eta1,...
                              r2old,r2,gam2,x2,eta2);
  end
  if all(sum(abs(r1-r1old).^2)./sum(abs(r1).^2) < tol^2)||(i>=nitMax)
    stop = true;
  end
  if stop 
      err = err(:,1:i); % trim
      if learnNoisePrec
          nitNoise = nitNoise(1:i); % trim
      else
          nitNoise = nan; 
      end
  end

  %-----measure time
  time(i) = toc(tstart);
end % while
time_iters = toc(tstart);

% Export Outputs
estFin.r1old = r1old;
estFin.r1 = r1;
estFin.gam1 = gam1;
estFin.x1 = x1;
estFin.gamwHat = gamwHat_(1:i);
estFin.nitNoise = nitNoise;
estFin.eta1 = eta1;
estFin.r2 = r2;
estFin.gam2 = gam2;
estFin.x2 = x2; % this is the main output: the final estimate of x
estFin.xvar2 = xvar2;
estFin.eta2 = eta2;
estFin.nit = i;
estFin.err = err;
estFin.gen = gen(1:i);
estFin.time_eig = time_eig;
estFin.time = time(1:i);

% Trim history
if saveHist
    iTrim = 1:floor(i/histIntvl);
    estHist.r1 = squeeze(estHist.r1(:,:,iTrim));
    estHist.gam1 = estHist.gam1(:,iTrim);
    estHist.x1 = squeeze(estHist.x1(:,:,iTrim));
    estHist.gamwHat = estHist.gamwHat(:,iTrim);
    estHist.eta1 = estHist.eta1(:,iTrim);
    estHist.r2 = squeeze(estHist.r2(:,:,iTrim));
    estHist.gam2 = estHist.gam2(:,iTrim);
    estHist.x2 = squeeze(estHist.x2(:,:,iTrim)); 
    estHist.xvar2 = squeeze(estHist.xvar2(:,:,iTrim));
    estHist.eta2 = estHist.eta2(:,iTrim);
    estHist.err = estHist.err(:,iTrim);
end

