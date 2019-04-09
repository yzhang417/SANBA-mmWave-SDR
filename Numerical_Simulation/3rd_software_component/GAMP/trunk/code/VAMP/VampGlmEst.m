function [x1,estFin,estHist] = VampGlmEst(denoiser,likelihood,A,opt)
% VampGlmEst:  VAMP for the Generalized Linear Model 
%
% DESCRIPTION
%------------
% This implements VAMP for the generalized linear model, i.e., estimating 
% a vector x from measurements 
%    y ~ p(y|A*x)
% where A is a known linear operator and p(y|z) is a known conditional pdf.
%
% For a detailed description, see https://arxiv.org/abs/1612.01186
%
%
% USAGE
% -----
% [x1,estFin,estHist] = VampGlmEst(denoiser,likelihood,A,[opt])
%
% denoiser: implements denoising of "r=x+AWGN(0,rvar)" to produce xhat.  
%   Several options:
%   1) handle to fxn of form [xhat,xvar] = denoiser(r,rvar) 
%   2) handle to fxn of form xhat = denoiser(r,rvar)
%   3) EstimIn object from GAMPmatlab
%
% likelihood: implements estimation of z with prior z~AWGN(p,pvar) from
%   measurements y to produce zhat.
%   1) handle to fxn of form [zhat,zvar] = likelihood(p,pvar)  
%   2) handle to fxn of form zhat = likelihood(p,pvar)
%   3) EstimOut object from GAMPmatlab
%
% A: linear operator.  Several options:
%   1) matrix of size MxN
%   2) fxn handle, in which case opt.N must specify the # of columns in A
%   3) LinTrans object from GAMPmatlab
%
% opt: options structure, left empty or initialized via opt=VampGlmOpt.
%   See VampGlmOpt.m for details...
%
%
% FAST OPERATORS
% --------------
% If M<N:
%   VampGlmOpt calls "[U,D]=eig(A*A')" unless U and d=diag(D) are specified.
%   Several options:
%   1) opt.U is an MxM matrix 
%   2) opt.U is a fxn handle for U and opt.Uh is a fxn handle for U' 
%   In all cases, opt.d must be an Mx1 vector
% If M>=N:
%   VampGlmOpt calls "[V,D]=eig(A'*A)" unless V and d=diag(D) are specified.
%   Several options:
%   1) opt.V is an NxN matrix 
%   2) opt.V is a fxn handle for V and opt.Vh is a fxn handle for V' 
%   In all cases, opt.d must be an Nx1 vector

     
% Process output format
if (nargout > 3),
    error('too many output arguments')
end
saveHist = (nargout >= 3);

% Get options
if (nargin < 3) || isempty(opt)
    opt = VampGlmOpt;   % default options
end
nitMax = opt.nitMax;  % max number of iterations 
tol = opt.tol;  % stopping tolerance
gamMin = opt.gamMin;  % minimum allowed precision
gamMax = opt.gamMax;  % maximum allowed precision
damp = opt.damp;  % damping parameter
dampGam = opt.dampGam; % damping parameter for precisions
if isempty(dampGam), dampGam = damp; end;
verbose = opt.verbose; % verbose output?
silent = opt.silent; % disable all warnings?
debug = opt.debug; % disable debugging?

% Process likelihood
if isa(likelihood,'function_handle')
    error('fxn handles not yet implemented for likelihood')
    % would need to add options opt.M & opt.L ...
elseif any(strcmp(superclasses(likelihood),'EstimOut')) 
    [M,L] = size(likelihood);
    if isempty(M)||isempty(L), error('size(likelihood) must return [M,L]'); end;
    % turn EstimOut object into a function handle
    fxnLike = @(phat,pvar) likelihood.estim(phat,pvar); 
else
    error('Second input (likelihood) must be either an EstimOut object, a fxn handle that accepts [phat,pvar] and produces [zhat,zvar], or a fxn handle that accepts [rhat,rvar] and produces only xhat.')
end 

% Process linear transform A, length N, eigenvectors U, eigenvalues d
time_eig = nan;
if isa(A,'numeric')

    % treat the case where A is an explicit matrix
    [Mtest,N] = size(A);

    if Mtest~=M, 
        error('Number of rows in 2nd & 3rd inputs (y and A) do not match'); 
    end;

    % if needed, compute eigendecomposition of A*A' or A'*A
    if M<N
        if isempty(opt.U)||isempty(opt.d)
            if verbose, 
                fprintf('Computing eigendecomposition of A*A''...\n'); 
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
        end
    else % M>=N
        if isempty(opt.V)||isempty(opt.d)
            if verbose, 
                fprintf('Computing eigendecomposition of A''*A...\n'); 
                pause(eps);
            end;
            tstart = tic;
              AhA = A'*A;
              [V,D] = eig(0.5*(AhA+AhA'));
              d = diag(D);
              clear AhA D;
            time_eig = toc(tstart);
        else
            V = opt.V;
            d = opt.d;
        end
    end

    % create function handles
    fxnA = @(x) A*x;
    fxnAh = @(z) A'*z;
    clear A;
    if M<N
      if isa(U,'function_handle')
        fxnU = U;
        if (~isempty(opt.Uh))&&isa(opt.Uh,'function_handle')
          fxnUh = opt.Uh;
        else
          error('Since opt.U is a fxn handle, opt.Uh must also be one')
        end
      else
        fxnU = @(z) U*z;
        fxnUh = @(z) U'*z;
        clear U;
      end
    else %M>=N
      if isa(V,'function_handle')
        fxnV = V;
        if (~isempty(opt.Vh))&&isa(opt.Vh,'function_handle')
          fxnVh = opt.Vh;
        else
          error('Since opt.V is a fxn handle, opt.Vh must also be one')
        end
      else
        fxnV = @(x) V*x;
        fxnVh = @(x) V'*x;
        clear V;
      end
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
    if debug
      try
        z = fxnA(zeros(N,L));
        if any(size(z)~=[M,L]), 
            error('output size of 3rd input (A) doesnt match 2nd input (likelihood)');
        end
        clear z;
      catch
        error('3rd input (A) doesnt accept inputs of correct size')
      end
    end

    % test fxnAh for size
    if debug
      try
        x = fxnAh(zeros(M,L));
        if any(size(x)~=[N,L]), 
            error('output size of opt.Ah doesnt match 2nd input (likelihood)');
        end
        clear x;
      catch
        error('opt.Ah doesnt accept inputs of correct size')
      end
    end

    % if needed, compute eigendecomposition of A*A' or A'*A
    if M<N
        if isempty(opt.U)||isempty(opt.d)
            if verbose, 
                fprintf('Computing eigendecomposition of A*A''...\n'); 
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
        end
    else % M>=N
        if isempty(opt.V)||isempty(opt.d)
            if verbose, 
                fprintf('Computing eigendecomposition of A''*A...\n'); 
                pause(eps);
            end;
            tstart = tic;
              try
                AhA = fxnAh(fxnA(speye(N)));
              catch
                if verbose
                    fprintf('fxnA doesnt support matrix argument!  Will slow down eigendecomposition...\n')
                end
                AhA = zeros(N);
                I = eye(N);
                for n=1:N, AhA(:,n) = fxnAh(fxnA(I(:,n))); end
              end
              [V,D] = eig(0.5*(AhA+AhA'));
              d = diag(D);
              clear AhA D;
            time_eig = toc(tstart);
        else
            V = opt.V;
            d = opt.d;
        end
    end

    % create function handles for U and U'
    if M<N
        if isa(U,'function_handle')
          fxnU = U;
          if (~isempty(opt.Uh))&&isa(opt.Uh,'function_handle')
            fxnUh = opt.Uh;
          else
            error('Since opt.U is a fxn handle, opt.Uh must also be one')
          end
        else
          fxnU = @(z) U*z;
          fxnUh = @(z) U'*z;
          clear U;
        end
    else % M>=N
        if isa(V,'function_handle')
          fxnV = V;
          if (~isempty(opt.Vh))&&isa(opt.Vh,'function_handle')
            fxnVh = opt.Vh;
          else
            error('Since opt.V is a fxn handle, opt.Vh must also be one')
          end
        else
          fxnV = @(x) V*x;
          fxnVh = @(x) V'*x;
          clear V;
        end
    end

else
    error('3rd input must be a matrix or fxn handle to a linear operator')
end
del = M/N;

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

% Process error-reporting functions
if isa(opt.fxnErr1,'function_handle')
    fxnErr1 = opt.fxnErr1;
else
    fxnErr1 = [];
end
if isa(opt.fxnErr2,'function_handle')
    fxnErr2 = opt.fxnErr2;
else
    fxnErr2 = [];
end

% Process stop function
if isa(opt.fxnStop,'function_handle')
    fxnStop = opt.fxnStop;
else
    fxnStop = [];
end

% Prepare for saving history 
if saveHist
    histIntvl = opt.histIntvl;
    nitSave = floor(nitMax/histIntvl);
    estHist.r1 = nan(N,L,nitSave);
    estHist.gam1x = nan(L,nitSave);
    estHist.p1 = nan(M,L,nitSave);
    estHist.gam1z = nan(L,nitSave);
    estHist.x1 = nan(N,L,nitSave);
    estHist.eta1x = nan(L,nitSave);
    estHist.z1 = nan(M,L,nitSave);
    estHist.eta1z = nan(L,nitSave);
    estHist.r2 = nan(N,L,nitSave);
    estHist.gam2x = nan(L,nitSave);
    estHist.p2 = nan(M,L,nitSave);
    estHist.gam2z = nan(L,nitSave);
    estHist.x2 = nan(N,L,nitSave); 
    estHist.eta2x = nan(L,nitSave);
    estHist.z2 = nan(M,L,nitSave); 
    estHist.eta2z = nan(L,nitSave);
    estHist.alf = nan(L,nitSave);
    estHist.err1 = nan(L,nitSave);
    estHist.err2 = nan(L,nitSave);
end

% Initialize VAMP
if ~isempty(opt.r1init)
    r1 = opt.r1init;
    if size(r1,2)==1, r1 = r1*ones(1,L); end;
else
    r1 = zeros(N,L); % default initialization
end
if ~isempty(opt.p1init)
    p1 = opt.p1init;    
    if size(p1,2)==1, p1 = p1*ones(1,L); end;
else
    p1 = zeros(M,L); % default initialization
end
gam1x = opt.gam1xinit;    
gam1z = opt.gam1zinit;
if size(gam1x,2)==1, gam1x = gam1x*ones(1,L); end;
if size(gam1z,2)==1, gam1z = gam1z*ones(1,L); end;
i = 0;
stop = false;

% Run VAMP
err1 = nan(L,nitMax);
err2 = nan(L,nitMax);
tstart = tic;
while ~stop

  %-----update counter
  i = i + 1;

  r1old = r1;
  p1old = p1;
  gam1xold = gam1x;
  gam1zold = gam1z;
  if i>1,       
      r2old = r2;
      p2old = p2;
      x1old = x1; 
      z2old = z2;
      gam2xold = gam2x;
      gam2zold = gam2z;
  end

  %-----first half of iteration
  [z1,zvar1] = fxnLike(p1,ones(M,1)*(1./gam1z));
  eta1z = 1./mean(zvar1,1); 
  gam2z = eta1z - gam1z; 
  p2 = bsxfun(@rdivide,bsxfun(@times,z1,eta1z)-bsxfun(@times,p1,gam1z),gam2z); 
  if ~silent
    if any(gam2z<gamMin), 
      bad = find(gam2z<gamMin);
      warning('gam2z=%g too small at column %d, iter=%d\n',...
              [gam2z(bad);bad;i*ones(size(bad))]);
      %gam2z(bad) = gam2zold(bad);
    end
  end
  if ~silent
    if any(gam2z>gamMax), 
      bad = find(gam2z>gamMax);
      warning('gam2z=%g too large at column %d, iter=%d\n',...
              [gam2z(bad);bad;i*ones(size(bad))]);
      %gam2z(bad) = gam2zold(bad);
    end
  end
  gam2z = min(max(gam2z,gamMin),gamMax);
  if i>1,       % apply damping
      gam2z = dampGam*gam2z + (1-dampGam)*gam2zold;
  end
  % NOTE: IN AWGN CASE gam2z=gamwHat AND p2=y

  [x1,xvar1] = fxnDenoise(r1,ones(N,1)*(1./gam1x));
  eta1x = 1./mean(xvar1,1); 
  if i>1,       % apply damping
      x1 = damp*x1 + (1-damp)*x1old;
  end
  gam2x = eta1x - gam1x; 
  r2 = bsxfun(@rdivide,bsxfun(@times,x1,eta1x)-bsxfun(@times,r1,gam1x),gam2x);
  if ~silent
    if any(gam2x<gamMin), 
      bad = find(gam2x<gamMin);
      warning('gam2x=%g too small at column %d, iter=%d\n',...
              [gam2x(bad);bad;i*ones(size(bad))]);
      %gam2x(bad) = gam2xold(bad);
    end
  end
  if ~silent
    if any(gam2x>gamMax), 
      bad = find(gam2x>gamMax);
      warning('gam2x=%g too large at column %d, iter=%d\n',...
              [gam2x(bad);bad;i*ones(size(bad))]);
      %gam2x(bad) = gam2xold(bad);
    end
  end
  gam2x = min(max(gam2x,gamMin),gamMax);

  %-----second half of iteration
  gam2Rat = gam2x./gam2z;
  inv_d_gam2Rat = 1./bsxfun(@plus,d,gam2Rat);
  alf = 1-(1/N)*(d'*inv_d_gam2Rat) + eps; % eps protects against alf==1
  beta = (1/del)*(1-alf);
  if M<N
    Ar2 = fxnA(r2);
    Up2Ar2_scaled = fxnUh(p2-Ar2).*inv_d_gam2Rat;
    x2 = r2 + fxnAh(fxnU(Up2Ar2_scaled));
    z2 = Ar2 + fxnU(bsxfun(@times,d,Up2Ar2_scaled));
  else
    Vr2Ap2 = fxnVh(bsxfun(@times,r2,gam2Rat)+fxnAh(p2));
    x2 = fxnV(Vr2Ap2.*inv_d_gam2Rat);
    z2 = fxnA(x2);
  end
  eta2x = gam2x./alf; % reported but not used
  eta2z = gam2z./beta; % reported but not used
  if i>1,       % apply damping
      z2 = damp*z2 + (1-damp)*z2old;
  end

  %-----record/report progress
  if ~isempty(fxnErr1)
      err1(:,i) = fxnErr1(x1,z1).'; % evaluate 1st error function
  end
  if ~isempty(fxnErr2)
      err2(:,i) = fxnErr2(x1,z1).'; % evaluate 2nd error function
  end

  if verbose
    if isempty(fxnErr1)&&isempty(fxnErr2)
      fprintf('i=%3i: eta2x/eta1x=%8.5f, eta2z/eta1z=%8.5f\n',...
              [i*ones(L,1),(eta2x./eta1x).',(eta2z./eta1z).'].')
    elseif isempty(fxnErr2)
      fprintf('i=%3i: eta2x/eta1x=%8.5f, eta2z/eta1z=%8.5f, err1=%8.5g\n',...
              [i*ones(L,1),(eta2x./eta1x).',(eta2z./eta1z).',err1(:,i)].')
    elseif isempty(fxnErr1)
      fprintf('i=%3i: eta2x/eta1x=%8.5f, eta2z/eta1z=%8.5f, err2=%8.5g\n',...
              [i*ones(L,1),(eta2x./eta1x).',(eta2z./eta1z).',err2(:,i)].')
    else
      fprintf('i=%3i: eta2x/eta1x=%8.5f, eta2z/eta1z=%8.5f, err1=%8.5g, err2=%8.5g\n',...
              [i*ones(L,1),(eta2x./eta1x).',(eta2z./eta1z).',err1(:,i),err2(:,i)].')
    end
  end

  %-----save history
  if saveHist && rem(i,histIntvl)==0
      iHist = i/histIntvl;
      estHist.r1(:,:,iHist) = r1;
      estHist.gam1x(:,iHist) = gam1x;
      estHist.p1(:,:,iHist) = p1;
      estHist.gam1z(:,iHist) = gam1z;
      estHist.x1(:,:,iHist) = x1;
      estHist.eta1x(:,iHist) = eta1x;
      estHist.z1(:,:,iHist) = z1;
      estHist.eta1z(:,iHist) = eta1z;
      estHist.r2(:,:,iHist) = r2;
      estHist.gam2x(:,iHist) = gam2x;
      estHist.p2(:,:,iHist) = p2;
      estHist.gam2z(:,iHist) = gam2z;
      estHist.x2(:,:,iHist) = x2;
      estHist.eta2x(:,iHist) = eta2x.';
      estHist.z2(:,:,iHist) = z2;
      estHist.eta2z(:,iHist) = eta2z.';
      estHist.alf(:,iHist) = alf.';
      estHist.err1(:,iHist) = err1(:,i).';
      estHist.err2(:,iHist) = err2(:,i).';
  end

  %-----continue second half of iteration
  r1 = bsxfun(@times,x2-bsxfun(@times,r2,alf),1./(1-alf));
  p1 = bsxfun(@times,z2-bsxfun(@times,p2,beta),1./(1-beta));
  gam1x = gam2x.*(1-alf)./alf;
  if ~silent
    if any(gam1x<gamMin),
      bad = find(gam1x<gamMin);
      warning('gam1x=%g too small at column %d, iter=%d\n',...
              [gam1x(bad);bad;i*ones(size(bad))]);
      %gam1x(bad) = gam1xold(bad);
    end
  end
  if ~silent
    if any(gam1x>gamMax), 
      bad = find(gam1x>gamMax);
      warning('gam1x=%g too large at column %d, iter=%d\n',...
              [gam1x(bad);bad;i*ones(size(bad))]);
      %gam1x(bad) = gam1xold(bad);
    end
  end
  gam1x = min(max(gam1x,gamMin),gamMax);
  gam1z = gam2z.*(1-beta)./beta;
  if ~silent
    if any(gam1z<gamMin), 
      bad = find(gam1z<gamMin);
      warning('gam1z=%g too small at column %d, iter=%d\n',...
              [gam1z(bad);bad;i*ones(size(bad))]);
      %gam1z(bad) = gam1zold(bad);
    end
  end
  if ~silent
    if any(gam1z>gamMax), 
      bad = find(gam1z>gamMax);
      warning('gam1z=%g too large at column %d, iter=%d\n',...
              [gam1z(bad);bad;i*ones(size(bad))]);
      %gam1z(bad) = gam1zold(bad);
    end
  end
  gam1z = min(max(gam1z,gamMin),gamMax);
  if i>1,       % apply damping
      gam1x = dampGam*gam1x + (1-dampGam)*gam1xold;
  end
  
  %-----stopping rule
  if (~isempty(fxnStop))&&(i>1)
    stop = fxnStop(i,err1(:,i),err2(:,i),...
                   r1old,r1,gam1x,x1,eta1x,...
                   p1old,p1,gam1z,z1,eta1z,...
                   r2old,r2,gam2x,x2,eta2x,...
                   p2old,p2,gam2z,z2,eta2z);
  end
  if all(sum(abs([r1;p1]-[r1old;p1old]).^2)./sum(abs([r1;p1]).^2) ...
         < tol^2)||(i>=nitMax)
      stop = true;
  end
  if stop 
      err1 = err1(:,1:i); % trim
      err2 = err2(:,1:i); % trim
  end
end % while
time_iters = toc(tstart);

% Export Outputs
estFin.r1old = r1old;
estFin.p1old = p1old;
estFin.r1 = r1;
estFin.gam1x = gam1x;
estFin.p1 = p1;
estFin.gam1z = gam1z;
estFin.x1 = x1; % this is the main output: the final estimate of x
estFin.xvar1 = xvar1;
estFin.eta1x = eta1x;
estFin.z1 = z1;
estFin.zvar1 = zvar1;
estFin.eta1z = eta1z;
estFin.r2 = r2;
estFin.gam2x = gam2x;
estFin.p2 = p2;
estFin.gam2z = gam2z;
estFin.x2 = x2; 
estFin.z2 = z2; 
estFin.eta2x = eta2x;
estFin.eta2z = eta2z;
estFin.nit = i;
estFin.err1 = err1;
estFin.err2 = err2;
estFin.time_eig = time_eig;
estFin.time_iters = time_iters;

% Trim history
if saveHist
    iTrim = 1:floor(i/histIntvl);
    estHist.r1 = estHist.r1(:,:,iTrim);
    estHist.gam1x = estHist.gam1x(:,iTrim);
    estHist.p1 = estHist.p1(:,:,iTrim);
    estHist.gam1z = estHist.gam1z(:,iTrim);
    estHist.x1 = estHist.x1(:,:,iTrim);
    estHist.eta1x = estHist.eta1x(:,iTrim);
    estHist.z1 = estHist.z1(:,:,iTrim);
    estHist.eta1z = estHist.eta1z(:,iTrim);
    estHist.r2 = estHist.r2(:,:,iTrim);
    estHist.gam2x = estHist.gam2x(:,iTrim);
    estHist.p2 = estHist.p2(:,:,iTrim);
    estHist.gam2z = estHist.gam2z(:,iTrim);
    estHist.x2 = estHist.x2(:,:,iTrim); 
    estHist.eta2x = estHist.eta2x(:,iTrim);
    estHist.z2 = estHist.z2(:,:,iTrim); 
    estHist.eta2z = estHist.eta2z(:,iTrim);
    estHist.alf = estHist.alf(:,iTrim);
    estHist.err1 = estHist.err1(:,iTrim);
    estHist.err2 = estHist.err2(:,iTrim);
end

