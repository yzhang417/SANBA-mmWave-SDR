function [mse,rvar] = ampSE(estInAvg,wvar,delta,rvar0,mse0,nit)
% ampSE:  AMP state evolution ... a basic implementation
%
% DESCRIPTION:
% ------------
% AMP aims to iteratively recover a vector x from measurements of the form
%
%   y = A*x + sqrt(wvar)*randn(M,1) 
%
% where A is a known MxN matrix.
%
% When A is large and iid Gaussian with mean zero and variance 1/M,
% and independent of x, the per-iteration MSE of AMP can be
% characterized by a scalar state evolution:
%
%   MSE = estInAvg.mse(rvar) 
%   rvar = wvar + MSE/delta 
% 
% where rvar is the MSE on the denoiser input and delta=M/N. 
%
% SYNTAX:
% -------
% [mse,rvar] = ampSE(estInAvg, wvar, delta, rvar0, mse0, nit)
%
% INPUTS:
% -------
% estInAvg:  An averaged input estimator derived from the EstimInAvg class
% wvar:  Variance of measurement noise
% delta:  Sampling ratio M/N
% rvar0:  Initial value of rvar [optional]
% mse0:  Initial value of mse [optional]
% nit:  Number of iterations [optional]
%
% OUTPUTS:
% --------
% mse:  MSE versus iteration 
% rvar:  rvar versus iteration 
%
% Note: this routine assumes that norm(A,'fro')^2 = N.  If this is not
% the case, then wvar needs to be scaled appropriately.


% handle inputs
if nargin<3
  error('need at least 3 input arguments')
end

if (nargin<5)||isempty(mse0)
  mse0 = mean(abs(estInAvg.x).^2);
end

if (nargin<4)||isempty(rvar0)
  rvar0 = 100*mse0;
end

if (nargin<5)||isempty(nit)
  nit = 20;
end

% declare arrays
rvar = nan(1,nit);
mse = nan(1,nit);

% save initialization
rvar(1) = rvar0;
mse(1) = mse0;

% run state evolution
for i = 2:nit
  mse(i) = estInAvg.mse(rvar(i-1));
  rvar(i) = wvar + mse(i)/delta;
end
