function mse = VampGlmSE(estInAvg,estOutAvg,d,N,del,opt)

if nargin<6
  nitMax = 20;
  gam1xinit = 1e-6;
  gam1zinit = 1e-6;
  altUpdate = false;
  damp = 1;
  dampGam = 1;
else
  nitMax = opt.nitMax;
  gam1xinit = opt.gam1xinit;
  gam1zinit = opt.gam1zinit;
  altUpdate = opt.altUpdate;
  damp = opt.damp;
  dampGam = opt.dampGam;
end

% run state-evolution
gam1x = gam1xinit; 
gam1z = gam1zinit; 
gam2x = gam1x; % used only when altUpdate==true
mse1x = 1./gam1x;
mse = nan(1,nitMax);
for i=1:nitMax
  if ~altUpdate || mod(i,2) % update on odd i
    % nonlinear output stage
    [mse1z,zvar] = estOutAvg.mse(1/gam1z);
    eta1z = 1/zvar;
    gam2z = eta1z - gam1z; % = gamw in AWGN case
  end
 
  if ~altUpdate || ~mod(i,2) % update on even i
    % nonlinear input stage
    [mse1x,xvar] = estInAvg.mse(1/gam1x);
    eta1x = 1/xvar;
    gam2x = eta1x - gam1x;

    if (altUpdate&&(i>2))||(~altUpdate&&(i>1))
      gam2x = dampGam*gam2x + (1-dampGam)*gam2xold; % damp
    end
  end

  % linear stage (note length(d)=min(M,N))
  alf = (1/N)*sum(d./(d+gam2x/gam2z)) - eps;
  eta2x = gam2x./(1-alf); % reported but not used
  eta2z = del*gam2z./alf; % reported but not used
  gam1x = gam2x.*alf./(1-alf);
  gam1z = gam2z.*(del-alf)./alf;

  % record mse
  mse(i) = mse1x;

  % prepare for next iteration
  mse1xold = mse1x;
  gam2xold = gam2x;
end %i
