%% Script comparing the phase retrieval performance of prGAMP and prVAMP with different types of iid measurement matrices
%The variable matrix_type chooses the matrix type

clear;clc;

%%Problem Settings
%Set up the signal
n = 16^2; % signal dimension [16^2]
p = 12*n; % number of phaseless measurements [12*n]
k = round(0.999*n); % signal sparsity
x_o = zeros(n,1); x_o(randperm(n,k)) = randn(k,1)+1i*randn(k,1);
xvar_nz = var(x_o(x_o~=0)); % variance of nonzero coefs
xvar = var(x_o); % variance of all coefs

%Set up the measurements
matrix_type = 3;
switch matrix_type
  case 1
    A = round(rand(p,n));% 0,1 Measurements
  case 2
    A = 2*round(rand(p,n))-1;% -1,1 Measurements
  case 3
    A = randn(p,n)+1i*randn(p,n);% Gaussian Measurements
end
sigma_w = 1e-8; % 1e-8
noise_vec = sigma_w*(1/sqrt(2)*randn(p,1)+1i*1/sqrt(2)*randn(p,1));
z_o = A*x_o;
y = abs(z_o+noise_vec);

%%Algorithmic Settings
%Set up the prior
xvar_nz_init = 10*xvar_nz; % iniital nonzero-coef variance (set >> xvar_nz !)
spars_init = 0.999; % initial sparsity rate (set near 1 !)
tuneDelayGamp = 25; % number of iterations to wait before tuning
tuneDelayVamp = 25; % number of iterations to wait before tuning
EstimIn = SparseScaEstim(...
           CAwgnEstimIn(0,xvar_nz_init,false,...
                        'autoTune',true,...
                        'mean0Tune',false,...
                        'counter',tuneDelayGamp),...
           spars_init,false,'autoTune',true,'counter',tuneDelayGamp);

%Set up the likelihood function
wvar_init = sigma_w^2;
EstimOut = ncCAwgnEstimOut(y,wvar_init*ones(p,1),false,false);

%Set up the initialization
xvar_init = xvar;
x_init = sqrt(xvar_init)*(1/sqrt(2)*randn(n,1)+1i*1/sqrt(2)*randn(n,1)); 

%%Simulation Settings
plot_on = true; % plot NMSE versus iteration? (slows things down!)

%%prGAMP
%Set up prGAMP
gampOpt = GampOpt();
gampOpt.legacyOut = false;
gampOpt.xhat0 = x_init;
gampOpt.xvar0 = xvar_init;
gampOpt.shat0 = eps*ones(p,1);
gampOpt.svar0 = eps;
gampOpt.pvarMin = eps;
gampOpt.adaptStep = false; 
gampOpt.step = .2; 
gampOpt.stepMax = .2; %1 means no damping
gampOpt.nit = 500;
gampOpt.tol = 1e-6;
gampOpt.uniformVariance = true;
gampOpt.pvarMin = 1e-11;
gampOpt.xvarMin = 1e-11;
%gampOpt.errFcn1=@(x,z,r,p) norm(x_o-disambig1Dfft(x,x_o))^2/norm(x_o)^2;
%gampOpt.errFcn2=@(x,z,r,p) norm(y-abs(z))^2/norm(y)^2;

%Run prGAMP
t0=cputime;
if ~plot_on
  gampFin = gampEst(EstimIn,EstimOut,A,gampOpt);
else
  [gampFin,~,gampHist] = gampEst(EstimIn,EstimOut,A,gampOpt);
end
x_prGAMP = gampFin.xhat(:);
t_prGAMP = cputime-t0;
prGAMP_NMSE = norm(x_o-disambig1Dfft(x_prGAMP,x_o))^2/norm(x_o)^2;

%Print results
display(['prGAMP Reconstruction NMSE=',num2str(prGAMP_NMSE),', time=',num2str(t_prGAMP),', sparsity_est=',num2str(EstimIn.p1)]);

%%prVAMP
% Reset prior since internal parameters changed while running prGAMP
EstimIn = SparseScaEstim(CAwgnEstimIn(0,xvar_nz_init,false,'autoTune',true,'mean0Tune',false,'counter',tuneDelayVamp),spars_init,0,'autoTune',true,'counter',tuneDelayVamp);

%Set up prVAMP
[U,S,V]=svd(A,'econ');
d=diag(S).^2;
vampOpt = VampGlmOpt;
vampOpt.nitMax = 500;
vampOpt.tol = 1e-6;
vampOpt.damp = 0.8; % try 0.8; 1 means no damping
vampOpt.dampGam = 0.5; % try 0.5; 1 means no damping
%vampOpt.dampConfig = [0,1,0,1,0,0, 1,0,0,0,1,0]; % original from VampGlmEst
%vampOpt.dampConfig = [0,1,1,1,0,0, 0,0,0,1,0,1]; % best from dampTest
vampOpt.dampConfig = [0,1,1,1,0,0, 0,0,0,0,0,1]; % best from dampTest
vampOpt.verbose = false;
vampOpt.U = U;
vampOpt.V = V;
vampOpt.d = d;
vampOpt.r1init = x_init;
vampOpt.gam1xinit = 1/xvar_init;
vampOpt.silent = true;
vampOpt.altUpdate = false;
%vampOpt.fxnErr1 = @(x,z) norm(x_o-disambig1Dfft(x,x_o))^2/norm(x_o)^2;
%vampOpt.fxnErr2 = @(x,z) norm(y-abs(z))^2/norm(y)^2;

%Run prVAMP
clear vampHist;
t0=cputime;
if ~plot_on
  [x_prVAMP,vampFin] = VampGlmEst2(EstimIn,EstimOut,A,vampOpt); 
else
  [x_prVAMP,vampFin,optFin,vampHist] = VampGlmEst2(EstimIn,EstimOut,A,vampOpt); % slow
end
t_prVAMP = cputime-t0;
x_prVAMP = x_prVAMP(:);
prVAMP_NMSE = norm(x_o-disambig1Dfft(x_prVAMP,x_o))^2/norm(x_o)^2;

% Print results
display(['prVAMP Reconstruction NMSE=',num2str(prVAMP_NMSE),', time=',num2str(t_prVAMP),', sparsity_est=',num2str(EstimIn.p1)]);

%%Plot results
if plot_on
  %Compute disambiguated error 
  vampErrX = sum(abs(disambig1Dfft(squeeze(vampHist.x1),x_o*ones(1,vampFin.nit))-x_o*ones(1,vampFin.nit)).^2,1)/norm(x_o)^2;
  vampErrZ = sum((abs(squeeze(vampHist.z1))-y*ones(1,vampFin.nit)).^2)/norm(y)^2;
  gampErrX = sum(abs(disambig1Dfft(gampHist.xhat,x_o*ones(1,gampFin.nit))-x_o*ones(1,gampFin.nit)).^2,1)/norm(x_o)^2;
  gampErrZ = sum((abs(gampHist.zhat)-y*ones(1,gampFin.nit)).^2)/norm(y)^2;
  %gampErrX = gampFin.err1';
  %gampErrZ = gampFin.err2';
  %vampErrX = vampFin.err1;
  %vampErrZ = vampFin.err2;

  %Plot GAMP and VAMP's NMSEs versus iteration
  figure(1); clf
  plot([1:vampFin.nit],10*log10(vampErrX),...
       [1:gampFin.nit],10*log10(gampErrX),...
       [1:vampFin.nit],10*log10(vampErrZ),'--',...
       [1:gampFin.nit],10*log10(gampErrZ),'--')
  axe = axis;
  axis([0,max(vampFin.nit,gampFin.nit),...
        min(10*log10([vampErrX,gampErrX]))-5,...
        min(axe(4),50)])
  legend('VAMP: x1','GAMP: x1','VAMP: |z1|','GAMP: |z1|')
  ylabel('NMSE [dB]')
  xlabel('iteration')
  title('GAMP and VAMP')
  grid on

  %Plot VAMP's precisions versus iteration
  figure(2); clf
  semilogy([1:vampFin.nit],vampHist.gam1x,...
           [1:vampFin.nit],vampHist.gam1z,...
           [1:vampFin.nit],vampHist.gam2x,...
           [1:vampFin.nit],vampHist.gam2z)
  legend('gam1x','gam1z','gam2x','gam2z','Location','Northwest')
  xlabel('iteration')
  ylabel('precision')
  title('VAMP')
  grid on
end
