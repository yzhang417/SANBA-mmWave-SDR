addpath('../main')
addpath('../stateEvo')

% handle random seed
if verLessThan('matlab','7.14')
  defaultStream = RandStream.getDefaultStream;
else
  defaultStream = RandStream.getGlobalStream;
end
if 1 % new RANDOM trial
  savedState = defaultStream.State;
  save random_state.mat savedState;
else % repeat last trial
  load random_state.mat
end
defaultStream.State = savedState;

% simulation parameters
L = 10; % # of measurement vectors [10]
SNRdB = 40; % [40]
N = 1024; % signal dimension [1024]
del = 0.5; % measurement rate M/N [0.5]
rho = 0.2; % normalized sparsity rate E{K}/M=(E{K}/N)/del [0.2]
svType = 'cond_num'; % in {'cond_num','spread','low_rank'}
cond_num = 1; % condition number [50 is limit for damped GAMP]
spread = 1; % amount to spread singular values (=1 means iid Gaussian A, =0 means frame) [6 is limit for GAMP]
low_rank = round(min(N,round(del*N))/2);
UType = 'I'; % in {'DFT','DCT','DHT','DHTrice','Haar','I'}
VType = 'Haar'; % in {'DFT','DCT','DHT','DHTrice','Haar','I'}
shuffle = true; % shuffle rows of V' ?
randsign = true; % randomly sign-flip columns of V' ?
isCmplx = false; % simulate complex-valued case?
plot_traj = true; % plot trajectory of each column?
median_on = true; % use median instead of mean?
runOracle = true; % calculate support oracle?
runEMBGAMP = true; % run GAMP?
verbose = false;

% algorithmic parameters
maxit = 50; % max iterations for VAMP
tol = min(1e-3,max(1e-6,10^(-SNRdB/10))); % stopping tolerance for VAMP
damp = 0.95; % damping for VAMP
denoiser = 'BG'; % in {'BG','DMM','MAPLaplace'}
learnPrior = false; % automatically tune the denoiser?
learnNoisePrec = false; % automatically tune the noise variance?

% other defaults
fixed_K = false; % used fixed sparsity K=E{K}=round(rho*M)?
Afro2 = N; % squared Frobenius norm of matrix
xvar0 = 1; % prior variance of x coefs
xmean1 = 0; % prior mean of non-zero x coefs

% setup
M = round(del*N);
beta = rho*M/N; % probability of a non-zero coef
xvar1 = xvar0/beta; % prior variance of non-zero x coefs
wvar = (Afro2/M)*10^(-SNRdB/10)*beta*(abs(xmean1)^2+xvar1); 
if (strcmp(UType,'DFT')||strcmp(VType,'DFT'))&&(~isCmplx)
  warning('setting isCmplx=true since complex-valued matrix')
  isCmplx = true;
elseif isCmplx&&(~strcmp(UType,'DFT'))&&(~strcmp(VType,'DFT'))
  warning('setting isCmplx=false since real-valued matrix')
  isCmplx = false;
end

% generate signal 
x = zeros(N,L);
for l=1:L
  if fixed_K
    supp = randperm(N,round(beta*N)); 
  else
    supp = find(rand(N,1)<beta); 
  end
  K = length(supp);
  %supp = 1:K; display('block support'); % for testing
  if isCmplx
    x(supp,l) = xmean1 + sqrt(0.5*xvar1)*randn(K,2)*[1;1j];
  else
    x(supp,l) = xmean1 + sqrt(xvar1)*randn(K,1);
  end
end
%x = abs(x); display('positive x')
%x =x(:,1)*ones(1,L); display('repeated x')

% generate noise 
if isCmplx
  w = sqrt(0.5*wvar)*(randn(M,L) + 1j*randn(M,L));
else
  w = sqrt(wvar)*randn(M,L);
end

% generate linear transform
switch svType
  case 'spread', svParam = spread;
  case 'cond_num', svParam = cond_num;
  case 'low_rank', svParam = low_rank;
end
mat = genMatSVD(M,N,UType,svType,svParam,VType,...
                'isCmplx',isCmplx,'Afro2',Afro2,...
                'shuffle',shuffle,'randsign',randsign,...
                'fxnHandles',true);
d = [mat.s.^2;zeros(M-length(mat.s),1)]; % need length(d)=M
A = mat.fxnA; Ah = mat.fxnAh;
U = mat.fxnU; Uh = mat.fxnUh;


% generate observation
z = A(x); 
SNRdB_test = 20*log10(norm(z(:))/norm(w(:)));
y = z + w;

% support-oracle performance bound
if runOracle
  I = speye(N);
  x0 = zeros(N,L);
  oracleNMSEdB = nan(L,1);
  for l=1:L
    supp = find(x(:,l)~=0);
    try
      A0 = A(I(:,supp)); % fast but not compatible with all A(.)
    catch
      K = length(supp); % slow but compatible with all A(.)
      A0 = zeros(M,K);
      for k=1:K, A0(:,k) = A([zeros(supp(k)-1,1);1;zeros(N-supp(k),1)]); end
    end
    a0 = A0*ones(length(supp),1);
    x0(supp,l) = xmean1 + (A0'*A0+(wvar/xvar1)*eye(length(supp)))\(A0'*(y(:,l)-a0*xmean1)); 
    oracleNMSEdB(l) = 20*log10(norm(x0(:,l)-x(:,l))/norm(x(:,l)));
  end
end

% establish denoiser
switch denoiser
case 'BG'
  if learnPrior
    betaInit = 1/N;
    xvar0init = xvar0;
    xvar1init = xvar0init/betaInit;
    tuneDim = 'col';
    if isCmplx
      if beta<1
        EstimIn = SparseScaEstim(CAwgnEstimIn(0,xvar1init,0,'autoTune',true,'tuneDim',tuneDim),betaInit,0,'autoTune',true,'tuneDim',tuneDim);
      elseif beta==1
        EstimIn = CAwgnEstimIn(0,xvar1init,0,'autoTune',true,'tuneDim',tuneDim);
      else
        error('invalid rho since rho>N/M')
      end
    else
      if beta<1
        EstimIn = SparseScaEstim(AwgnEstimIn(0,xvar1init,0,'autoTune',true,'tuneDim',tuneDim),betaInit,0,'autoTune',true,'tuneDim',tuneDim);
      elseif beta==1,
        EstimIn = AwgnEstimIn(0,xvar1init,0,'autoTune',true,'tuneDim',tuneDim);
      else
        error('invalid rho since rho>N/M')
      end
    end
  else
    if isCmplx
      EstimIn = SparseScaEstim(CAwgnEstimIn(xmean1,xvar1),beta);
    else
      EstimIn = SparseScaEstim(AwgnEstimIn(xmean1,xvar1),beta);
    end
  end
case 'DMM'
  alpha = 1.5;
  debias = false;
  EstimIn = SoftThreshDMMEstimIn(alpha,'debias',debias);
  if learnPrior, 
    warning('learnPrior not implemented for SoftThreshDMM'); 
  end;
case 'MAPLaplace'
  lam = 1/sqrt(wvar);
  if learnPrior,
    EstimIn = SoftThreshEstimIn(lam,0,'autoTune',true,'counter',10) 
  else
    EstimIn = SoftThreshEstimIn(lam);
  end
otherwise
  error('unknown denoiser')
end

% setup VAMP
vampOpt = VampSlmOpt;
vampOpt.nitMax = maxit;
vampOpt.tol = tol;
vampOpt.damp = damp;
vampOpt.learnNoisePrec = learnNoisePrec;
if learnNoisePrec
  vampOpt.learnNoisePrec = true;
else
  vampOpt.learnNoisePrec = false;
  vampOpt.NoisePrecInit = 1/wvar; 
end
vampOpt.learnNoisePrecAlg = 'EM';
vampOpt.learnNoisePrecNit = 100;
vampOpt.learnNoisePrecTol = 0.01;
vampOpt.verbose = verbose;
vampOpt.fxnErr = @(x2) 10*log10( sum(abs(x2-x).^2,1)./sum(abs(x).^2,1) ); 
vampOpt.Ah = Ah; vampOpt.d = d; vampOpt.N = N;
vampOpt.U = U; vampOpt.Uh = Uh;

% run VAMP
[~,vampEstFin] = VampSlmEst(EstimIn,y,A,vampOpt);
vampNMSEdB_ = vampEstFin.err; 
vampNit = vampEstFin.nit;

% run VAMP state evolution
estInAvg = EstimInAvg(EstimIn,x);
gam1init = 1/(beta*(abs(xmean1).^2+xvar1) + wvar*sum(d)/sum(d.^2)); % VampSlm default
%gam1init = 1000 % to test replica, try starting from near-perfect initialization
vampSeNMSE = VampSlmSE(estInAvg,d,N,wvar,vampNit,gam1init,damp)./(beta*(abs(xmean1)^2+xvar1));

% setup and run EM-BG-GAMP
gampNit = 0;
if runEMBGAMP
  Agamp = FxnhandleLinTrans(M,N,A,Ah,Afro2/(M*N));
  clear optEM optGAMP;
  optEM.heavy_tailed = false;
  optEM.robust_gamp = true;
  optGAMP.removeMean = false;
  tstart = tic;
  [~,EMfin,gampEstHist,~,optGAMPfin] = EMBGAMP(y,Agamp,optEM,optGAMP);
  time_gamp = toc(tstart);
  gampNit = length(gampEstHist.it);
  gampNMSEdB_ = nan(L,gampNit);
  for l=1:L
    gampNMSEdB_(l,:) = 10*log10(sum(abs(gampEstHist.xhat((l-1)*N+[1:N],:)-x(:,l)*ones(1,gampNit)).^2,1)/norm(x(:,l))^2);
  end
  %figure(2); clf; gampShowHist(gampEstHist,optGAMPfin,x); % debug GAMP
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot results
vstr = denoiser; if learnPrior || learnNoisePrec, vstr = ['EM-',vstr]; end;
gstr = 'BG-AMP'; if learnPrior || learnNoisePrec, gstr = ['EM-',gstr]; end;
if plot_traj
  figure(1); clf;
  % plot VAMP
  handy = semilogx(1:vampNit,vampNMSEdB_,'b.-');
  set(handy(1),'Displayname',[vstr,'-VAMP']);
  % plot GAMP
  if runEMBGAMP
    hold on;
      handy = [handy, semilogx(1:gampNit,gampNMSEdB_,'r.-')];
    hold off;
    set(handy(1,end),'Displayname',gstr);
  end
  % plot support oracle
  if runOracle
    ax = gca; ax.ColorOrderIndex = 1; % use same colors
    hold on; 
      handy = [handy, semilogx([1;max(vampNit,gampNit)],oracleNMSEdB*[1,1],'k--')]; 
    hold off;
    set(handy(1,end),'Displayname','oracle');
  else
    oracleNMSEdB = inf; % used below in axis 
  end
  % legend
  legend(handy(1,:)); 
  if median_on
    ylabel('median NMSE [dB]')
  else
    ylabel('average NMSE [dB]')
  end
  xlabel('iterations')
  grid on
  axis([1,max(vampNit,gampNit),5*floor(min([vampNMSEdB_(:);oracleNMSEdB])/5),1])
end % plot_traj

% plot state evolution
figure(2); clf;
if median_on
  vampNMSE_avg = median(10.^(vampNMSEdB_/10),1);
  oracleNMSE_avg = median(10.^(oracleNMSEdB/10),1); 
  if runEMBGAMP, gampNMSE_avg = median(10.^(gampNMSEdB_/10),1); end
else
  vampNMSE_avg = mean(10.^(vampNMSEdB_/10),1);
  oracleNMSE_avg = mean(10.^(oracleNMSEdB/10),1); 
  if runEMBGAMP, gampNMSE_avg = mean(10.^(gampNMSEdB_/10),1); end
end
plot(1:vampNit,vampNMSE_avg,'+-','Displayname',[vstr,'-VAMP']);
set(gca,'YScale','log','XScale','log')
hold on;
  semilogx(1:vampNit,vampSeNMSE,'o-','Displayname','VAMP SE'); 
  if runEMBGAMP 
    semilogx(1:gampNit,gampNMSE_avg,'x-','Displayname',gstr); 
    if runOracle
      semilogx([1,gampNit],[1,1]*oracleNMSE_avg,'-.','Displayname','oracle');
    end
  else
    if runOracle
      semilogx([1,vampNit],[1,1]*oracleNMSE_avg,'-.','Displayname','oracle');
    end
  end
hold off;
legend(gca,'show')
if runEMBGAMP
  axis([1,gampNit,10^floor(log10(min([vampNMSE_avg,oracleNMSE_avg,vampSeNMSE]))),1])
else
  axis([1,vampNit,10^floor(log10(min([vampNMSE_avg,oracleNMSE_avg,vampSeNMSE]))),1])
end
grid on
xlabel('iteration')
if median_on
  ylabel('median NMSE [dB]')
else
  ylabel('average NMSE [dB]')
end
