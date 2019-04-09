addpath('../main')
addpath('../stateEvo')
addpath('../classification')
addpath('../phase')

% handle random seed
if verLessThan('matlab','7.14')
  defaultStream = RandStream.getDefaultStream;
else
  defaultStream = RandStream.getGlobalStream;
end
if 0 % new RANDOM trial
  savedState = defaultStream.State;
  save random_state.mat savedState;
else % repeat last trial
  load random_state.mat
end
defaultStream.State = savedState;

% simulation parameters
likeType = 'AWGN'; % in {'AWGN','Probit','Phaseless'}
switch likeType % choose which demo to run
 case 'AWGN'
  isCmplx = false; % simulate complex-valued case?
  N = 1024; % signal dimension [1024]
  del = 0.5; % measurement rate M/N [0.5]
  beta = 0.2; % sparsity rate K/N [0.2]
  SNRdB = 40; % [40]
  L = 100;  % # of measurement vectors [100]
 case 'Probit'
  isCmplx = false; % must set =false 
  N = 2048; % signal dimension [2048 or larger to match state evo]
  del = 4.0;  % measurement rate M/N [4.0]
  beta = 1/4; % sparsity rate K/N [1/4]
  SNRdB = 40; % [40]
  L = 100;  % # of measurement vectors [100 to match state evo]
 case 'Phaseless'
  isCmplx = true; % must set =true 
  N = 256; % signal dimension [256]
  del = 12.0; % measurement rate M/N [12.0]
  beta = 0.1 % sparsity rate K/N [1.0]
  SNRdB = 40; % [40]
  L = 1  % # of measurement vectors [10]
 otherwise
  error('unrecognized likeType')
end
svType = 'spread'; % in {'cond_num','spread','low_rank'}
cond_num = 1e0; % condition number
spread = 1; % amount to spread singular values (=1 means iid Gaussian A, =0 means frame)
low_rank = round(min(N,round(del*N))/8);
UType = 'Haar'; % in {'DFT','DCT','DHT','DHTrice','Haar','I'}
VType = 'Haar'; % in {'DFT','DCT','DHT','DHTrice','Haar','I'}
shuffle = true; % shuffle rows of V' ?
randsign = true; % randomly sign-flip columns of V' ?
plot_traj = true; % plot trajectory of each column?
plot_sig = true; % plot signal values for each column?
runOracle = true; % calculate support oracle in AWGN case?
runGAMP = true; % run GAMP?
showMedian = true; % show median instead of mean

% algorithmic parameters
maxIt = 100; % max iterations for VAMP
tol = min(1e-3,max(1e-6,10^(-SNRdB/10))); % stopping tolerance for VAMP
denoiser = 'BG'; % in {'BG','DMM','MAPLaplace'}
learnPrior = false; % automatically tune the denoiser?
learnLike = false; % automatically tune the likelihood?
altUpdate = false; % alternate updates of x and z in VAMP?
switch likeType
  case 'AWGN'
    damp = 1.0; % damping parameter for mean term
    dampGam = 1.0; % damping parameter for precision term 
  case 'Probit'
    damp = 1.0; % damping parameter for mean term
    dampGam = 1.0; % damping parameter for precision term 
  case 'Phaseless'
    damp = 0.8; % damping parameter for mean term
    dampGam = 0.5; % damping parameter for precision term 
end

% other defaults
fixed_K = true; % used fixed sparsity K=E{K}=round(rho*M)?
Afro2 = N; % squared Frobenius norm of matrix
xvar0 = 1; % prior variance of x elements
xmean1 = 0; % prior mean of non-zero x coefs

% setup
M = round(del*N);
xvar1 = xvar0/beta; % prior variance of non-zero x coefs
wvar = (Afro2/M)*10^(-SNRdB/10)*beta*(abs(xmean1)^2 + xvar1); 
if (strcmp(UType,'DFT')||strcmp(VType,'DFT'))&&(~isCmplx)
  warning('setting isCmplx=true since complex-valued matrix')
  isCmplx = true;
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
%x = x(:,1)*ones(1,L); display('repeated x'); % for testing

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
A = mat.fxnA; Ah = mat.fxnAh;
U = mat.fxnU; Uh = mat.fxnUh;
V = mat.fxnV; Vh = mat.fxnVh;
d = mat.s.^2; 

% crease noisy observations
z = A(x); 
SNRdB_test = 20*log10(norm(z(:))/norm(w(:)));
switch likeType
  case 'AWGN'
    y = z + w;
  case 'Probit'
    if isCmplx
      error('Set isCmplx=false for Probit likelihood')
    else
      y = ((z+w)>0);
    end
  case 'Phaseless'
    if ~isCmplx
      error('Set isCmplx=true for Phaseless likelihood')
    else
      y = abs(z+w);
    end
  otherwise
    error('unsupported likeType')
end

% support-oracle performance bound for AWGN case
if strcmp(likeType,'AWGN')&&runOracle
  I = speye(N);
  x0 = zeros(N,L);
  oracleNMSEdB = nan(1,L);
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

% establish input denoiser
switch denoiser
case 'BG'
  if learnPrior
    betaInit = beta/10  
    xvar0init = xvar0;
    xvar1init = xvar0init/betaInit 
    tuneDim = 'joint';
    if isCmplx
      EstimIn = SparseScaEstim(CAwgnEstimIn(0,xvar1init,0,'autoTune',true,'mean0Tune',false,'tuneDim',tuneDim),betaInit,0,'autoTune',true,'tuneDim',tuneDim);
    else
      EstimIn = SparseScaEstim(AwgnEstimIn(0,xvar1init,0,'autoTune',true,'mean0Tune',false,'tuneDim',tuneDim),betaInit,0,'autoTune',true,'tuneDim',tuneDim);
    end
  else
    if isCmplx
      EstimIn = SparseScaEstim(CAwgnEstimIn(xmean1,xvar1),beta);
    else
      EstimIn = SparseScaEstim(AwgnEstimIn(xmean1,xvar1),beta);
    end
  end
  if strcmp(likeType,'Phaseless')
    % important tricks for phase retrieval
    EstimIn.estim1.var0 = 10*xvar1; % increase initial prior variance
    EstimIn.p1 = 0.999; % set initial prior sparsity to just under 1
    EstimIn.autoTune = true; % use autoTuning to find true value
    EstimIn.estim1.autoTune = true; % use autoTuning to find true value
    EstimIn.counter = 25; % don't turn on autoTuning until things stabilize
    EstimIn.estim1.counter = 25; % don't turn on autoTuning until things stabilize
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

% establish likelihood
tuneDim = 'joint';
switch likeType
  case 'AWGN'
    if learnLike
      wvarInit = 0.01 * norm(y,'fro')^2 / (M*L)  % SNR ~= -20 dB
      if isCmplx,
        EstimOut = CAwgnEstimOut(y,wvarInit,false,'autoTune',true,...
                        'tuneMethod','EM','tuneDamp',1,'tuneDim',tuneDim);
      else
        EstimOut = AwgnEstimOut(y,wvarInit,false,'autoTune',true,...
                        'tuneMethod','EM','tuneDamp',1,'tuneDim',tuneDim);
      end
    else
      if isCmplx,
        EstimOut = CAwgnEstimOut(y,wvar);
      else
        EstimOut = AwgnEstimOut(y,wvar);
      end
    end
  case 'Probit'
    if learnLike, 
      wvarInit = wvar*100 
      EstimOut = ProbitEstimOut(y,0,wvarInit,false,'autoTune',true,...
                        'tuneMethod','ML','tuneDim',tuneDim);
    else
      EstimOut = ProbitEstimOut(y,0,wvar);
    end
  case 'Phaseless'
    if learnLike, 
      error('not yet supported!')
    else
      EstimOut = ncCAwgnEstimOut(y,wvar,0,false);
    end
end

% establish debiasing/disambiguation 
switch likeType
  case 'AWGN'
    debias = @(xhat,x) xhat; % trivial
  case 'Probit'
    debias = @(xhat,x) bsxfun(@times, xhat, sum(conj(xhat).*x,1)./sum(abs(xhat).^2,1));
  case 'Phaseless'
    debias = @(xhat,x) disambig1Dfft(xhat,x);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup VAMP
vampOpt = VampGlmOpt;
vampOpt.nitMax = maxIt;
vampOpt.tol = tol; 
vampOpt.damp = damp;
vampOpt.dampGam = dampGam;
%         [z1,gam2z,p2,x1,gam2x,r2, z2,gam1z,p1,x2,gam1x,r1]
%vampOpt.dampConfig = [0,1,0,1,0,0, 1,0,0,0,1,0]; % original from VampGlmEst
%vampOpt.dampConfig = [0,1,1,1,0,0, 0,0,0,0,1,1]; % from dampTest
vampOpt.dampConfig = [0,1,1,1,0,0, 0,0,0,0,0,1]; % from dampTest
%vampOpt.dampConfig = [0,1,0,0,1,1, 0,1,1,0,0,0]; % symmetric
vampOpt.verbose = false;
vampOpt.fxnErr1 = @(x1,z1) 10*log10( sum(abs(x1-x).^2,1)./sum(abs(x).^2,1) );
vampOpt.fxnErr2 = @(x1,z1) 10*log10( sum(abs(x-debias(x1,x)).^2,1)./sum(abs(x).^2,1) );
vampOpt.Ah = Ah; vampOpt.d = d; vampOpt.N = N;
vampOpt.U = U; vampOpt.Uh = Uh; vampOpt.V = V; vampOpt.Vh = Vh;
vampOpt.altUpdate = altUpdate;
vampOpt.silent = true; % suppress warnings about negative precisions

% initialize VAMP 
zvar = beta*(abs(xmean1)^2+xvar1)*Afro2/M;
zvarTest = norm(z,'fro')^2/numel(z); % for testing
if 1
  % matches state evolution but either uses z or sets phat=0
  pvar = 1e0*zvar; % set pvar between 0 and zvar; note pvar=zvar gives phat=0
  phat = (1-pvar/zvar)*z ...
         + sqrt((1-pvar/zvar)*pvar/2)*(randn(M,L)+(1j^isCmplx)*randn(M,L));
         %+ pvar/zvar*A(beta*xmean1*ones(N,L));
else
  % uses random phat, but not consistent with state evolution
  pvar = 1e1*zvar; 
  phat = sqrt(pvar/2)*(randn(M,L)+(1j^isCmplx)*randn(M,L));
end
if 1
  % better match to state-evolution but uses true x
  rvar = 1e2*xvar0; % set rvar between 0 and inf; note rvar=0 gives rhat=x
  rhat = x + sqrt(rvar/2)*(randn(N,L)+(1j^isCmplx)*randn(N,L));
else
  % doesn't require true x but not consistent with state evolution
  rvar = 2*xvar0; % set rvar >> xvar0
  rhat = sqrt((rvar-xvar0)/2)*(randn(N,L)+(1j^isCmplx)*randn(N,L));
end
vampOpt.p1init = phat;
vampOpt.r1init = rhat;
if isCmplx, % inform SparseScaEstim that input is complex-valued
  vampOpt.r1init(1) = vampOpt.r1init(1) + eps*1i; 
end
vampOpt.gam1zinit = 1./pvar;
vampOpt.gam1xinit = 1./rvar;

% run VAMP
if 1 
  % standard approach
  [x1,vampEstFin] = VampGlmEst2(EstimIn,EstimOut,A,vampOpt); 
  vampNMSEdB_ = vampEstFin.err1;
  vampNMSEdB_debiased_ = vampEstFin.err2;
  vampNit = vampEstFin.nit;
else
  % demonstrate stopping and warm-starting; should be identical to standard approach 
  vampOpt1 = vampOpt;
  vampOpt1.nitMax = 10; % only run a few iterations
  [x1,vampEstFin1,vampOptFin1] = VampGlmEst2(EstimIn,EstimOut,A,vampOpt1);
  vampOpt2 = vampOptFin1.warmStart(vampEstFin1,'nitMax',vampOpt.nitMax); % warm start
  [x2,vampEstFin2] = VampGlmEst2(EstimIn,EstimOut,A,vampOpt2);
  vampNMSEdB_ = [vampEstFin1.err1,vampEstFin2.err1];
  vampNMSEdB_debiased_ = [vampEstFin1.err2,vampEstFin2.err2];
  vampNit = vampEstFin1.nit + vampEstFin2.nit;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup VAMP state evolution
estInAvg = EstimInAvg(EstimIn,x);
if strcmp(likeType,'AWGN')
  estOutAvg = EstimOutAvg2(EstimOut,z); % monte-carlo given z
  %clear estOutAvg; estOutAvg.mse = @(pvar,pvarTrue) deal( 1/(1/wvar+1/pvar), 1/(1/wvar+1/pvar) );  % can do this in closed-form for AWGN
else
  estOutAvg = EstimOutAvg2(EstimOut,z); % monte-carlo given z
end

% run VAMP state evolution
vampSeNMSE = VampGlmSE(estInAvg,estOutAvg,d,N,M/N,vampOpt)/(beta*(abs(xmean1)^2+xvar1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print learned parameters
if learnPrior || strcmp(likeType,'Phaseless')
  beta
  betaEstimate = EstimIn.p1(1,:)
  xvar1
  xvar1estimate = EstimIn.estim1.var0(1,:)
end
if learnLike
  wvar
  switch likeType
    case 'AWGN'
      wvarEstimate = EstimOut.wvar(1) 
    case 'Probit'
      wvarEstimate = EstimOut.Var(1) 
    case 'Phaseless'
      wvarEstimate = EstimOut.var0(1) 
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup and run GAMP
gampNit = 0;
if runGAMP
  %Agamp = MatrixLinTrans(A);
  Agamp = FxnhandleLinTrans(M,N,A,Ah,Afro2/(M*N));
  switch likeType
    case 'AWGN'
      optGAMP = GampOpt('legacyOut',false,'uniformVariance',true,'adaptStepBethe',true,'step',1.0,'stepIncr',1.05,'stepWindow',2,'stepMax',1.0,'tol',tol/10,'nit',500,'xvar0',beta*(abs(xmean1)^2+xvar1));
    case 'Probit'
      optGAMP = GampOpt('legacyOut',false,'uniformVariance',true,'adaptStepBethe',true,'step',0.2,'stepIncr',1.05,'stepWindow',1,'stepMax',0.5,'stepMin',0.02,'tol',tol/100,'nit',500,'xvar0',beta*(abs(xmean1)^2+xvar1));
    case 'Phaseless'
      optGAMP = GampOpt('legacyOut',false,'uniformVariance',true,'adaptStep',false,'step',0.2,'stepMax',0.2,'tol',tol/10,'nit',500,'xvar0',beta*(abs(xmean1)^2+xvar1));
  end
  % reset prior if needed, since VAMP may have changed it
  if learnPrior
    EstimIn.p1 = betaInit;
    EstimIn.estim1.var0 = xvar1init;
    EstimIn.estim1.mean0 = 0;
  end
  % reset likelihood if needed, since VAMP may have changed it
  if learnLike
    switch likeType
      case 'AWGN'
        EstimOut.wvar = wvarInit;
        EstimOut.tuneMethod = 'ML';
        EstimOut.tuneDim = 'col'; % seems to be important
      case 'Probit'
        EstimOut.Var = wvarInit;
        EstimOut.tuneMethod = 'ML';
        warning('NEED TO SET tuneDim=col')
      case 'Phaseless'
        EstimOut.var0 = wvarInit;
    end
  end
  % some tricks for phase retrieval
  if strcmp(likeType,'Phaseless')
    EstimIn.p1 = 0.999; % set initial prior sparsity to just under 1
    EstimIn.estim1.var0 = 10*xvar1; % increase initial prior variance
    EstimIn.autoTune = true; % use autoTuning to find true value
    EstimIn.estim1.autoTune = true; % use autoTuning to find true value
    EstimIn.counter = 25; % don't turn on autoTuning until things stabilize
    EstimIn.estim1.counter = 25; % don't turn on autoTuning until things stabilize
  end
  % initialize
  [xhat,xvar] = EstimIn.estim(rhat,rvar);
  optGAMP.xhat0 = xhat;
  optGAMP.xvar0 = xvar;
  % run GAMP
  tstart = tic;
  [gampEstFin,optGampFin,gampEstHist] = gampEst(EstimIn,EstimOut,Agamp,optGAMP);
  time_gamp = toc(tstart);
  gampNit = gampEstFin.nit;
  gampNMSEdB_ = nan(L,gampNit);
  gampNMSEdB_debiased_ = nan(L,gampNit);
  for l=1:L
    xhat_ = gampEstHist.xhat((l-1)*N+[1:N],:);
    gampNMSEdB_(l,:) = 10*log10(sum(abs(xhat_-x(:,l)*ones(1,gampNit)).^2,1)/norm(x(:,l))^2);
    gampNMSEdB_debiased_(l,:) = 10*log10(sum(abs( debias(xhat_,x(:,l)*ones(1,gampNit))-x(:,l)*ones(1,gampNit)).^2,1)/norm(x(:,l))^2);
  end
  figure(4); clf; gampShowHist(gampEstHist,optGampFin,x); % debug GAMP
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot results
figure(1); clf;
  if showMedian
    vampNMSE_avg = median(10.^(vampNMSEdB_debiased_/10),1);
    if runGAMP, gampNMSE_avg = median(10.^(gampNMSEdB_debiased_/10),1); end;
  else
    vampNMSE_avg = mean(10.^(vampNMSEdB_debiased_/10),1);
    if runGAMP, gampNMSE_avg = mean(10.^(gampNMSEdB_debiased_/10),1); end;
  end
  plot(1:vampNit,vampNMSE_avg,'.-','Displayname','VAMP');
  set(gca,'YScale','log','XScale','log')
  hold on;
    if runGAMP, semilogx(1:gampNit,gampNMSE_avg,'.-','Displayname','GAMP'); end;
    semilogx(1:maxIt,vampSeNMSE,'k--','Displayname','VAMP-SE');
  hold off;
  axe = axis;
  axis([1,axe(2),10^floor(log10(max(min([vampNMSE_avg,vampSeNMSE]),1e-15))),100])
  legend(gca,'show')
  grid on
  xlabel('iteration')
  if showMedian
    ylabel('median debiased NMSE')
  else
    ylabel('avg debiased NMSE')
  end
  tit_str = [likeType,...
             ', N=',num2str(N),', del=',num2str(del),', beta=',num2str(beta),...
             ', damp=',num2str(damp),', dampG=',num2str(dampGam),...
             ', SNRdB=',num2str(SNRdB),', L=',num2str(L)];
  title(tit_str);

if plot_traj
figure(2); clf;
subplot(211) 
  % plot VAMP
  handy = plot(1:vampNit-1,vampNMSEdB_(:,2:end).','b.-'); % 1st iteration is trivial
  set(handy(1),'Displayname','VAMP') 
  % plot GAMP
  if runGAMP
    hold on; 
      handy = [handy, plot(1:gampNit-1,gampNMSEdB_(:,2:end),'r.-')];
      set(handy(1,end),'Displayname','GAMP') 
    hold off; 
  end
  % plot support oracle
  if strcmp(likeType,'AWGN')&&runOracle
    ax = gca; ax.ColorOrderIndex = 1; % use same colors
    hold on; 
      handy = [handy, plot([1;max(vampNit,gampNit)],[1;1]*oracleNMSEdB,'k--')]; 
      set(handy(1,end),'Displayname','oracle') 
    hold off; 
  end
  legend(handy(1,:))
  ylabel('NMSE [dB]')
  xlabel('iterations')
  grid on

subplot(212) 
  % plot VAMP
  handy = plot(1:vampNit-1,vampNMSEdB_debiased_(:,2:end).','b.-'); % 1st iteration trivial
  set(handy(1),'Displayname','VAMP') 
  % plot GAMP
  if runGAMP
    hold on; 
      handy = [handy, plot(1:gampNit-1,gampNMSEdB_debiased_(:,2:end),'r.-')];
      set(handy(1,end),'Displayname','GAMP') 
    hold off; 
  end
  % plot support oracle
  if strcmp(likeType,'AWGN')&&runOracle
    ax = gca; ax.ColorOrderIndex = 1; % use same colors
    hold on; 
      handy = [handy, plot([1;max(vampNit,gampNit)],[1;1]*oracleNMSEdB,'k--')]; 
      set(handy(1,end),'Displayname','oracle') 
    hold off; 
  end
  legend(handy(1,:))
  ylabel('debiased NMSE [dB]')
  xlabel('iterations')
  grid on
end % plot_traj

if plot_sig
  figure(3); clf;
  l = 1; % choose which column to plot
  if runGAMP, subplot(211); end;
    stem(real(x1(:,l)))
    hold on; 
      stem(real(debias(x1(:,l),x(:,l))),'x'); 
      stem(real(x(:,l)),'--'); 
    hold off;
    legend('VAMP','VAMP debiased','true')
    if L>1, title(['column ',num2str(l),' of ',num2str(L)]); end;
    xlabel('coefficient index')
    grid on;
    if isCmplx, ylabel('real part'); end
  if runGAMP,
  subplot(212)
    xg = gampEstFin.xhat;
    stem(real(xg(:,l)))
    hold on; 
      stem(real(debias(xg(:,l),x(:,l))),'x'); 
      stem(real(x(:,l)),'--'); 
    hold off;
    legend('GAMP','GAMP debiased','true')
    xlabel('coefficient index')
    grid on;
    if isCmplx, ylabel('real part'); end
  end
end
