% Simulation parameters
ntest = 100;                    % Number of simulations per value of nz
nzTest = [50 75 100 150 200]';  % Values of nz to test
nmeasTest = length(nzTest);
simGLasso = true;    % Set to simulate GLasso.  Otherwise load from file

% Methods to test other than GLasso
GAMP_METH = 1;      % Group AMP
LS_METH = 2;        % Least squares
GLASSO_METH = 3;    % Group lasso
GOMP_METH = 4;      % Group OMP
AMP_METH = 5;       % AMP
methStr = {'gamp', 'lmmse', 'glasso', 'gomp', 'amp'};
methTest0 = [LS_METH GAMP_METH GOMP_METH AMP_METH];  
nmeth = length(methTest0);   % Total number of methods

% Initialize vectors depending on whether GLasso is being simulated
mseMethTot = cell(nmeasTest,1);
if (simGLasso)
    mseGam = cell(nmeasTest,1);
    mseOptGam = zeros(ntest,nmeasTest);
    gamOpt = zeros(nmeasTest,1);
else
    load data/grpSparseSim;
end
saveDat = true; 

% Loop over measurements to test
for imeas = 1:nmeasTest
    
    nz = nzTest(imeas);
    
    if simGLasso
        % Compute range of gamma values to test
        if (imeas == 1)
            gamRange = 10;
            gamNom = 0.01;
            ngam = 11;
        else
            gamNom = gamOpt(imeas-1);
            gamRange = 2;
            ngam = 3;
        end
        gamTest = logspace(log10(gamNom/gamRange), log10(gamNom*gamRange), ngam);        
        mseGami = zeros(ntest,ngam);
        
        % Loop over gamma values
        for igam = 1:ngam
            
            % Get gamma value
            gam = gamTest(igam);
            
            % Run group sparsity test
            methTest = [GLASSO_METH];
            rng(42); % Set seed for repeatability
            grpSparseTest;
            
            % Save value
            mseGami(:,igam) = mseMeth;
            fprintf(1,'nz=%d gam=%f mse=%f\n', nz, gam, median(mseGami(:,igam)));
        end
        
        % Find the optimal gamma
        [mm,im] = min(median(mseGami));
        gamOpt(imeas) = gamTest(im);
        mseOptGam(:,imeas) = mseGami(:,im);
        mseGam{imeas} = mseGami;
        
        % Save the results
        %mseMeth{imeas}
        if saveDat
            save data/grpSparseSim mseGam mseOptGam gamOpt nzTest;
        end
    end
    
    % Simulate the other estimators
    methTest = methTest0;
    rng(42); % Set seed for repeatability
    grpSparseTest;
    mseMethTot{imeas} = mseMeth;
    
end

if saveDat
    save data/grpSparseSim_omp_amp mseMethTot mseGam mseOptGam gamOpt nzTest methTest0;
end
