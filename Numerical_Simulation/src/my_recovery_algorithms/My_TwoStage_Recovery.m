%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2019 Yi Zhang and The University of Texas at Austin 
%  
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including 
% without limitation the rights to use, copy, modify, merge, publish, 
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the 
% following conditions:
% 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you use this code or any (modified) part of it in any publication,
% please cite:
%
% Yi Zhang, Kartik Patel, Sanjay Shakkottai, and Robert W. Heath Jr.. 2019. 
% Side-information-aided Non-coherent Beam Alignment Design for Millimeter 
% Wave Systems. In MobiHoc '19: The Twentieth ACM International Symposium 
% on Mobile Ad Hoc Networking and Computing, July 02-05, 2019, Catania, 
% Italy. ACM, New York, NY, USA, 10 pages.
%
% Author: Yi Zhang
% Contact email: yi.zhang.cn@utexas.edu 
% Last modified: Apr. 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function description:
% This function uses the proposed two-stage algorithm to solve the 
% compressive phase retrieval problem. For details, see
%   Yi Zhang, Kartik Patel, Sanjay Shakkottai, and Robert W. Heath Jr.. 2019. 
%   Side-information-aided Non-coherent Beam Alignment Design for Millimeter 
%   Wave Systems. In Proceedings of ACM the Twentieth International Symposium
%   on Mobile Ad Hoc Networking and Computing (MobiHoc'19). ACM, New York,
%   NY, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% measurements: received noised measurements (column vector).
% measurementMat_Input: sensing matrix.
% s: level of sparsity.
% noise_power: power of noise.
% plot_flag: whether show the recovery performance of the sparse channel 
% vector.
% Method: a structure array that groups related parameters on the recovery
% algorithms that need to be tested.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output arguments:
% recoveredSig_PLOMP: recovered sparse signal with OMP in the second stage
% (column vector).
% recoveredSig_PLGAMP: recovered sparse signal with GAMP in the second 
% stage (column vector).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [recoveredSig_PLOMP, recoveredSig_PLGAMP] = ...
            My_TwoStage_Recovery(measurements, measurementMat_Input, s, noise_power, plot_flag, Method)
    %% Problem parameters
    [m, n] = size(measurementMat_Input);
    mCS_required = round(1.75*s*log(n/s));   % No. of CS measurements

    %% Creat P and C Matrix Via SVD Decomposition
    [U,S,V] = svd(measurementMat_Input,0);
    dS = diag(S);
    mCS = min(mCS_required,length(dS))-1;
    Var_ratio = sum(dS(1:mCS))/sum(dS);
    EigThreshold = 0.80;
    while Var_ratio < EigThreshold && mCS < length(dS)
        mCS = mCS + 1;
        Var_ratio = sum(dS(1:mCS))/sum(dS);
    end
    while round(1.75*mCS*log(mCS)) < m && mCS < length(dS)
    %while mCS < 4*m && mCS < length(dS)
        mCS = mCS + 1;
    end
    Var_ratio = sum(dS(1:mCS))/sum(dS);
    U = U(:,1:mCS);
    S = S(1:mCS,1:mCS);
    V = V(:,1:mCS);
    P = U*sqrt(S);
    C = sqrt(S)*V';

    %% Measurement Matrix
    measurementMat.C = C;
    measurementMat.P = P;

    %% Print out parameters to standard output
    if plot_flag
        fprintf('\n------- Two-Stage Recovery-------\n');
        fprintf('Energy Captured By SVD Approximation: %f\n', Var_ratio);
        fprintf('Problem dimension - %d\n', n );
        fprintf('Sparsity - %d\n', s );
        fprintf('No. of measurements - %d\n', m );
        fprintf('CS problem dimension - %d\n', mCS );
        fprintf('\n------- Two-Stage Recovery-------\n');
    end

    %% Step I with PhaseLift
    % PhaseLift
    % Solve the Phase Retrieval problem using m measurements, obtaining an mCS length intermediate vector
    % Use TFOCS to solve the problem
    % Set TFOCS parameters
    opts.maxIts = 4e3;
    opts.tol = max(noise_power*1e-3,1e-10);
    opts.tol = 1e-10;
    opts.restart = 200;
    opts.printEvery = 0;
    if n > 2000
        opts.largescale = true;
    end

    % Regularization parameter
    lambda = 0.05;

    % Initial guess
    initGuess = zeros(mCS);

    % TFOCS requires special initialization of the measurement model
    % We use the Phase Retrieval matrix, P
    tfocsLinOp = initializeLinopPR( {[mCS, mCS], m}, measurementMat.P ); 

    recoveredMat = solver_TraceLS( tfocsLinOp, measurements, ...
                                                lambda, initGuess, opts );

    % Intermediate solution
    % Choose intermediate solution to be (scaled) leading eigenvector
    % Use eigs for large problems ...
    try
        [intSoln, eVal] = eig(recoveredMat);
    catch
        [intSoln, eVal] = eigs(recoveredMat);
    end
    intSoln = sqrt(eVal(end))*intSoln(:,end);
    intSoln_PL = intSoln;
    
    %% Step 2 with OMP
    if Method.PLOMP
        A_OMP = measurementMat.C;
        recoveredSig_PLOMP = OMP(A_OMP, intSoln_PL, 1e-12);
    else
        recoveredSig_PLOMP = zeros(n,1);
    end

    %% Step 2 with EMBGAMP
    if Method.PLGAMP
        A_Gamp = measurementMat.C;
        optEM.SNRdB = 10*log10(1/noise_power);
        optEM.maxEMiter = 200;
        optEM.cmplx_in = true;
        optEM.cmplx_out = true;
        optEM.lambda = s/n;
        optEM.learn_lambda = false;
        optEM.robust_gamp = true;
        try
            recoveredSig_PLGAMP = EMBGAMP(intSoln_PL,A_Gamp,optEM);
        catch
            fprintf('\nFail of EMBGAMP in step II of two-stage recovery\n');
            recoveredSig_PLGAMP = OMP(A_Gamp, intSoln_PL, 1e-12);
        end
    else
        recoveredSig_PLGAMP = zeros(n,1);
    end
end


