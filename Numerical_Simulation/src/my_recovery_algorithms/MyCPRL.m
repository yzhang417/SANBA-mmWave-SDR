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
% This function uses Compressive Phase Retrieval via Lifting (CPRL) solve 
% the typical compressive phase retrieval problem.
% It uses the CVX optimization software package. For details, see
%   Compressive Phase Retrieval From Squared Output Measurements Via 
%   Semidefinite Programming
%   Ohlsson, Henrik; Yang, Allen Y.; Dong, Roy; Shankar Sastry, S.
%   eprint arXiv:1111.6323
%   2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% measurements: received noised measurements (column vector).
% measurementMat: sensing matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output arguments:
% recoveredSig: recovered sparse signal (column vector).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function recoveredSig = MyCPRL(measurements, measurementMat)
%% Problem parameters
[m, n] = size(measurementMat);
A = measurementMat;
y = measurements;

%% CPRL
%   Solve the Sparse Phase Retrieval problem using
%       min      || b - A(X) ||_1 + mu*|| X ||_1
%   subject to          X >= 0
% Note: We are using a variant of PhaseLift for noisy data
% Standard PhaseLift problem formulation is
%       min             Tr(X)
%   subject to         A(X) = b
%                       X >= 0
% 
% For noisy data, [1] recommends solving
%       min         || b - A(X) ||_1
%   subject to          X >= 0
% The CPRL formulation we solve is based on this implementation of
% PhaseLift
%
% [1]   Solving Quadratic Equations via PhaseLift when There Are About As 
%       Many Equations As Unknowns
%       Emmanuel J. Candes and Xiaodong Li
%       http://arxiv.org/abs/1208.6247v2

% Regularization parameter
% Change this with added noise level
mu = 5e-2;
cvx_begin sdp quiet
    cvx_solver mosek
    %cvx_solver_settings('MSK_DPAR_INTPNT_CO_TOL_PFEAS', 1.0e-4)
    %cvx_solver_settings('MSK_DPAR_INTPNT_CO_TOL_DFEAS', 1.0e-4)
    %cvx_solver_settings('MSK_DPAR_INTPNT_CO_TOL_REL_GAP', 1.0e-6);
    cvx_precision low
    variable Z(n,n) hermitian
    minimize norm( diag(A*Z*A')-y,1) + mu*norm(Z(:),1)
    subject to
        Z >= 0;
cvx_end

try
    [recoveredSig_raw, eVal] = eig(Z);
catch
    [recoveredSig_raw, eVal] = eigs(Z);
end
recoveredSig = sqrt(eVal(end))*recoveredSig_raw(:,end);

