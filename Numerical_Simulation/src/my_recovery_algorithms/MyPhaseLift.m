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
% This function uses PhaseLift via TFOCS optimization package
% to solve the compressive phase retrieval problem. For details, see
%   Emmanuel J.Candès, Yonina C.Eldar, Thomas Strohmer, and Vladislav
%   Voroninski. 2013. Phase Retrieval via Matrix Completion. SIAM Journal 
%   on Imaging Sciences 6, 1 (2013), 199-225.
%   and 
%   Emmanuel J. Candès, Thomas Strohmer, and Vladislav Voroninski. 2013. 
%   PhaseLift: Exact and Stable Signal Recovery from Magnitude Measurements 
%   via Convex Programming. Communications on Pure and Applied Mathematics 
%   66, 8 (2013), 1241-1274.
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


function recoveredSig = MyPhaseLift(measurements, measurementMat)
%% Problem parameters
[m, n] = size(measurementMat);
A = measurementMat;
y = measurements;

%% PhaseLift
%   Solve the Phase Retrieval problem using m measurements, using PhaseLift
%
%   We solve
%       min      0.5 || b - A(X) ||_2^2 + lambda*trace(X)
%   subject to          X >= 0

% Set TFOCS parameters
opts.maxIts = 4e3;
opts.tol = 1e-10;
opts.restart = 200;
opts.printEvery = 0;
if n > 2000
    opts.largescale = true;
end

% Regularization parameter
lambda = 5e-2;

% Initial guess
initGuess = zeros(n);

% TFOCS requires special initialization of the measurement model
tfocsLinOp = initializeLinopPR( {[n, n], m}, measurementMat ); 

recoveredMat = solver_TraceLS( tfocsLinOp, measurements, ...
                                            lambda, initGuess, opts );

% The above SDP recovers the matrix X = xx*; we need to extract x
% Since X is low-rank (ideally rank-1), choose solution to be (scaled) 
% leading eigenvector                                        
[recoveredSig, eVal] = eig(recoveredMat);
recoveredSig = sqrt(eVal(end))*recoveredSig(:,end);

