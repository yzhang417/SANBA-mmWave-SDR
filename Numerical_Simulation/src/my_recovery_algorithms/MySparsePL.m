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
% This function uses the idea from the following paper to solve 
% the compressive phase retrieval problem.
% For details, see
%   P. Netrapalli, P. Jain and S. Sanghavi, "Phase Retrieval Using 
%   Alternating Minimization," in IEEE Transactions on Signal Processing, 
%   vol. 63, no. 18, pp. 4814-4826, Sept.15, 2015.
% The algorithm proposed in the above paper eliminates the components, 
% of the target sparse vector, which is to be zero with high probability. 
% The performance is not good.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% measurements: received noised measurements (column vector).
% measurementMat: sensing matrix.
% opt_AltMin: parameters required for the alternative minimization
% algorithm.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output arguments:
% recoveredSig: recovered sparse signal (column vector).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function recoveredSig = MySparsePL(measurements, measurementMat, opt_AltMin)
%% Problem parameters
[m, n] = size(measurementMat);
A = measurementMat;
y = sqrt(measurements);

%% SparseAltMinPhase

% Elimination of certain columns
k = ceil(0.05*n);
rep_y = repmat(y,1,n);
Corr = sum(abs(A.*rep_y));
[Corr_Sorted, ind] = sort(Corr,'descend');
ind = ind(1:k);
A = A(:,ind);

% % Partition of A
% c = 10;
% t0 = c*log(1/opt_AltMin.opt_epslion);
% sub_m = max(1,floor(m/t0));
% t0 = floor(m/sub_m);
% rest = m - t0*sub_m;
% if rest > 0.5 * sub_m
%     t0 = t0 + 1;
% end
% 
% % Initialize the initial point
% S = 0;
% for i=1:sub_m
%     S = S + y(i).^2 + A(i,:)'*A(i,:);
% end
% [initial_Guess, eVal] = eig(S);
% initial_Guess = sqrt(eVal(end))*initial_Guess(:,end);
% xt = initial_Guess;
% 
% % Iteration
% for t=1:1:t0-1
%     range = t*sub_m+1:(t+1)*sub_m;
%     if t == t0-1
%         range = (t*sub_m+1):m;
%     end
%     Az = A(range,:)*xt;
%     C = diag(Az/norm(Az));
%     xt = pinv(A(range,:))*C*y(range,:);
% end

% Output result
recoveredSig = zeros(n,1);
xt = MyPhaseLift(y.^2, A);
recoveredSig(ind) = xt;


