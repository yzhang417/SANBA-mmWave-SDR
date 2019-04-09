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
% This function uses the GAMP to solve the compressive phase retrieval 
% problem directly. The performance is not good. For details of the GAMP 
% package, please visit https://sourceforge.net/projects/gampmatlab/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% measurements: received noised measurements (column vector).
% measurementMat: sensing matrix.
% opt_pr: parameters required for prGAMP4 algorithm in GAMP package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output arguments:
% recoveredSig: recovered sparse signal (column vector).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function recoveredSig = MyPRGAMP(measurements, measurementMat, opt_pr)
%% Problem parameters
[m, n] = size(measurementMat);
A = measurementMat;
y = measurements;

%% PRGAMP
try
    recoveredSig = prGAMP4(sqrt(y), A, opt_pr);
catch
    fprintf('\nFail of PRGAMP in MyPRGAMP\n');
    recoveredSig = zeros(n,1);
end

