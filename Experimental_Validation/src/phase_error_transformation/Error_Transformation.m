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
% This function is called by phase_error_calculation.m to calculate phase
% error. In particualr, it removes the trigonometric ambiguity. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [error_array, L] = Error_Transformation(sol, error_min, error_max, Remove_Ambiguity)
    error = asin(sol);
    error_array1 = error + 2*pi*(-3:1:3);
    error_array1 = error_array1(find(error_array1 >= error_min & error_array1 <= error_max));
    error_array2 = (pi-error) + 2*pi*(-3:1:3);
    error_array2 = error_array2(find(error_array2 >= error_min & error_array2 <= error_max));
    error_array = [error_array1 error_array2];
    L = length(error_array);
    
    %remove ambiguity of error
    for l=1:L
        error_temp = error_array(l);
        p = Remove_Ambiguity.CheckPoint(1);
        q = Remove_Ambiguity.CheckPoint(2);
        sign_Ref = Remove_Ambiguity.mp-Remove_Ambiguity.mq;
        a_1_1 = Remove_Ambiguity.a_1_1;
        a_n_p = Remove_Ambiguity.a_n_p;
        a_n_q = Remove_Ambiguity.a_n_q;
        sign_To_Check = abs(a_1_1+a_n_p*exp(1j*(error_temp+(p-1)*pi/2))) - abs(a_1_1+a_n_q*exp(1j*(error_temp+(q-1)*pi/2)));
        if sign_Ref*sign_To_Check < 0
            error_array(l) = inf;
        end
    end
    error_array = error_array(find(error_array<inf));
    L = length(error_array);
end