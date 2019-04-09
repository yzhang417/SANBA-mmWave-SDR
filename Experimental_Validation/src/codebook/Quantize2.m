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
% This function compensates the phase error (with four methods) and then
% provides the codebooks for a desired beam pattern. It is a generalized
% version of [Experimental_Validation/src/codebook/Quantize.m]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% pattern: desired beam pattern.
% calibrate: coarse calibrate result.
% phase_error: fine calibration method 1.
% Fine_calibrate: fine calibration method 2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output arguments:
% codebook1: codebook for desired pattern without calibration.
% codebook2: codebook for desired pattern with coarse calibration.
% codebook3: codebook for desired pattern with fine calibration method 1.
% codebook4: codebook for desired pattern with fine calibration method 2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [codebook1, codebook2, codebook3, codebook4] = Quantize2(pattern, calibrate, phase_error, Fine_calibrate)
    %% Non calibration
    codebook1 = pattern;
    len = length(pattern);
    phi = [0, pi/2, pi, -pi/2, -pi];
    quant = [0 1 2 3 2];
    for i=1:1:len
       [v, ind] = min( abs(angle(pattern(i)) - phi) );
       codebook1(i) = quant(ind);
    end

    %% Coarse calibration
    codebook2 = codebook1 + calibrate';
    codebook2 = mod(codebook2,4);
    
    %% Fine calibration method 1
    codebook3 = pattern;
    len = length(pattern);
    for i=1:1:len
       phi = [0, pi/2, pi, 3*pi/2] + phase_error(i);  
       phi = wrapToPi(phi);
       for j = 1:length(phi)
           while phi(j) < 0 
               phi(j) = phi(j) + 2*pi;
           end
       end
       quant = [0 1 2 3];   
       target1 = angle(pattern(i));
       while target1<0
           target1 = target1 + 2*pi;
       end
       target2 = target1 + 2*pi;
       [v1, ind1] = min( abs(target1 - phi) );
       [v2, ind2] = min( abs(target2 - phi) );
       if v1<v2
           codebook3(i) = quant(ind1);
       else
           codebook3(i) = quant(ind2);
       end
    end
      
    %% Fine calibration method 2
    codebook4 = pattern;
    len = length(pattern);
    for i=1:1:len
       phi = Fine_calibrate(:,i)';  
       phi = wrapToPi(phi);
       for j = 1:length(phi)
           while phi(j) < 0 
               phi(j) = phi(j) + 2*pi;
           end
       end
       quant = [0 1 2 3];   
       target1 = angle(pattern(i));
       while target1<0
           target1 = target1 + 2*pi;
       end
       target2 = target1 + 2*pi;
       [v1, ind1] = min( abs(target1 - phi) );
       [v2, ind2] = min( abs(target2 - phi) );
       if v1<v2
           codebook4(i) = quant(ind1);
       else
           codebook4(i) = quant(ind2);
       end
    end    
end

