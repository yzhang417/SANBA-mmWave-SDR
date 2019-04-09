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
% This function permits calculating and plotting the recovery performance
% of the phase retrieval algorithm.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% TrueSig: true target sparse vector (column vector).
% recoveredSig: recovered sparse vector (column vector).
% Method_Name: name of the method (algorithm) used (character array).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Plot_Recovery_Performance(TrueSig,recoveredSig,Method_Name)
    % Compute error 
    err_2 = norm( TrueSig - recoveredSig )/ norm(TrueSig);
    err_dB = 10*log10( norm(TrueSig-recoveredSig,2)^2 / norm(TrueSig,2)^2 );

    % Print errors to standard output
    fprintf( '\nError in recovery (2-norm) is %3.3e \n', err_2 );
    fprintf( 'Error in recovery (dB) is %3.3e \n', err_dB );

    % Plot reconstructions
    % ... real part
    figure; subplot(1,4,1); 
    stem( real(TrueSig), 'r-+', 'linewidth', 2 ); hold on
    stem( real(recoveredSig), 'b-o', 'linewidth', 2 );
    xlabel 'k'; ylabel 'Re(x[k])'; grid; title 'Real Part'
    legend( 'True', 'Recovered' )

    % ... imaginary part
    subplot(1,4,2); 
    stem( imag(TrueSig), 'r-+', 'linewidth', 2 ); hold on
    stem( imag(recoveredSig), 'b-o', 'linewidth', 2 );
    xlabel 'k'; ylabel 'Im(x[k])'; grid; title 'Imaginary Part'
    legend( 'True', 'Recovered' )

    % ... Absolute value
    subplot(1,4,3); 
    stem( abs(TrueSig), 'r-+', 'linewidth', 2 ); hold on
    stem( abs(recoveredSig), 'b-o', 'linewidth', 2 );
    xlabel 'k'; ylabel 'Abs(x[k])'; grid; title 'Absolute Value'
    legend( 'True', 'Recovered' )

    % ... Angle value
    subplot(1,4,4); 
    stem( angle(TrueSig), 'r-+', 'linewidth', 2 ); hold on
    stem( angle(recoveredSig), 'b-o', 'linewidth', 2 );
    xlabel 'k'; ylabel 'angle(x[k])'; grid; title 'Angle Value'
    legend( 'True', 'Recovered' )
    
    % Add title
    title(Method_Name);
end