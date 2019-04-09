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
% Script description:
% This script numerically demonstrates the measured codebook. Also, the 
% UART commands corresponding to the codebook is printed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ULA Parameters
ULA.M = 12;
ULA.lambda = 3*10^8/(60.48*10^9);
ULA.d = 3.055*10^(-3); %It is measured from the board that the spacing is around 3mm

%% Phase_Array_1: R means in first quadrant
%ID: 1419-01
Phvec50R = [1 4 1 2 4 4 3 4 2 1 2 3]-1; %Ok
Phvec40R = [1 1 1 2 2 2 2 2 1 4 2 3]-1; %Ok
Phvec30R = [1 1 4 1 2 3 3 4 2 2 4 1]-1; %Ok
Phvec20R = [1 1 4 2 3 4 4 1 1 1 3 1]-1; %Ok
Phvec10R = [1 2 4 2 4 4 1 3 3 3 2 4]-1; %Ok
Phvec0   = [1 2 4 3 1 2 3 2 2 3 2 4]-1; %Ok
Phvec10L = [1 3 1 4 2 3 1 4 1 2 2 1]-1; %Ok
Phvec20L = [1 3 1 4 2 4 2 2 3 1 1 4]-1; %Ok
Phvec30L = [1 3 4 1 3 1 4 4 1 4 4 4]-1; %Ok
Phvec40L = [1 4 4 4 4 2 2 2 3 2 3 3]-1; %Ok
Phvec50L = [1 4 4 1 1 4 3 4 1 4 2 2]-1; %Ok

%
calibrate = Phvec0;
Codebook_Pattern = ...
[ Phvec50L' Phvec40L' Phvec30L' Phvec20L' Phvec10L' Phvec0' ...
  Phvec10R' Phvec20R' Phvec30R' Phvec40R' Phvec50R'];
Theta_disired_d = -50:10:50;

%
Comparing_Test_To_Theory_BeamPattern(Theta_disired_d, Codebook_Pattern, calibrate, ULA);

%
fprintf('Phase Array I CodeBook\n')
Phvec_ID_1 = ... 
[phvec(Phvec50R);...
 phvec(Phvec40R);...
 phvec(Phvec30R);...
 phvec(Phvec20R);...
 phvec(Phvec10R);...
 phvec(Phvec0);...
 phvec(Phvec10L);...
 phvec(Phvec20L);...
 phvec(Phvec30L);...
 phvec(Phvec40L);...
 phvec(Phvec50L)]
fprintf('\n\n\n')


%% Receiver: R means the angle in the first quadrant
%ID: 0900-02
Phvec50L = [1 4 3 4 1 4 4 4 2 1 1 2]-1; %OK
Phvec40L = [1 3 4 4 4 2 2 2 3 3 2 3]-1; %OK
Phvec30L = [1 3 4 4 3 1 1 4 1 4 4 4]-1; %OK
Phvec20L = [1 3 1 1 3 1 4 2 4 2 1 1]-1; %OK
Phvec10L = [1 2 4 4 2 3 1 4 4 2 1 1]-1; %OK
Phvec0   = [1 2 4 3 1 2 4 2 2 4 2 1]-1; %Ok
Phvec10R = [1 2 1 4 4 1 2 4 4 2 4 2]-1; %Ok 
Phvec20R = [1 1 4 2 4 4 1 2 1 2 4 3]-1; %Ok
Phvec30R = [1 1 4 2 3 2 4 4 3 4 1 2]-1; %Ok
Phvec40R = [1 1 1 2 2 2 2 2 1 1 2 3]-1; %Ok
Phvec50R = [1 4 4 1 1 4 4 4 2 2 2 3]-1; %Ok

%
calibrate = Phvec0;
Codebook_Pattern = ...
[ Phvec50L' Phvec40L' Phvec30L' Phvec20L' Phvec10L' Phvec0' ...
  Phvec10R' Phvec20R' Phvec30R' Phvec40R' Phvec50R'];
Theta_disired_d = -50:10:50;

%
Comparing_Test_To_Theory_BeamPattern(Theta_disired_d, Codebook_Pattern, calibrate, ULA);

fprintf('Phase Array I \n')
Phvec_ID_2 = ... 
[phvec(Phvec50R);...
 phvec(Phvec40R);...
 phvec(Phvec30R);...
 phvec(Phvec20R);...
 phvec(Phvec10R);...
 phvec(Phvec0);...
 phvec(Phvec10L);...
 phvec(Phvec20L);...
 phvec(Phvec30L);...
 phvec(Phvec40L);...
 phvec(Phvec50L)]
fprintf('\n\n\n')


%% Historical records
% Phvec50R = [1 4 3 4 1 4 3 4 1 4 1 2]-1; %Ok
% Phvec40R = [1 4 2 2 3 2 1 1 2 1 2 1]-1; %Ok
% Phvec30R = [1 3 3 3 3 1 4 3 4 3 3 2]-1; %Ok
% Phvec20R = [1 3 1 1 2 1 3 2 3 2 1 1]-1; %Ok
% Phvec10R = [1 2 4 2 4 1 2 4 3 4 2 4]-1; %OK

% Phvec50R = [1 4 2 3 4 3 2 2 4 3 4 4]-1; Phvec50R = [1 4 3 4 1 4 3 4 1 4 1 2]-1;
% Phvec40R = [1 4 3 2 3 2 1 1 2 1 2 2]-1; Phvec40R = [1 3 2 2 3 2 1 1 2 1 2 1]-1; Phvec40R = [1 4 2 2 3 2 1 1 2 1 2 1]-1;
% Phvec30R = [1 3 3 3 3 1 4 3 4 3 3 2]-1;
% Phvec20R = [1 3 1 1 2 1 3 2 3 2 1 1]-1;
% Phvec10R = [1 2 4 2 4 1 2 4 3 4 2 4]-1;
% Phvec0      = [0 1 3 2 0 1 2 1 1 2 1 3];
% Phvec10L   = [1 2 1 3 4 1 2 4 4 4 3 1]-1;
% Phvec20L   = [1 1 3 1 3 3 4 1 4 4 3 4]-1;
% Phvec30L   = [1 1 4 1 2 3 3 4 2 2 4 1]-1;
% Phvec40L   = [1 1 1 2 2 2 2 2 1 4 2 3]-1;  Phvec40L   = [1 1 4 1 2 2 1 2 4 4 1 2]-1;
% Phvec50L   = [1 4 1 2 1 4 4 4 2 1 3 3]-1;  %Phvec50L = [1 4 4 4 4 3 2 3 1 4 1 2]-1;  %Phvec50L  = [1 4 4 4 4 4 3 3 1 4 1 2]-1;

