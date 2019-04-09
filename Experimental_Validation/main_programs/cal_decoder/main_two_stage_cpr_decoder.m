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
% This script is specifically for checking details of the performance of 
% the proposed non-coherent beam alignment algorithm. It is a simplified
% version of [Experimental_Validation/src/general_receiver/
% /Decoding_ResultVisulization.m].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization and load of data
global Raw_corruptSignal
load('prmQPSKReceiver.mat');
load('prmPAControl.mat');

%%
r = 7;
c = prmPAControl.Number_State_To_Test;
Number_Test_Beg = 1;
Number_Test_End = 1;
TotalBERMER = zeros(r,c,Number_Test_End-Number_Test_Beg+1);

%% Decoding and visuliazation
TotalBER_MER_Power_Overall_Training = zeros(r,c,prmPAControl.nthTxBeam,Number_Test_End-Number_Test_Beg+1);
for nthTest = Number_Test_Beg : 1 : Number_Test_End
    if prmPAControl.Program_ID == 6
        NTRxBeams = prmPAControl.nthRxBeam;
    else
        NTRxBeams = prmPAControl.nthTxBeam;
    end
    for nthTRxBeam = 1:1:NTRxBeams
        load(['Raw_corruptSignal' num2str(nthTest) '_' num2str(nthTRxBeam) '.mat'],'Raw_corruptSignal');
        clear run_offline_decoder_mex %#ok<UNRCH>
        BERMER = run_offline_decoder_mex(prmQPSKReceiver, prmPAControl, 0);
        save(['BERMER' num2str(nthTest) '_' num2str(nthTRxBeam) '.mat'],'BERMER');
        TotalBER_MER_Power_Overall_Training(:,:,nthTRxBeam,nthTest) = BERMER;
    end
end
BER_MER_Power_Overall_Training = mean(TotalBER_MER_Power_Overall_Training,4);
save(['BER_MER_Power_Overall_Training' num2str(nthTest) '.mat'],'BER_MER_Power_Overall_Training');



