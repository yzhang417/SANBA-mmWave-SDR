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
% This main script runs the decoder to decode the corrupted signal and 
% visualize the result.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization and load of data
global Raw_corruptSignal
global BERMER_Non_Count
load('prmQPSKReceiver.mat');
load('prmPAControl.mat');
r = 7;
c = prmPAControl.Number_State_To_Test;
Number_Test_Beg = 1;
Number_Test_End = 1;
TotalBERMER = zeros(r,c,Number_Test_End-Number_Test_Beg+1);

%% Decoding and visualization
if prmPAControl.Program_ID ~= 6 && prmPAControl.Program_ID ~= 8.1
    for nthTest = Number_Test_Beg : 1 : Number_Test_End
        fprintf('\n');
        disp('---------------------------------------------')
        fprintf('Decoding %d-th Test\n', nthTest);
        disp('---------------------------------------------')
        load(['Raw_corruptSignal' num2str(nthTest) '.mat'],'Raw_corruptSignal');
        clear run_Offline_SDRuQPSKDecoder_mex %#ok<UNRCH>
        tic
        BERMER = run_Offline_SDRuQPSKDecoder_mex(prmQPSKReceiver, prmPAControl, 0);
        toc
        save(['BERMER' num2str(nthTest) '.mat'],'BERMER');
        save(['BERMER_Non_Count' num2str(nthTest) '.mat'],'BERMER_Non_Count');
        save(['Baseband_y' num2str(nthTest) '_' num2str(nthTest) '.mat'],'Baseband_y');

        %% Post-processing and visulization    
        switch prmPAControl.Program_ID
            case {0, 1, 3, 5, 2.1, 4.1, 8, 9, 9.1, 9.2}
                %Result without removing boundary frames
                Result_Visulization(prmPAControl, BERMER);

%                 %Result with removing boundary frames by averaging BER_Non_Count
%                 BERMER_Removing = Removing_Boundary_Frames(prmPAControl, BERMER_Non_Count, Raw_corruptSignal, 0.2);
%                 Result_Visulization(prmPAControl, BERMER_Removing);
% 
%                 %Result with removing boundary frames by redecoding
%                 clear run_Offline_SDRuQPSKDecoder_mex %#ok<UNRCH>
%                 BERMER_Redecoding_Removing = run_Offline_SDRuQPSKDecoder_mex(prmQPSKReceiver, prmPAControl, 0.3);
%                 Result_Visulization(prmPAControl, BERMER_Redecoding_Removing);

            case {2, 4}
                %Result_Visulization(prmPAControl, BERMER);
                load(['Raw_corruptSignal' num2str(nthTest) '.mat'],'Raw_corruptSignal');
                load(['BERMER' num2str(nthTest) '.mat'],'BERMER');
                BERMER_Non_Count_Sync = Result_Frame_Synchronization2(prmQPSKReceiver, prmPAControl, BERMER_Non_Count, Raw_corruptSignal);
                load(['Raw_corruptSignal_Sync.mat'],'Raw_corruptSignal_Sync');
                %Result without removing boundary frames
                Remove_Propotional1 = 0;
                BERMER_Sync_Removing1 = Removing_Boundary_Frames(prmPAControl, BERMER_Non_Count_Sync, Raw_corruptSignal_Sync, Remove_Propotional1);
                Result_Visulization(prmPAControl, BERMER_Sync_Removing1);

%                 %Result with removing boundary frames
%                 Remove_Propotional2 = 0.2;
%                 BERMER_Sync_Removing2 = Removing_Boundary_Frames(prmPAControl, BERMER_Non_Count_Sync, Raw_corruptSignal_Sync, Remove_Propotional2);
%                 Result_Visulization(prmPAControl, BERMER_Sync_Removing2);
% 
%                 %Result with removing boundary frames
%                 Remove_Propotional3 = 0.4;
%                 BERMER_Sync_Removing3 = Removing_Boundary_Frames(prmPAControl, BERMER_Non_Count_Sync, Raw_corruptSignal_Sync, Remove_Propotional3);
%                 Result_Visulization(prmPAControl, BERMER_Sync_Removing3);   

%                 %Load sync-data and reprocess it: Result with removing boundary frames
%                 load(['Raw_corruptSignal_Sync.mat'],'Raw_corruptSignal_Sync');
%                 Raw_corruptSignal = Raw_corruptSignal_Sync;
%                 tic
%                 clear run_Offline_SDRuQPSKDecoder_mex %#ok<UNRCH>
%                 BERMER_Redecoded_After_Sync_Removing = run_Offline_SDRuQPSKDecoder_mex(prmQPSKReceiver, prmPAControl, 0.3);
%                 toc
%                 Result_Visulization(prmPAControl, BERMER_Redecoded_After_Sync_Removing);
        end
    end
else
    TotalBER_MER_Power_Overall_Training = zeros(r,c,prmPAControl.nthTxBeam,Number_Test_End-Number_Test_Beg+1);
    for nthTest = Number_Test_Beg : 1 : Number_Test_End
        fprintf('\n');
        disp('---------------------------------------------')
        fprintf('Decoding %d-th Test\n', nthTest);
        disp('---------------------------------------------')
        if prmPAControl.Program_ID == 6
            NTRxBeams = prmPAControl.nthRxBeam;
        else
            NTRxBeams = prmPAControl.nthTxBeam;
        end
        for nthTRxBeam = 1:1:NTRxBeams
            disp('---------------------------------------------------------------------')
            if prmPAControl.Program_ID == 8.1 
                fprintf('Decoding best receiver beam pattern for %d-th transmiter beam pattern\n', nthTRxBeam);
            else
                fprintf('Decoding best transmitter beam pattern for %d-th receiver beam pattern\n', nthTRxBeam);
            end
            disp('---------------------------------------------------------------------')
            load(['Raw_corruptSignal' num2str(nthTest) '_' num2str(nthTRxBeam) '.mat'],'Raw_corruptSignal');
            tic
            clear run_Offline_SDRuQPSKDecoder_mex %#ok<UNRCH>
            BERMER = run_Offline_SDRuQPSKDecoder_mex(prmQPSKReceiver, prmPAControl, 0);
            toc
            save(['BERMER' num2str(nthTest) '_' num2str(nthTRxBeam) '.mat'],'BERMER');
            save(['BERMER_Non_Count' num2str(nthTest) '_' num2str(nthTRxBeam)  '.mat'],'BERMER_Non_Count');
            save(['Baseband_y' num2str(nthTest) '_' num2str(nthTRxBeam) '.mat'],'Baseband_y');
            TotalBER_MER_Power_Overall_Training(:,:,nthTRxBeam,nthTest) = BERMER;
            
            %% Post-processing and visualizzation    
            if prmPAControl.Program_ID == 6
                load(['Raw_corruptSignal' num2str(nthTest) '_' num2str(nthTRxBeam) '.mat'],'Raw_corruptSignal');
                load(['BERMER' num2str(nthTest) '_' num2str(nthTRxBeam) '.mat'],'BERMER');
                BERMER_Non_Count_Sync = Result_Frame_Synchronization2(prmQPSKReceiver, prmPAControl, BERMER_Non_Count, Raw_corruptSignal);

                %Load sync-data and reprocess it: Result with removing boundary frames
                load(['Raw_corruptSignal_Sync.mat'],'Raw_corruptSignal_Sync');
                Raw_corruptSignal = Raw_corruptSignal_Sync;
                tic
                clear run_Offline_SDRuQPSKDecoder_mex %#ok<UNRCH>
                BERMER_Redecoded_After_Sync_Removing = run_Offline_SDRuQPSKDecoder_mex(prmQPSKReceiver, prmPAControl, 0.30);
                toc
                TotalBER_MER_Power_Overall_Training(:,:,nthTRxBeam,nthTest) = BERMER_Redecoded_After_Sync_Removing(:,1:11);

%                 %Load sync-data and reprocess it to verify the data sync
%                 load(['Raw_corruptSignal_Sync.mat'],'Raw_corruptSignal_Sync');
%                 Raw_corruptSignal = Raw_corruptSignal_Sync;
%                 clear run_Offline_SDRuQPSKDecoder_mex %#ok<UNRCH>
%                 BERMER_Redecoded_After_Sync = run_Offline_SDRuQPSKDecoder_mex(prmQPSKReceiver, prmPAControl, 0);
%                 BERMER_Non_Count_Re_Sync = Result_Frame_Synchronization2(prmQPSKReceiver, prmPAControl, BERMER_Non_Count, Raw_corruptSignal);
            end
        end
    end
    BER_MER_Power_Overall_Training = mean(TotalBER_MER_Power_Overall_Training,4);
    save(['BER_MER_Power_Overall_Training' num2str(nthTest) '.mat'],'BER_MER_Power_Overall_Training');
    Result_Visulization(prmPAControl, BER_MER_Power_Overall_Training);
end
    


