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
% This function helps to synchronize the frame to the correct beam pattern 
% used. It is used in the Program ID 2/4/6/ where the Tx phased array
% configuration is not at the Rx side.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% prmQPSKReceiver: a structure array that groups the related 
% parameters on the QPSK receiver including frame format, modulation 
% scheme and etc. With the generated structure array, the baseband post 
% signal processing is set.
% prmPAControl: a structure array that groups the related parameters
% on phased arrays including available codebooks, receiving power gain,
% antenna switching on/off settings and etc. With prmPAControl, the RF 
% configuration is set.
% BERMER_Non_Count: a structure array that groups the related evaluation 
% metrics which includes 'BER' (1st-row), 'MER in dB' (4-th row), 
% '95% MER in dB' (5-th row), 'Received signal strength indication (RSSI)' 
% (7-th row). BERMER_Non_Count does not perform the average operation, 
% which allows BER and MER analysis for each received frame separately.
% Raw_corruptSignal: all the collected raw corrupted signal.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output arguments:
% BERMER_Non_Count_Sync: synchronized BERMER_Non_Count
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function BERMER_Non_Count_Sync = Result_Frame_Synchronization(prmQPSKReceiver, prmPAControl, BERMER_Non_Count, Raw_corruptSignal) 
    %% Take the part where there is value
    BERMER_Non_Count = BERMER_Non_Count(:,:,1:prmPAControl.Number_State_To_Test*prmPAControl.Number_of_Frames_Per_Test_Rx_Buffered);
    Raw_corruptSignal_Sync = Raw_corruptSignal;
    Raw_corruptSignal = Raw_corruptSignal(:,:,1:prmPAControl.Number_State_To_Test);
    
    %% Reshape the BERMER_Non_Count
    Total_Frame = prmPAControl.Number_of_Frames_Per_Test_Tx*prmPAControl.Number_State_To_Test;
    BERMER_Non_Count_Reshape = reshape(BERMER_Non_Count,7,Total_Frame);
    index_phvec_state = BERMER_Non_Count_Reshape(7,:);   
    i = 2; j = 3;
    check = sum(BERMER_Non_Count(:,i,j) - BERMER_Non_Count_Reshape(:,(j-1)*prmQPSKReceiver.RxBufferedFrames+i));
    if check ~=0 
        disp('Error in reshaping in BERMER_Non_Count');
    end

    %% Find Boundary of different types of test states
    index_phvec_state_back = [index_phvec_state(2:end) 0];
    index_phvec_state_front = [0 index_phvec_state(1:end-1)];
    diff1 = index_phvec_state-index_phvec_state_back;
    diff2 = index_phvec_state-index_phvec_state_front;
    Set_index_repeat_begin = find(diff1==0);
    Set_index_change_begin = find(diff2==1);
    Potential_Boundary = intersect(Set_index_repeat_begin,Set_index_change_begin);
    Potential_Boundary_Copy = Potential_Boundary;
    for t = 1:1:length(Potential_Boundary)
        i = Potential_Boundary(t);
        index_phvec_state_front = index_phvec_state(i-3:i-1);
        if i+3 <= length(index_phvec_state)
            index_phvec_state_back = index_phvec_state(i:i+3);
        else
            index_phvec_state_back = index_phvec_state(i:end);
        end
        if mean(index_phvec_state_front) ~= index_phvec_state(i-1) || ...
                mean(index_phvec_state_back) ~= index_phvec_state(i)
            Potential_Boundary_Copy(t)=0;
        end
    end
    index_Boundary = find(Potential_Boundary_Copy~=0);
    index_Frame_Begin = Potential_Boundary_Copy(index_Boundary);
    Shift_Number_Selected = mode(mod(index_Frame_Begin,prmPAControl.Number_of_Frames_Per_Test_Tx)-1);
    if isnan(Shift_Number_Selected)
        Shift_Number_Selected = 0;
        disp('Error in synchronization');
    end
    
    %% Verify the boundary and shift the original BERMER_Non_Count data set
    index_phvec_state_shift = circshift(index_phvec_state,-Shift_Number_Selected);
    Reshape_Presentation = reshape(index_phvec_state_shift,prmPAControl.Number_of_Frames_Per_Test_Tx,prmPAControl.Number_State_To_Test);
    Reshape_Presentation_Temp = Reshape_Presentation;
    Check_Synchronization_Performance1 = mode(Reshape_Presentation_Temp);
    Reshape_Presentation_Temp(find(Reshape_Presentation_Temp==100)) = NaN;
    Check_Synchronization_Performance1 = round(nanmean(Reshape_Presentation_Temp));
    BERMER_Non_Count_Reshape_Shift = circshift(BERMER_Non_Count_Reshape,-Shift_Number_Selected,2);
    [r, c] = size(BERMER_Non_Count_Reshape_Shift);
    l = c/prmPAControl.Number_of_Frames_Per_Test_Tx;
    new_c = c/l;
    BERMER_Non_Count_New = reshape(BERMER_Non_Count_Reshape_Shift,r,new_c,l);
    BERMER_Non_Count_Index_Part = BERMER_Non_Count_New(7,:,:);
    Check_Synchronization_Performance2 = permute((mode(BERMER_Non_Count_Index_Part,2)),[1 3 2]);
    BERMER_Non_Count_Index_Part(find(BERMER_Non_Count_Index_Part==100)) = NaN;
    Check_Synchronization_Performance2 = permute(round(nanmean(BERMER_Non_Count_Index_Part,2)),[1 3 2]);
    Correct_Order = 1:prmPAControl.Number_State_To_Test;
    Common_Number = zeros(1,prmPAControl.Number_State_To_Test);
    for shift = 1:1:prmPAControl.Number_State_To_Test
        Diff = Correct_Order - circshift(Check_Synchronization_Performance2,shift);
        Common_Number(shift) = length(find(Diff==0));
    end
    Common_Number;
    [~, Big_Shift] = max(Common_Number);
    Big_Shift = Big_Shift;
    if isnan(Big_Shift)
        Big_Shift = 0;
        disp('Error in synchronization');
    end
    BERMER_Non_Count_Final = circshift(BERMER_Non_Count_New,Big_Shift,3);
    save BERMER_Non_Count_Final BERMER_Non_Count_Final
    BERMER_Non_Count_Sync = ... 
        reshape(BERMER_Non_Count_Final, 7, ...
                prmQPSKReceiver.RxBufferedFrames, ...
                prmPAControl.Number_of_Frames_Per_Test_Rx_Buffered*prmPAControl.Number_State_To_Test);

            
    %% Shift the original Raw_corruptSignal data set
    SamplesPerFrame = prmQPSKReceiver.FrameSize*prmQPSKReceiver.Upsampling*prmQPSKReceiver.RxBufferedFrames;
    Total_Frame = prmPAControl.Number_of_Frames_Per_Test_Rx_Buffered*prmPAControl.Number_State_To_Test;
    Raw_corruptSignal_Reshape = reshape(Raw_corruptSignal,SamplesPerFrame,Total_Frame);
    i = 2; j = 3;
    Check = sum(Raw_corruptSignal(:,i,j)...
        -Raw_corruptSignal_Reshape(:,(j-1)*prmPAControl.Number_of_Frames_Per_Test_Rx_Buffered+i));
    if Check ~=0 
        disp('Error in reshaping in Raw_corruptSignal');
    end
    Shift_Number_Selected_Raw_corruptSignal = round(Shift_Number_Selected/prmQPSKReceiver.RxBufferedFrames);
    Raw_corruptSignal_Reshape_Shift = circshift(Raw_corruptSignal_Reshape,-Shift_Number_Selected_Raw_corruptSignal,2);
    [r, c] = size(Raw_corruptSignal_Reshape_Shift);
    l = c/prmPAControl.Number_of_Frames_Per_Test_Rx_Buffered;
    new_c = c/l;
    Raw_corruptSignal_New = reshape(Raw_corruptSignal_Reshape_Shift,r,new_c,l);
    Raw_corruptSignal_Sync_Short = circshift(Raw_corruptSignal_New,Big_Shift,3);
    Raw_corruptSignal_Sync(:,:,1:prmPAControl.Number_State_To_Test) = Raw_corruptSignal_Sync_Short;
    save Raw_corruptSignal_Sync Raw_corruptSignal_Sync;
end


