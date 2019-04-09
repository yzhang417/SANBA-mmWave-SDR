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
% This function removes the boundary frames to ensure that the frames to be
% processed are well transmitted by the expected beam pattern
% configuration.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% prmPAControl: a structure array that groups the related parameters
% on phased arrays including available codebooks, receiving power gain,
% antenna switching on/off settings and etc. With prmPAControl, the RF 
% configuration is set.
% BERMER_Non_Count: a structure array that groups the related evaluation 
% metrics which includes 'BER' (1st-row), 'MER in dB' (4-th row), 
% '95% MER in dB' (5-th row), 'Received signal strength indication (RSSI)' 
% (7-th row). BERMER_Non_Count does not perform the average operation, 
% which allows BER and MER analysis for each received frame separately.
% Raw_corruptSignal: raw signal collected.
% Remove_Propotional: proportion of the frames that are
% removed to ensure that the counted frame is correctly transmitted by the
% intended array configuration.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output arguments:
% BERMER_Removing: a structure array that groups the related evaluation 
% metrics which includes 'BER' (1st-row), 'MER in dB' (4-th row), 
% '95% MER in dB' (5-th row), 'Received signal strength indication (RSSI)' 
% (7-th row). The boundary frames are removed and the average operation is
% performed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function BERMER_Removing = Removing_Boundary_Frames(prmPAControl, BERMER_Non_Count, Raw_corruptSignal, Remove_Propotional)    
    %% Take the part where there is value
    BERMER_Non_Count = BERMER_Non_Count(:,:,1:prmPAControl.Number_State_To_Test*prmPAControl.Number_of_Frames_Per_Test_Rx_Buffered);
    Raw_corruptSignal = Raw_corruptSignal(:,:,1:prmPAControl.Number_State_To_Test);    

    %% Calculate the power
    %NBF_Remove: Number_of_Boundary_Frame_To_Remocw_For_Power_Calculation
    NBF_Remove = prmPAControl.Number_of_Frames_Per_Test_Rx_Buffered * Remove_Propotional;
    Begininig_Frame = 1 + NBF_Remove;
    Ending_Frame = prmPAControl.Number_of_Frames_Per_Test_Rx_Buffered - NBF_Remove;    
    BERMER_Non_Count_Power = zeros(1,prmPAControl.Number_State_To_Test);
    for i =1:1:prmPAControl.Number_State_To_Test
        BERMER_Non_Count_Power(i) = rms(rms(Raw_corruptSignal(:,Begininig_Frame:Ending_Frame,i)));
    end
    
    %% Average of BERMER_Non_Count Matrix
    Beginning_Frame = 1 + Remove_Propotional * prmPAControl.Number_of_Frames_Per_Test_Tx;
    Ending_Frame = prmPAControl.Number_of_Frames_Per_Test_Tx - ...
         Remove_Propotional*prmPAControl.Number_of_Frames_Per_Test_Tx;
    Total_Frame = prmPAControl.Number_of_Frames_Per_Test_Tx*prmPAControl.Number_State_To_Test;
    BERMER_Non_Count_Reshape = reshape(BERMER_Non_Count,7,Total_Frame);
    [r, c] = size(BERMER_Non_Count_Reshape);
    l = c/prmPAControl.Number_of_Frames_Per_Test_Tx;
    new_c = c/l;
    BERMER_Non_Count_After_Reshape = reshape(BERMER_Non_Count_Reshape,r,new_c,l);
    BERMER_Non_Count_After_Reshape_Removing = BERMER_Non_Count_After_Reshape(:,:,:);
    BERMER_Removing = zeros(7,prmPAControl.Number_State_To_Test);
    for i=1:1:prmPAControl.Number_State_To_Test
        BERMER_Non_Count_After_Reshape_Removing_i = BERMER_Non_Count_After_Reshape_Removing(:,Beginning_Frame:Ending_Frame,i);
        BERMER_Removing(1,i) = mean(BERMER_Non_Count_After_Reshape_Removing_i(1,:));
        BERMER_Removing(2,i) = sum(BERMER_Non_Count_After_Reshape_Removing_i(2,:));
        BERMER_Removing(3,i) = sum(BERMER_Non_Count_After_Reshape_Removing_i(3,:));
        BERMER_Removing(4,i) = 10*log10(mean(10.^(BERMER_Non_Count_After_Reshape_Removing_i(4,:)./10)));
        %BERMER_Removing(4,i) = mean(BERMER_Non_Count_After_Reshape_Removing_i(4,:));
        BERMER_Removing(5,i) = min(BERMER_Non_Count_After_Reshape_Removing_i(5,:));
        BERMER_Removing(6,i) = prctile(BERMER_Non_Count_After_Reshape_Removing_i(6,:),0.95);
        BERMER_Removing(7,i) = BERMER_Non_Count_Power(i);
    end    
end