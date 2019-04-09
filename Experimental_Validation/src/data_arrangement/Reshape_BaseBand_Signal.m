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
% This function reshapes the baseband signal and allows visualizing the
% demodulated signal. This function is used to test and analyze the 
% coherence of the received signal.
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
% Baseband_y: demodulated baseband signal.
% show_base_band_y: whether to show the demodulated baseband signal.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output arguments:
% Baseband_y_Final: reshaped demodulated basedband signal
% RSSI: RSSI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Baseband_y_Final, RSSI] = Reshape_BaseBand_Signal(prmQPSKReceiver, prmPAControl, Baseband_y, show_base_band_y) 
    %% Take the part where there is value
    Baseband_y = Baseband_y(:,:,1:prmPAControl.Number_State_To_Test*prmPAControl.Number_of_Frames_Per_Test_Rx_Buffered);
    
    %% Reshape the BERMER_Non_Count
    Total_Frame = prmPAControl.Number_of_Frames_Per_Test_Tx*prmPAControl.Number_State_To_Test;
    Baseband_y_Reshape = reshape(Baseband_y,prmQPSKReceiver.FrameSize,Total_Frame); 
    i = 2; j = 3;
    check = sum(Baseband_y(:,i,j) - Baseband_y_Reshape(:,(j-1)*prmQPSKReceiver.RxBufferedFrames+i));
    if check ~=0 
        disp('Error in reshaping in Baseband_y');
    end
    Baseband_y_Final = reshape(Baseband_y_Reshape,prmQPSKReceiver.FrameSize,prmPAControl.Number_of_Frames_Per_Test_Tx,prmPAControl.Number_State_To_Test);  
    
    i = 2; j = 3;
    check = sum(Baseband_y_Final(:,i,j) - Baseband_y_Reshape(:,(j-1)*prmPAControl.Number_of_Frames_Per_Test_Tx+i));
    if check ~=0 
        disp('Error in reshaping in Baseband_y');
    end
    
    %% Calculate Mean
    RSSI = zeros(prmPAControl.Number_State_To_Test,1);
    Ave_y = zeros(prmQPSKReceiver.FrameSize,prmPAControl.Number_State_To_Test);
    Beg = 1;
    Ave_y = Ave_y(Beg:end,:);
    
    %% Average performance
    if show_base_band_y
        figure
    end
    for i=1:1:prmPAControl.Number_State_To_Test
        Baseband_y_temp = Baseband_y_Final(Beg:end,:,i);
        sum_packet = sum(abs(Baseband_y_temp));
        pos = find(sum_packet == 0);
        Baseband_y_temp(:,pos) = [];
        sum_packet = sum(abs(Baseband_y_temp));
        pos = find(sum_packet >= 1000);
        Baseband_y_temp(:,pos) = [];
        Effective_Packet = length(Baseband_y_temp)
        RSSI(i) = rms(rms(Baseband_y_temp));        
        Ave_y(:,i) = mean(mean(Baseband_y_temp,3),2);
        if show_base_band_y
            subplot(prmPAControl.Number_State_To_Test,4,(i-1)*4+1)
            bar((real(Ave_y(Beg:end,i))));
            title('Real')
            subplot(prmPAControl.Number_State_To_Test,4,(i-1)*4+2)
            bar((imag(Ave_y(Beg:end,i))));
            title('Imag')
            subplot(prmPAControl.Number_State_To_Test,4,(i-1)*4+3)
            bar((rad2deg(angle(Ave_y(Beg:end,i)))));
            title('Angle')
            subplot(prmPAControl.Number_State_To_Test,4,(i-1)*4+4)
            bar(abs(abs(Ave_y(Beg:end,i))));
            title('Module')
        end
    end
    if show_base_band_y
        figure
        bar(RSSI)
        title('Average of All Packets')
    end
    
    %% Random frame
    RSSI_Random = zeros(prmPAControl.Number_State_To_Test,1);
    Ave_y_Random = zeros(prmQPSKReceiver.FrameSize,prmPAControl.Number_State_To_Test);
    Beg = 1;
    Ave_y_Random = Ave_y_Random(Beg:end,:);
    if show_base_band_y
        figure
    end
    for i=1:1:prmPAControl.Number_State_To_Test
        Baseband_y_temp = Baseband_y_Final(Beg:end,:,i);
        sum_packet = sum(abs(Baseband_y_temp));
        pos = find(sum_packet == 0);
        Baseband_y_temp(:,pos) = [];
        sum_packet = sum(abs(Baseband_y_temp));
        pos = find(sum_packet >= 1000);
        Baseband_y_temp(:,pos) = [];
        Random_Frame = randi([1 length(Baseband_y_temp(1,:))],1,2);
        Baseband_y_temp = Baseband_y_temp(Beg:end,Random_Frame);
        RSSI_Random(i) = rms(rms(Baseband_y_temp));  
        Ave_y_Random(1:length(Baseband_y_temp),i) = mean(mean(Baseband_y_temp,3),2);
        if show_base_band_y
            subplot(prmPAControl.Number_State_To_Test,4,(i-1)*4+1)
            bar((real(Ave_y_Random(Beg:end,i))));
            title('Real')
            subplot(prmPAControl.Number_State_To_Test,4,(i-1)*4+2)
            bar((imag(Ave_y_Random(Beg:end,i))));
            title('Imag')
            subplot(prmPAControl.Number_State_To_Test,4,(i-1)*4+3)
            bar((rad2deg(angle(Ave_y_Random(Beg:end,i)))));
            title('Phase')
            subplot(prmPAControl.Number_State_To_Test,4,(i-1)*4+4)
            bar(abs(abs(Ave_y_Random(Beg:end,i))));
            title('Module')
        end
    end
    if show_base_band_y
        figure
        bar(RSSI_Random)
        title('Random Packet')
    end
    
    %% Original Frame
    if show_base_band_y
        figure
        load('Baseband_y_original.mat');
        subplot(4,1,1)
        bar((real(Baseband_y_original(Beg:end,i))));
        title('Real')
        subplot(4,1,2)
        bar((imag(Baseband_y_original(Beg:end,i))));
        title('Imag')
        subplot(4,1,3)
        bar((rad2deg(angle(Baseband_y_original(Beg:end,i)))));
        title('Phase')
        subplot(4,1,4)
        bar(abs(abs(Baseband_y_original(Beg:end,i))));
        title('Module')
    end
end


