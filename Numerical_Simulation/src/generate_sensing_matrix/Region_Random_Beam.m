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
% This function generates random precoder and combiner matrix within a 
% particular searching region. The generated beams are supposed to randomly
% cover the spatial domain.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% Mt: number of measurements at Tx side.
% Mr: number of measurements at Rx side.
% ULA: a structure array that groups related parameters on ULA.
% H: a structure array that groups related parameters on mmWave channel,
% please refer to Generate_Channel.m for details of its fields.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output arguments:
% F: precoder matrix (each column represents one beam pattern codebook).
% W: combiner matrix (each column represents one beam pattern codebook).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [F, W] = Region_Random_Beam(Mt,Mr,ULA,H)
    % Parameter fetching
    lambda = ULA.lambda;
    d = ULA.d;
    Nt = ULA.Nt;
    Nr = ULA.Nr;
    Phase_Bit = ULA.Phase_Bit;        

    % Resolution
    NQt = ULA.NQt;
    NQr = ULA.NQr;
    NQt = Nt;
    NQr = Nr;
    small_gain = 0.01;

    % Quantization of spatial domain
    partition_T = linspace(-1,1,NQt+1);
    partition_T = partition_T(1:end-1);
    partition_R = linspace(-1,1,NQr+1);
    partition_R = partition_R(1:end-1);
    AoD_Virtual_quantized_array = 2*pi*d/lambda*partition_T;
    AoA_Virtual_quantized_array = 2*pi*d/lambda*partition_R;

    % Dictionary matrix (approximating channel)
    A_Tx = zeros(Nt,NQt);
    for u = 1:1:NQt
        A_Tx_temp = 1/sqrt(Nt)*transpose(exp(-1i*AoD_Virtual_quantized_array(u).*(0:1:Nt-1)));
        A_Tx(:,u) = A_Tx_temp;
    end
    A_Rx = zeros(Nr,NQr);
    for v = 1:1:NQr
        A_Rx_temp = 1/sqrt(Nr)*transpose(exp(-1i*AoA_Virtual_quantized_array(v).*(0:1:Nr-1)));
        A_Rx(:,v) = A_Rx_temp;
    end

    % Application of position information    
    AoDA_range = [H.AoD_range H.AoA_range];  
    AoDA_range_radian = deg2rad(AoDA_range);
    AoDA_range_Virtual = 2*pi*d/lambda*sin(AoDA_range_radian);
    Quan_Pos_Err_Position_Info = zeros(2,4);
    for i=1:1:2
        [Err, Pos] = min(abs(AoD_Virtual_quantized_array - AoDA_range_Virtual(i)));
        Quan_Pos_Err_Position_Info(:,i) = [Err Pos];
    end
    for i=3:1:4
        [Err, Pos] = min(abs(AoA_Virtual_quantized_array - AoDA_range_Virtual(i)));
        Quan_Pos_Err_Position_Info(:,i) = [Err Pos];
    end 
 
    %Tx
    posTx = Quan_Pos_Err_Position_Info(2,1):1:Quan_Pos_Err_Position_Info(2,2);
    randgainTx = abs(randn(length(posTx),Mt))*20+5; 
    meanMax = mean(max(randgainTx));
    meanRestGain = (sum(sum(randgainTx))-meanMax*Mt)/Mt;
    [v, r] = max(randgainTx);
    for i=1:Mt
        randgainTx(r(i),i) = meanMax;
        temp = randgainTx([1:r(i)-1,r(i)+1:length(posTx)],i);
        temp = temp*meanRestGain/sum(temp);
        randgainTx([1:r(i)-1,r(i)+1:length(posTx)],i) = temp;
    end
    for i=1:Mt
        randgainTx(:,i) = circshift(randgainTx(:,i),(i-1)*round(length(posTx)/Mt)-r(i));
    end
    F_beam_space = ones(NQt,Mt)*small_gain; 
    F_beam_space(posTx,:) = randgainTx;
    F = pinv(A_Tx')*F_beam_space;
    F = Quantize_PS(F,Phase_Bit);
    
    %Rx
    posRx = Quan_Pos_Err_Position_Info(2,3):1:Quan_Pos_Err_Position_Info(2,4);
    randgainRx = abs(randn(length(posRx),Mr))*20+8; 
    W_beam_space = ones(NQr,Mr)*small_gain;
    W_beam_space(posRx,:) = randgainRx;
    W = pinv(A_Rx')*W_beam_space;
    W = Quantize_PS(W,Phase_Bit);
end

