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
% This function generates directional precoder and combiner matrix. The
% generated beam patterns uniformly cover the spatial domain rather than
% angular domain.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% Mt: number of measurements at Tx side.
% Mr: number of measurements at Rx side.
% ULA: a structure array that groups related parameters on ULA.
% H: a structure array that groups related parameters on mmWave channel,
% please refer to Generate_Channel.m for details of its fields.
% Rank_Eliminated: number of rank to be eliminated for both F and W. This
% is motivated by that low rank is preferred by our proposed two-stage
% algorithm.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output arguments:
% F: precoder matrix (each column represents one beam pattern codebook).
% W: combiner matrix (each column represents one beam pattern codebook).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [F, W] = Directional_Beam(Mt,Mr,ULA,H,Rank_Eliminated)
    % Parameter fetching
    lambda = ULA.lambda;
    d = ULA.d;
    Nt = ULA.Nt;
    Nr = ULA.Nr;
    Phase_Bit = ULA.Phase_Bit;        

    % Resolution
    NQt = 20*ULA.NQt;
    NQr = 20*ULA.NQr;
    small_gain = 0.05;

    % Number of indenpendent pattern
    Rank_Eliminated = min(Rank_Eliminated,Mt-3);
    Rank_Eliminated = max(Rank_Eliminated,0);
    Mt = Mt - Rank_Eliminated;
    Mr = Mr - Rank_Eliminated;

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

    % Tx   
    posTx = Quan_Pos_Err_Position_Info(2,1):1:Quan_Pos_Err_Position_Info(2,2);
    n_sub_grid_t = ceil(length(posTx)/Mt);
    n_overlap = n_sub_grid_t*Mt - length(posTx);
    n_overlap_left = ceil(n_overlap/2);
    n_overlap_right = floor(n_overlap/2);
    overlap_index = [1:n_overlap_left, Mt-n_overlap_right:Mt-1];
    gainTx = ones(length(posTx),Mt)*small_gain;
    start_point = 1;
    for i=1:1:Mt
        gainTx(start_point:(start_point+n_sub_grid_t-1),i) = 1;
        if ismember(i,overlap_index)
            start_point = start_point+n_sub_grid_t - 1;
        else
            start_point = start_point+n_sub_grid_t;
        end
    end   
    F_beam_space = zeros(NQt,Mt);
    F_beam_space(posTx,:) = gainTx;
    F = pinv(A_Tx')*F_beam_space;
    F = F/norm(F);
    F = Quantize_PS(F,Phase_Bit);
    
    % Rx
    posRx = Quan_Pos_Err_Position_Info(2,3):1:Quan_Pos_Err_Position_Info(2,4);
    n_sub_grid_r = ceil(length(posRx)/Mr);
    n_overlap = n_sub_grid_r*Mr - length(posRx);
    n_overlap_left = ceil(n_overlap/2);
    n_overlap_right = floor(n_overlap/2);
    overlap_index = [1:n_overlap_left, Mr-n_overlap_right:Mr-1];
    gainRx = ones(length(posRx),Mr)*small_gain;
    start_point = 1;
    for i=1:1:Mr
        gainRx(start_point:(start_point+n_sub_grid_r-1),i) = 1;
        if ismember(i,overlap_index)
            start_point = start_point+n_sub_grid_r - 1;
        else
            start_point = start_point+n_sub_grid_r;
        end
    end    
    W_beam_space = zeros(NQr,Mr);
    W_beam_space(posRx,:) = gainRx;
    W = pinv(A_Rx')*W_beam_space;
    W = W/norm(W);
    W = Quantize_PS(W,Phase_Bit);  

    % Randomly add correlated beam pattern
    if Rank_Eliminated > 0
        Fc_ind = datasample(1:Mt,min(Rank_Eliminated*2,Mt),'Replace',false);
        Wc_ind = datasample(1:Mr,min(Rank_Eliminated*2,Mr),'Replace',false);
        for i=1:Rank_Eliminated
            F_Correlated = (F(:,Fc_ind(i))+F(:,Fc_ind(i+1)));
            F = [F F_Correlated];
            W_Correlated = (W(:,Wc_ind(i))+W(:,Wc_ind(i+1)));
            W = [W W_Correlated];
        end  
    end
end