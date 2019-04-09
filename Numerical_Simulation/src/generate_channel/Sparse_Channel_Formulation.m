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
% This function formulates the sparse representation of the channel.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% ULA: a structure array that groups related parameters on ULA.
% L: number of the dominant path.
% H: a structure array that groups related parameters on mmWave channel,
% please refer to Generate_Channel.m for details of its fields.
% Show_leakeage_Flag: whether to show the quantization of the generated
% channel.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output arguments:
% Sparse_Channel_Representation: a structure array that groups related 
% parameters on quantized sparse mmWave channel.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Sparse_Channel_Representation = Sparse_Channel_Formulation(ULA,L,H,Show_leakeage_Flag)
    %% Parameter fetching
    lambda = ULA.lambda;
    d = ULA.d;
    Nt = ULA.Nt;
    Nr = ULA.Nr;
    NQt = ULA.NQt;
    NQr = ULA.NQr;
    
    %% Sparse formulation of mmWave channel
    % Quantization of Spatial Domain
    partition_T = linspace(-1,1,NQt+1);
    partition_T = partition_T(1:end-1);
    partition_R = linspace(-1,1,NQr+1);
    partition_R = partition_R(1:end-1);
    AoD_Virtual_quantized_array = 2*pi*d/lambda*partition_T;
    AoA_Virtual_quantized_array = 2*pi*d/lambda*partition_R;
    
    % Dictionary matrix (approximating Channel)
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
    
    %% Determine non-zero position of sparse vector z
    H_AoD_Virtual = 2*pi*d/lambda*sin(H.AoD_array_radian);
    H_AoA_Virtual = 2*pi*d/lambda*sin(H.AoA_array_radian);
    Quan_Pos_Err_h_array = zeros(2,L,2);
    for l=1:1:L
        [Err, Pos] = min(abs(AoD_Virtual_quantized_array - H_AoD_Virtual(l)));
        Quan_Pos_Err_h_array(:,l,1) = [Err, Pos];
        [Err, Pos] = min(abs(AoA_Virtual_quantized_array - H_AoA_Virtual(l)));
        Quan_Pos_Err_h_array(:,l,2) = [Err, Pos];
    end
    
    % Generate sparse vector z
    PosZ = zeros(L,1);
    for l=1:1:L
        PosZ(l) = (Quan_Pos_Err_h_array(2,l,1)-1)*NQr + Quan_Pos_Err_h_array(2,l,2);
    end
    z = zeros(NQr*NQt,1);
    z(PosZ) = H.h_array;
      
    %% Application of position information    
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

    % Reduce the dimension of sparse vector Z and Dictionary matrix AD by
    % applying position information
    pos_index = [];
    for i=Quan_Pos_Err_Position_Info(2,1):1:Quan_Pos_Err_Position_Info(2,2)
        pos_index_temp = (i-1)*NQr+Quan_Pos_Err_Position_Info(2,3):(i-1)*NQr+Quan_Pos_Err_Position_Info(2,4);
        pos_index = [pos_index pos_index_temp];
    end
    z_reduced = z(pos_index);
    AD_reduced = zeros(Nr*Nt,length(z_reduced));
    i = 1;
    for u = Quan_Pos_Err_Position_Info(2,1):Quan_Pos_Err_Position_Info(2,2)
        for v = Quan_Pos_Err_Position_Info(2,3):Quan_Pos_Err_Position_Info(2,4)
        	AD_reduced(:,i) = kron(conj(A_Tx(:,u)),A_Rx(:,v));
            i = i + 1;
        end
    end   
    AD = AD_reduced;
    z_full = z;
    z = z_reduced;
            
    %% Plot leakage
    Leakage = A_Rx'*H.H_Matrix*A_Tx;
    z_leakage = vec(Leakage);
    z_leakage_reduced = z_leakage(pos_index);
    With_Leakage_Case = abs(Leakage);
    Non_Leakage_Case = abs(reshape(z_full,ULA.NQr,ULA.NQt));
    if Show_leakeage_Flag
        figure;
        bar3(With_Leakage_Case);hold on;
        bar3(Non_Leakage_Case,'r');
        grid on
        title('Leakage');
    end
    
    %% Output result
    Sparse_Channel_Representation.AoDA_range = AoDA_range;
    Sparse_Channel_Representation.AD = AD;
    Sparse_Channel_Representation.z = z;
    Sparse_Channel_Representation.z_leakage = z_leakage;
    Sparse_Channel_Representation.z_leakage_reduced = z_leakage_reduced;
    Sparse_Channel_Representation.A_Rx = A_Rx;
    Sparse_Channel_Representation.AoD_Virtual_quantized_array = AoD_Virtual_quantized_array;
    Sparse_Channel_Representation.AoA_Virtual_quantized_array = AoA_Virtual_quantized_array;
    Sparse_Channel_Representation.Quan_Pos_Err_h_array = Quan_Pos_Err_h_array;
    Sparse_Channel_Representation.Quan_Pos_Err_Position_Info = Quan_Pos_Err_Position_Info;
end

