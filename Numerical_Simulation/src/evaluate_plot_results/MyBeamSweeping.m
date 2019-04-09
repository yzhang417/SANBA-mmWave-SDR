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
% This function performs beam sweeping to identify the best precoder and
% combiner pair. It serves as an important benchmarking strategy for the 
% proposed non-coherent beam training algorithm.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% Mtr: a structure array that groups Mt and Mr, where Mt denotes the number
% of measurements used in the Tx side while Mr denotes the number of
% measurements used in the Rx side.
% ULA: a structure array that groups the related parameters on ULA.
% Sparse_Channel_Representation: a structure array that groups related 
% parameters on quantized sparse mmWave channel. Please refer to
% Sparse_Channel_Formulation.m for the details of its fields.
% L: Number of the dominant path.
% H: a structure array that groups related parameters on mmWave channel,
% please refer to Generate_Channel.m for details of its fields.
% SNR: SNR in dB.
% Pattern_Method_Number: it permits to choose the beam patterns that used
% in beam sweeping. When Pattern_Method_Number equals to 2, it performs
% beam sweeping with directional beam pattern that covers the angular
% domain uniformly. When Pattern_Method_Number equals to 1, it allows
% flexible configuration of different beam patterns mode by using function
% Generate_Sensing_Matrix, whose details are specified in 
% Generate_Sensing_Matrix.m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output arguments:
% f_bsweep: best precoder obtained by beam sweeping (column vector).
% w_bsweep: best combiner obtained by beam sweeping (column vector).
% Est_AoD: estimated AoD in degree.
% Est_AoA: estimated AoA in degree.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [f_bsweep, w_bsweep, Est_AoD, Est_AoA] = MyBeamSweeping(Mtr,ULA,Sparse_Channel_Representation,L,H,SNR,Pattern_Method_Number)
    %% Parameters fectching
    Mt = Mtr.Mt;
    Mr = Mtr.Mr;
    Nt = ULA.Nt;
    Nr = ULA.Nr;
    d = ULA.d;
    lambda = ULA.lambda;
    switch Pattern_Method_Number
        case 1
            Sensing_Matrix = Generate_Sensing_Matrix("Directional_Beam_Angular",Mt,Mr,ULA,Sparse_Channel_Representation,0,L,H);
        case 2
            AoD_range = H.AoD_range;
            AoA_range = H.AoA_range;
            Partition_AoD_range = linspace(AoD_range(1),AoD_range(2),Mt+1);
            Partition_AoA_range = linspace(AoA_range(1),AoA_range(2),Mr+1);
            Testing_AoD = (Partition_AoD_range(1:end-1)+Partition_AoD_range(2:end))/2;
            Testing_AoA = (Partition_AoA_range(1:end-1)+Partition_AoA_range(2:end))/2;
            F = zeros(Nt,Mt);
            W = zeros(Nr,Mr);
            for i=1:1:Mt
                F(:,i) = transpose(exp(-1i*2*pi*d/lambda*sind(Testing_AoD(i)).*(0:1:Nt-1)));
            end
            for i=1:1:Mr
                W(:,i) = transpose(exp(-1i*2*pi*d/lambda*sind(Testing_AoA(i)).*(0:1:Nr-1)));
            end
            F = Quantize_PS(F,ULA.Phase_Bit);
            W = Quantize_PS(W,ULA.Phase_Bit);
            FW = kron(transpose(F),W');
            AD = Sparse_Channel_Representation.AD;
            A = FW*AD;
            Sensing_Matrix.F = F;
            Sensing_Matrix.W = W;  
            Sensing_Matrix.FW = FW;
            Sensing_Matrix.AD = AD;        
            Sensing_Matrix.measurementMat = A; 
    end
    [measurements, ~] = Generate_Measurement(Sensing_Matrix, SNR, H, 1);
    [v,p] = max(measurements.measurements_norm_square);
    
    [vSorted, ind] = sort(measurements.measurements_norm_square,'descend');
    p = ind(1:L);
    p = p(1);
    ind_F = ceil(p/Mr);
    ind_W = p-(ind_F-1)*Mr;
    F_Set = Sensing_Matrix.F;
    W_Set = Sensing_Matrix.W;
    f_bsweep = F_Set(:,ind_F);
    w_bsweep = W_Set(:,ind_W);
    
    switch Pattern_Method_Number
        case 1
            % Estimated AoD
            angle_range = -90:0.005:90;
            gain = zeros(length(angle_range),3);
            ind = 1;
            for angle = angle_range
                Channel = transpose(exp(-1i*2*pi*d/lambda*sind(angle).*(0:1:Nt-1)));
                gain(ind,:) = abs(f_bsweep'*Channel);
                ind = ind + 1;
            end
            [maxgain, p] = max(gain);
            Est_AoD = angle_range(p);
            Est_AoD = sort(Est_AoD,'descend');
            % Estimated AoA
            ind = 1;
            for angle = angle_range
                Channel = transpose(exp(-1i*2*pi*d/lambda*sind(angle).*(0:1:Nr-1)));
                gain(ind,:) = abs(w_bsweep'*Channel);
                ind = ind + 1;
            end
            [maxgain, p] = max(gain);
            Est_AoA = angle_range(p);  
            Est_AoA = sort(Est_AoA,'descend');
        case 2
            Est_AoD = Testing_AoD(ind_F);
            Est_AoA = Testing_AoA(ind_W);
    end
end

