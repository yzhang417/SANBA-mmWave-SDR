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
% This function generates the sensing matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% Method: a structure array that groups related parameters on the
% activation setting of the different methods (algorithms) to test.
% Mt: number of measurements at Tx side.
% Mr: number of measurements at Rx side.
% ULA: a structure array that groups related parameters on ULA.
% Sparse_Channel_Representation: a structure array that groups related 
% parameters on quantized sparse mmWave channel. Please refer to 
% Sparse_Channel_Formulation.m for details of its fields.
% Show_Beam_Pattern_Flag: whether to plot the beam pattern associated with
% the generated sensing matrix.
% L: Number of the dominant path.
% H: a structure array that groups related parameters on mmWave channel,
% please refer to Generate_Channel.m for details of its fields.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output arguments:
% Sensing_Matrix: a structure array that groups related parameters sensing
% matrix, please refer to Generate_Sensing_Matrix.m for details of its
% fields.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Sensing_Matrix = Generate_Sensing_Matrix(Method,Mt,Mr,ULA,Sparse_Channel_Representation,Show_Beam_Pattern_Flag,L,H)
    %% Parameter fetching
    lambda = ULA.lambda;
    d = ULA.d;
    Nt = ULA.Nt;
    Nr = ULA.Nr;
    Phase_Bit = ULA.Phase_Bit;        
    AD = Sparse_Channel_Representation.AD;
    
    %% Generate beam patterns 
    switch Method
        case 'Random_Phase_State'                   % random gaussian sensing matrix without constraints on phase array  
            Np = Phase_Bit^2;
            randt = randi([0 Np-1],Nt,Mt);
            F = exp(1j*randt*pi/Np)/sqrt(Nt);
            randr = randi([0 Np-1],Nr,Mr);
            W = exp(1j*randr*pi/Np)/sqrt(Nr);

        case 'Gaussian_Infinite_Bits_Phase'         % random gaussian sensing matrix with infinite-bit phase shifters   
            F = (randn(Nt,Mt) + 1i*randn(Nt,Mt));
            F = exp(1j*angle(F))/sqrt(Nt);
            W = (randn(Nr,Mr) + 1i*randn(Nr,Mr));
            W = exp(1j*angle(W))/sqrt(Nr);

        case 'Gaussian_Finite_Bits_Phase'           % random gaussian sensing matrix with finite-bit phase shifters
            F = (randn(Nt,Mt) + 1i*randn(Nt,Mt));
            F = Quantize_PS(F,Phase_Bit);
            W = (randn(Nr,Mr) + 1i*randn(Nr,Mr));
            W = Quantize_PS(W,Phase_Bit);

        case 'Gaussian_Finite_Bits_Phase_Low_Rank'  % random gaussian sensing matrix with finite-bit phase shifters and with rank reduced operation
            F = (randn(Nt,Mt) + 1i*randn(Nt,Mt));
            F = Quantize_PS(F,Phase_Bit);
            rk = round(Mt*0.4);
            for i=1:1:rk
                coef = randn(Mt-rk,1);
                coef = coef/norm(coef);
                F(:,i) = F(:,rk+1:Mt)*coef;
            end
            %F = Quantize_PS(F,Phase_Bit);
            W = (randn(Nr,Mr) + 1i*randn(Nr,Mr));
            W = Quantize_PS(W,Phase_Bit);    
            rk = round(Mr*0.4);
            for i=1:1:rk
                coef = randn(Mt-rk,1);
                coef = coef/norm(coef);
                W(:,i) = W(:,rk+1:Mr)*coef;
            end            
            %W = Quantize_PS(W,Phase_Bit);
            FW = kron(transpose(F),W');
            [ui,~]=size(unique(FW,'row'));
            fprintf('Rank of FW %d \n\n', rank(FW));
            fprintf('Unique row of FW %d \n\n', ui);
       
        case 'Region_Random_Beam'                   % random beam pattern within a certain region
            [F, W] = Region_Random_Beam(Mt,Mr,ULA,H);
            
        case 'Directional_Random_Beam'              % directional beam pattern with random gain in spatial domain
            [F, W] = Directional_Random_Beam(Mt,Mr,ULA,H);

        case 'Directional_Beam'                     % directional beam pattern with uniform gain in spatial domain
            Rank_Eliminated = 0;
            [F, W] = Directional_Beam(Mt,Mr,ULA,H,Rank_Eliminated);
            
        case 'Directional_Beam_Angular'             % directional beam pattern with uniform gain in angular domain
            [F, W] = Directional_Beam_Angular(Mt,Mr,ULA,H);            
            
        case 'Directional_Beam_Low_Rank'            % directional beam pattern with uniform gain in spatial domain and with reduced rank operation
            Rank_Eliminated = 2;
            [F, W] = Directional_Beam(Mt,Mr,ULA,H,Rank_Eliminated);
               
        case 'ProtoType_Pattern'                    % using real beam pattern that is used by the phased array hardware
            load('prmCode.mat');
            load('BeamformerQuant_F.mat');          % This mat file stores the beamformer
            F = BeamformerQuant_F(:,:,prmCode.Ind);
            load('BeamformerQuant_W.mat');          % This mat file stores the combiner
            W = BeamformerQuant_W(:,:,prmCode.Ind);
        otherwise
            ;
    end
    if Show_Beam_Pattern_Flag
        show_beam_pattern(lambda, d, F);
        show_beam_pattern(lambda, d, W);
    end
    
    %% Output result
    if rank(F) < min(Mt,Nt)
        fprintf('Repeat beam pattern with rank of F is %d \n', rank(F));
    end
    if rank(W) < min(Mr,Nr)
        fprintf('Repeat beam pattern with rank of W is %d \n', rank(W));
    end        
    FW = kron(transpose(F),W');
    A = FW*AD;
    Sensing_Matrix.F = F;
    Sensing_Matrix.W = W;  
    Sensing_Matrix.FW = FW;
    Sensing_Matrix.AD = AD;        
    Sensing_Matrix.measurementMat = A;    
end


