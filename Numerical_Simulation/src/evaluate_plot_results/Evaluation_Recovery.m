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
% This function helps to plot the simulation results of the proposed
% algorithm and the related benchmarking algorithms.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% recoveredSig: X-axis label (character array).
% H: a structure array that groups related parameters on mmWave channel,
% please refer to Generate_Channel.m for details of its fields.
% Sparse_Channel_Representation: a structure array that groups related 
% parameters on quantized sparse mmWave channel. Please refer to
% Sparse_Channel_Formulation.m for the details of its fields.
% ULA: a structure array that groups the related parameters on ULA.
% SNR: SNR in dB.
% Mtr: a structure array that groups Mt and Mr, where Mt denotes the
% number.
% of measurements used in the Tx side while Mr denotes the number of
% measurements used in the Rx side.
% noise_power: power of noise.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output arguments:
% Evaluation_result: a structure array that groups the results with 
% different evaluation metrics.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Evaluation_result = Evaluation_Recovery(recoveredSig,L,H,Sparse_Channel_Representation,ULA,SNR,Mtr,noise_power)
    %% Parameter fetch
    Quan_Pos_Err_Position_Info = Sparse_Channel_Representation.Quan_Pos_Err_Position_Info;
    AoD_Virtual_quantized_array = Sparse_Channel_Representation.AoD_Virtual_quantized_array;
    AoA_Virtual_quantized_array = Sparse_Channel_Representation.AoA_Virtual_quantized_array;
    Quan_Pos_Err_h_array = Sparse_Channel_Representation.Quan_Pos_Err_h_array;
    lambda = ULA.lambda;
    d = ULA.d;
    Nt = ULA.Nt;
    Nr = ULA.Nr;
    Phase_Bit = ULA.Phase_Bit;
    
    %% Extract the L maximum values of z
    [recoveredSig_sorted, ind] = sort(recoveredSig,'descend');
    AoDA_recovered = ind(1:L);
    AoDA_recovered = sort(AoDA_recovered,'ascend');
    
    %% Estimated AoD and AoA
    AoD_Estimated = zeros(1,L);
    AoA_Estimated = zeros(1,L);
    AoD_True_Quantized = zeros(1,L);
    AoA_True_Quantized = zeros(1,L);
    AoD_True = H.AoD_array;
    AoA_True = H.AoA_array;
    Avai_AoD_Index = Quan_Pos_Err_Position_Info(2,1):Quan_Pos_Err_Position_Info(2,2);
    Avai_AoA_Index = Quan_Pos_Err_Position_Info(2,3):Quan_Pos_Err_Position_Info(2,4);
    Num_Avai_AoD = length(Avai_AoD_Index);
    Num_Avai_AoA = length(Avai_AoA_Index);  
    for l=1:L
        ind_AoD = ceil(AoDA_recovered(l)/Num_Avai_AoA);
        ind_AoA = AoDA_recovered(l) - (ind_AoD-1)*Num_Avai_AoA;     
        
        AoD_Estimated_Virtual = AoD_Virtual_quantized_array(Avai_AoD_Index(ind_AoD));
        AoD_Estimated_Virtual = AoD_Estimated_Virtual/(2*pi*d/lambda);
        AoD_Estimated_In_Deg = rad2deg(asin(AoD_Estimated_Virtual));
        AoD_Estimated(l) = AoD_Estimated_In_Deg;
        
        AoA_Estimated_Virtual = AoA_Virtual_quantized_array(Avai_AoA_Index(ind_AoA));
        AoA_Estimated_Virtual = AoA_Estimated_Virtual/(2*pi*d/lambda);
        AoA_Estimated_In_Deg = rad2deg(asin(AoA_Estimated_Virtual));
        AoA_Estimated(l) = AoA_Estimated_In_Deg;
                     
        AoD_True_Quantized_Virtual = AoD_Virtual_quantized_array(Quan_Pos_Err_h_array(2,l,1));
        AoD_True_Quantized_Virtual = AoD_True_Quantized_Virtual/(2*pi*d/lambda);
        AoD_True_Quantized_In_Deg = rad2deg(asin(AoD_True_Quantized_Virtual));
        AoD_True_Quantized(l) = AoD_True_Quantized_In_Deg;
 
        AoA_True_Quantized_Virtual = AoA_Virtual_quantized_array(Quan_Pos_Err_h_array(2,l,2));
        AoA_True_Quantized_Virtual = AoA_True_Quantized_Virtual/(2*pi*d/lambda);
        AoA_True_Quantized_In_Deg = rad2deg(asin(AoA_True_Quantized_Virtual));
        AoA_True_Quantized(l) = AoA_True_Quantized_In_Deg;
    end    
        
    % Re-order the AoD and AoA for correctly computing the error of AoD and AOA
    [AoD_True, order_True] = sort(AoD_True,'descend');
    AoA_True(order_True) = AoA_True(order_True);
    AoD_True_Quantized = AoD_True_Quantized(order_True);
    AoA_True_Quantized = AoA_True_Quantized(order_True);
    [AoD_Estimated, order_Estimated] = sort(AoD_Estimated,'descend');
    AoA_Estimated = AoA_Estimated(order_Estimated);
       
    % Compute the AoD and AoA estimation error
    AoD_Err_to_True_Quantized = mean(abs(AoD_Estimated-AoD_True_Quantized));
    AoA_Err_to_True_Quantized = mean(abs(AoA_Estimated-AoA_True_Quantized));
    AoD_Err_to_True = mean(abs(AoD_Estimated-AoD_True)); 
    AoA_Err_to_True = mean(abs(AoA_Estimated-AoA_True));
    AoDA_Err_Quantized = mean([AoD_Err_to_True_Quantized AoA_Err_to_True_Quantized]);
    if Nt ==1 
        AoDA_Err = AoA_Err_to_True;
    elseif Nr ==1 
        AoDA_Err = AoD_Err_to_True;
    else
        AoDA_Err = mean([AoD_Err_to_True AoA_Err_to_True]);
    end
    
    %% Calculate the AoD and AoA Error In Quantization Level (Success Rate of Recovery)
    if L == 1
        Num_Quantization_Error_AoD = abs(Quan_Pos_Err_h_array(2,l,1) - Avai_AoD_Index(ind_AoD));
        Num_Quantization_Error_AoA = abs(Quan_Pos_Err_h_array(2,l,2) - Avai_AoA_Index(ind_AoA));
        if Nt ==1 
            Num_Quantization_Error = Num_Quantization_Error_AoA;
        elseif Nr ==1 
            Num_Quantization_Error = Num_Quantization_Error_AoD;
        else
            Num_Quantization_Error = ( Num_Quantization_Error_AoD + Num_Quantization_Error_AoA )/2;
        end
    else
        Num_Quantization_Error = 0;
    end
    
    %% Compute the normalized mean squred error of channel matrix H and array response vector
    % MSE of Array response vector (ARV)
    ATx_Perfect = zeros(Nt,L);
    ARx_Perfect = zeros(Nr,L);
    ATx_Estimated = zeros(Nt,L);
    ARx_Estimated = zeros(Nr,L);
    for l=1:1:L
        ATx_Perfect(:,l) = transpose(exp(-1i*2*pi/lambda*d*sind(AoD_True(l)).*(0:1:Nt-1)))*1/sqrt(Nt);
        ARx_Perfect(:,l) = transpose(exp(-1i*2*pi/lambda*d*sind(AoA_True(l)).*(0:1:Nr-1)))*1/sqrt(Nr);
        ATx_Estimated(:,l) = transpose(exp(-1i*2*pi/lambda*d*sind(AoD_Estimated(l)).*(0:1:Nt-1)))*1/sqrt(Nt);
        ARx_Estimated(:,l) = transpose(exp(-1i*2*pi/lambda*d*sind(AoA_Estimated(l)).*(0:1:Nr-1)))*1/sqrt(Nr);
    end
    
    if Nt ==1 
        MSE_ARM = norm(ARx_Perfect-ARx_Estimated,'fro')^2/norm(ARx_Perfect,'fro')^2;
    elseif Nr ==1 
        MSE_ARM = norm(ATx_Perfect-ATx_Estimated,'fro')^2/norm(ATx_Perfect,'fro')^2;
    else
        MSE_ARM = (norm(ATx_Perfect-ATx_Estimated,'fro')^2/norm(ATx_Perfect,'fro')^2 + ...
           norm(ARx_Perfect-ARx_Estimated,'fro')^2/norm(ARx_Perfect,'fro')^2)/2;
    end
        
    % MSE of Channel Matrix
    % Correct the global phase factor
    vecH_Estimated = Sparse_Channel_Representation.AD*recoveredSig;
    phaseFac = exp( 1i* angle( (vecH_Estimated'*H.vecH)/(H.vecH'*H.vecH) ) );
    vecH_Estimated = vecH_Estimated*phaseFac;
    H_estimated = reshape(vecH_Estimated,Nr,Nt);
    MSE_H = norm(H_estimated-H.H_Matrix,'fro')^2/norm(H.H_Matrix,'fro')^2;

    %% Spectrum Efficiency
    Num_RFchains = 1;
    % Unconstrained SVD with Perfect CSI
    [U,S,V] = svd(H.H_Matrix);
    
    % With Perfect CSI
    if Num_RFchains == 1 || min(ULA.Nr,ULA.Nt) == 1 || ULA.Nr == ULA.Nt
        w_Unconstrained = U(:,1);
        f_Unconstrained = V(:,1);
        f_PerCSI = Quantize_PS(f_Unconstrained,Phase_Bit);
        w_PerCSI = Quantize_PS(w_Unconstrained,Phase_Bit);
    else
        w_Unconstrained = U(:,1:L);
        f_Unconstrained = repmat(V(:,1),1,L);
        f_PerCSI = HybridPrecoding(f_Unconstrained,ULA.Nt,Num_RFchains,Phase_Bit);
        w_PerCSI = HybridPrecoding(w_Unconstrained,ULA.Nr,Num_RFchains,Phase_Bit);
    end
    
    % Estimated Channel
    [U,S,V] = svd(H_estimated);
    if Num_RFchains == 1 || min(ULA.Nr,ULA.Nt) == 1 || ULA.Nr == ULA.Nt
        f_Est_Infinite_Bit = V(:,1);
        w_Est_Infinite_Bit = U(:,1);
        f_Est = Quantize_PS(f_Est_Infinite_Bit,Phase_Bit);
        w_Est = Quantize_PS(w_Est_Infinite_Bit,Phase_Bit);
    else
        f_Est_Infinite_Bit = repmat(V(:,1),1,L);
        w_Est_Infinite_Bit = U(:,1:L); 
        f_Est = HybridPrecoding(f_Est_Infinite_Bit,ULA.Nt,Num_RFchains,Phase_Bit);
        w_Est = HybridPrecoding(w_Est_Infinite_Bit,ULA.Nr,Num_RFchains,Phase_Bit);
    end

    % Picking up the single dominant path with the estimated Channel
    f_Infinite_PS_Bits = transpose(exp(-1i*2*pi/lambda*d*sind(AoD_Estimated(l)).*(0:1:Nt-1)))*1/sqrt(Nt);
    w_Infinite_PS_Bits = transpose(exp(-1i*2*pi/lambda*d*sind(AoA_Estimated(l)).*(0:1:Nr-1)))*1/sqrt(Nr);
    f_Est_AoD = Quantize_PS(f_Infinite_PS_Bits,Phase_Bit); 
    w_Est_AoA = Quantize_PS(w_Infinite_PS_Bits,Phase_Bit);
    
    % Beam sweeping with directional beam pattern
    [f_bsweep1, w_bsweep1, Est_AoD1, Est_AoA1] = MyBeamSweeping(Mtr,ULA,Sparse_Channel_Representation,L,H,SNR,1);
    [f_bsweep2, w_bsweep2, Est_AoD2, Est_AoA2] = MyBeamSweeping(Mtr,ULA,Sparse_Channel_Representation,L,H,SNR,2); 
    if Num_RFchains == 1 || min(ULA.Nr,ULA.Nt) == 1 || ULA.Nr == ULA.Nt
        f_bsweep1 = Quantize_PS(f_bsweep1,Phase_Bit);
        w_bsweep1 = Quantize_PS(w_bsweep1,Phase_Bit);
        f_bsweep2 = Quantize_PS(f_bsweep2,Phase_Bit);
        w_bsweep2 = Quantize_PS(w_bsweep2,Phase_Bit);
    else
        f_bsweep1 = HybridPrecoding(f_bsweep1,ULA.Nt,Num_RFchains,Phase_Bit);
        w_bsweep1 = HybridPrecoding(w_bsweep1,ULA.Nr,Num_RFchains,Phase_Bit);
        f_bsweep2 = HybridPrecoding(f_bsweep2,ULA.Nt,Num_RFchains,Phase_Bit);
        w_bsweep2 = HybridPrecoding(w_bsweep2,ULA.Nr,Num_RFchains,Phase_Bit);
    end
    
    AoD_Err_Sweep1 = mean(abs(Est_AoD1-AoD_True));
    AoA_Err_Sweep1 = mean(abs(Est_AoA1-AoA_True));
    AoD_Err_Sweep2 = mean(abs(Est_AoD2-AoD_True));
    AoA_Err_Sweep2 = mean(abs(Est_AoA2-AoA_True));
    if Nt ==1 
        AoDA_Err_Sweep1 = AoA_Err_Sweep1;
        AoDA_Err_Sweep2 = AoA_Err_Sweep2;
    elseif Nr ==1 
        AoDA_Err_Sweep1 = AoD_Err_Sweep1;
        AoDA_Err_Sweep2 = AoD_Err_Sweep2;
    else
        AoDA_Err_Sweep1 = (AoD_Err_Sweep1+AoA_Err_Sweep1)/2;
        AoDA_Err_Sweep2 = (AoD_Err_Sweep2+AoA_Err_Sweep2)/2;
    end

    %% Data saved for future use (prototyping)
    save w_Est_Infinite_Bit w_Est_Infinite_Bit
    save f_Est_Infinite_Bit f_Est_Infinite_Bit
    save AoD_Estimated AoD_Estimated
    save AoA_Estimated AoA_Estimated
    
    %% Gain and nosie power
    Unconstrained_gain = w_Unconstrained'*H.H_Matrix*f_Unconstrained;   
    PerCSI_gain = w_PerCSI'*H.H_Matrix*f_PerCSI;
    Est_gain = w_Est'*H.H_Matrix*f_Est;
    Est_gain_AoDA = w_Est_AoA'*H.H_Matrix*f_Est_AoD;
    BSweep_gain1 = w_bsweep1'*H.H_Matrix*f_bsweep1;
    BSweep_gain2 = w_bsweep2'*H.H_Matrix*f_bsweep2;
    PowerNoise = noise_power;
    
    % Compute Rate
    Rate_Unconstrained = abs(log2(det(eye(L)+Unconstrained_gain*Unconstrained_gain'/PowerNoise)));
    Rate_PerCSI = abs(log2(det(eye(L)+PerCSI_gain*PerCSI_gain'/PowerNoise)));
    Rate_Est = abs(log2(det(eye(L)+Est_gain*Est_gain'/PowerNoise)));
    Rate_Est_AoDA = abs(log2(det(eye(L)+Est_gain_AoDA*Est_gain_AoDA'/PowerNoise)));
    Rate_BSweep1 = abs(log2(det(eye(L)+BSweep_gain1*BSweep_gain1'/PowerNoise)));
    Rate_BSweep2 = abs(log2(det(eye(L)+BSweep_gain2*BSweep_gain2'/PowerNoise)));
    
    if max(Rate_BSweep2,Rate_BSweep1) > Rate_Unconstrained
        stop = 1;
    end
    
    % Output Result
    AoD_Estimated = AoD_Estimated;
    AoA_Estimated = AoA_Estimated;
    Evaluation_result=[Num_Quantization_Error;...
                       AoDA_Err;...
                       MSE_ARM;...
                       MSE_H;...
                       Rate_Unconstrained;...
                       Rate_PerCSI;...
                       Rate_Est;...
                       Rate_Est_AoDA;...
                       Rate_BSweep1;...
                       Rate_BSweep2;...
                       AoDA_Err_Sweep1;...
                       AoDA_Err_Sweep2];
     save AoDA_Err AoDA_Err;
end

