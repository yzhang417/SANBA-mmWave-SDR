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
% This function generates the measurements according to the generated
% channel and sensing matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% Sensing_Matrix: a structure array that groups related parameters sensing
% matrix, please refer to Generate_Sensing_Matrix.m for details of its
% fields.
% SNR: the SNR before beam training in dB.
% H: a structure array that groups related parameters on mmWave channel.
% Add_Noise_Flag: whether the Gaussian noise is introduced in the Rx side.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output arguments:
% measurements: a structure array that groups related parameters on the 
% received measurements with different transformations. 
% noise_power: the power of noise calculated according to SNR.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [measurements, noise_power] = Generate_Measurement(Sensing_Matrix, SNR, H, Add_Noise_Flag)
    Calibration_Error = 0; % Whether introduce calibration error
    if Calibration_Error
        Calibration_Noise = exp(-1i.*unifrnd(0,90,length(Sensing_Matrix.F(:,1)),length(Sensing_Matrix.F(1,:))).*pi/180);
        F = Sensing_Matrix.F.*Calibration_Noise;
        W = Sensing_Matrix.W.*Calibration_Noise;
        FW = kron(transpose(F),W');
    else
        FW = Sensing_Matrix.FW;
    end
    measurementMat = Sensing_Matrix.measurementMat;
    vecH = H.vecH;
    if Add_Noise_Flag
        % Number of measurements
        m = length(measurementMat(:,1));

        % Signal power is normalized to one already
        signal_power = 1;
        
        % Noise power
        noise_power = signal_power/(10^(SNR/10));

        % IID Additive Gaussian Noise
        [Nrx, Mrx] = size(Sensing_Matrix.W);
        [Ntx, Mtx] = size(Sensing_Matrix.F);
        noiseMatrix = sqrt(noise_power) * (randn(Nrx,Mrx)+1j*randn(Nrx,Mrx)) * (1/sqrt(2));
        noiseVec = [];
        for ntx=1:Mtx
            noiseSubVec = diag(Sensing_Matrix.W' * noiseMatrix);
            noiseVec = [noiseVec;noiseSubVec];
        end       
        
        % Output noisy measurement
        measurements_with_perfect_phase = FW*vecH + noiseVec;
        measurements_norm_square = abs(measurements_with_perfect_phase).^2;
        measurements_with_noisy_phase = measurements_with_perfect_phase .* ... 
            (randn(length(measurements_with_perfect_phase),1)+1j*randn(length(measurements_with_perfect_phase),1))*1/sqrt(2);
        ISNR = (abs(FW*vecH).^2./abs(noiseVec).^2);
    else
        % Output noiseless measurement
        noise_power = 10^(-10);
        ISNR = inf;
        measurements_norm_square = abs(FW*vecH).^2;   
        measurements_with_perfect_phase = FW*vecH;
        measurements_with_noisy_phase = measurements_with_perfect_phase .* ... 
            (randn(length(measurements_with_perfect_phase),1)+1j*randn(length(measurements_with_perfect_phase),1))*1/sqrt(2);    
    end
    measurements.measurements_norm_square = measurements_norm_square;
    measurements.measurements_with_perfect_phase = measurements_with_perfect_phase;
    measurements.measurements_with_noisy_phase = measurements_with_noisy_phase;
    measurements.ISNR = ISNR; %Instanteous SNR
end

