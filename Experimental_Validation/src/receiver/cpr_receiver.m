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
% This function performs three major tasks illustrated as below:
% Task I: running the Rx USRP to receive the signal over the air and 
% fetching the digital data stream from the Rx USRP. To this end, the
% necessary Tx/Rx phased array configurations are performed. 
% Task II: after collecting the measurements, this function runs the 
% proposed two-stage non-coherent beam alignment algorithm to estimate the 
% AoA and AoD. 
% Task III: to further justify the correct implementation and evaluate the 
% algorithm's performance. This function tests the estimated best 
% beamformer and combiner pair over the air and compare its RSSI result 
% with the beam sweeping strategy. The collected data is saved in 
% [Experimental_Validation\data\val_cpr_data].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% AoD_True: true AoD, measured manually.
% AoA_True: true AoA, measured manually.
% compileIt: whether to compile the run_usrp_rx.m again for Codegen.
% num_codebook_to_train: number of beam patterns to be probed at each side.
% G: quantization level.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output arguments:
% There is no direct output argument from the function, but the following
% variable is store for future demonstration (plotted by functions in 
% the folder [Experimental_Validation\main_programs\plot_results\].
% BeamformerQuant_W: quantized best combiner codebook
% BeamformerQuant_F: quantized best beamformer codebook
% cpr_result:
%    cpr_result.AoA_Estimated: estimated AoA by non-coherent algorithm.
%    cpr_result.AoD_Estimated: estimated AoD by non-coherent algorithm.
%    cpr_result.AoA_Estimated_Sweep: estimated AoA by beam sweeping.
%    cpr_result.AoD_Estimated_Sweep: estimated AoD by beam sweeping.
%    cpr_result.Training_measurement = RSSI of beam sweeping results.
%    cpr_result.Testing_measurement = RSSI of non-coherent result.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function cpr_receiver(AoD_True, AoA_True, compileIt, calibration_method_id, num_codebook_to_train, G)
    %% Initialization 
    address = '192.168.20.2';
    platform = 'N200/N210/USRP2';
    SerialPortTx = ['/dev/tty.usbmodem3441'];
    SerialPortRx = ['/dev/tty.usbmodem221'];
    Program_ID = 8.1;

    %% System and receiver parameters
    prmQPSKReceiver = sdruqpskreceiver_init(platform);
    prmQPSKReceiver.Platform = platform;
    prmQPSKReceiver.Address = address;  
    phase_array_id_rx = 2; 
    phase_array_id_tx = 3;

    %% Set Tx Rx training beam patterns
    M = num_codebook_to_train;
    switch M
        case 6
            Searching_Area = 45;
        case 7
            Searching_Area = 47.5;
        case 8
            Searching_Area = 50;
    end
    Angle_range = [-Searching_Area/2 Searching_Area/2];
    Partition_Angle_range = linspace(Angle_range(1),Angle_range(2),M+1);
    Test_Angles = (Partition_Angle_range(1:end-1)+Partition_Angle_range(2:end))/2;
    
    %%
    [prmPAControl, type_prmPAControl] = phase_array_control_init(prmQPSKReceiver, SerialPortTx, SerialPortRx, Program_ID, 'num_codebook_to_train', num_codebook_to_train);
    type_portion_of_boundary_frame_to_remove = coder.newtype('double',[1 1]);
    %
    Output_tx = get_beam_pattern(Test_Angles,phase_array_id_tx,false);
    BeamformerQuant_F = Output_tx.BeamformerQuant_F;
    save('BeamformerQuant_F.mat','BeamformerQuant_F');
    %
    Output_rx = get_beam_pattern(Test_Angles,phase_array_id_rx,false);
    BeamformerQuant_W = Output_rx.BeamformerQuant_W;
    save('BeamformerQuant_W.mat','BeamformerQuant_W');
    %
    switch calibration_method_id
        case 1
            employed_codebook_set_tx = Output_tx.Phvec_command_no_calibration;
            employed_codebook_set_rx = Output_rx.Phvec_command_no_calibration;
        case 2
            employed_codebook_set_tx = Output_tx.Phvec_command_coarse;
            employed_codebook_set_rx = Output_rx.Phvec_command_coarse;
        case 3
            employed_codebook_set_tx = Output_tx.Phvec_command_fine_m1;
            employed_codebook_set_rx = Output_rx.Phvec_command_fine_m1;
        case 4
            employed_codebook_set_tx = Output_tx.Phvec_command_fine_m2;
            employed_codebook_set_rx = Output_rx.Phvec_command_fine_m2;
    end
    prmPAControl.CodeBook_Tx(1:num_codebook_to_train,:) = employed_codebook_set_tx;
    prmPAControl.CodeBook_Rx(1:num_codebook_to_train,:) = employed_codebook_set_rx;

    %% Global storage for corrupted signal and calculated BERMER
    global Raw_corruptSignal
    SamplesPerFrame = prmQPSKReceiver.FrameSize * prmQPSKReceiver.Upsampling * prmQPSKReceiver.RxBufferedFrames;
    Raw_corruptSignal = complex(zeros(...
        SamplesPerFrame,...
        int32(prmPAControl.Number_of_Frames_Per_Test_Rx_Buffered),...
        prmPAControl.Maximum_Number_State_To_Test));

    %% Clearance
    clear run_usrp_rx_mex %#ok<UNRCH>
    clear Control_Phase_Array_External_mex %#ok<UNRCH>
    clear run_offline_decoder_mex %#ok<UNRCH>

    %% Compilation
    if compileIt
        codegen('run_usrp_rx', 'PATR.c', 'PATR.h', '-args', {coder.Constant(prmQPSKReceiver), type_prmPAControl},...
                                                   '-I', char(strcat(pwd,'/Experimental_Validation/src/array_control')));%#ok<UNRCH>
        codegen('Control_Phase_Array_External', 'PATR.c', 'PATR.h', '-args', {type_prmPAControl},...
                                                   '-I', char(strcat(pwd,'/Experimental_Validation/src/array_control')));%#ok<UNRCH>
        codegen('run_offline_decoder', '-args', {coder.Constant(prmQPSKReceiver), type_prmPAControl, type_portion_of_boundary_frame_to_remove});%#ok<UNRCH>       
    end

    %% Data collection
    Number_Test = 1;
    for nthTest = 1 : Number_Test
        disp('-------------------------------------------------------------')
        fprintf('-----------%d-th collection of received raw signal------------\n',nthTest);
        disp('-------------------------------------------------------------')
        for nthTxBeam = 1:1:prmPAControl.Number_CodeBook_Tx_to_train
            prmPAControl.nthTxBeam = nthTxBeam;
            %clear Control_Phase_Array_External_mex %#ok<UNRCH>
            Control_Phase_Array_External_mex(prmPAControl);
            clear run_usrp_rx_mex %#ok<UNRCH>
            run_usrp_rx_mex(prmQPSKReceiver, prmPAControl);
            fprintf('\n-----------Saving data (%d-th transmitter beaming pattern applied)-----------  \n', nthTxBeam);
            save(['Raw_corruptSignal' num2str(nthTest) '_' num2str(nthTxBeam) '.mat'],'Raw_corruptSignal');
            fprintf('\n-----------Saved data (%d-th transmitter beaming pattern applied)----------- %d \n', nthTxBeam)
        end
    end
    save('prmQPSKReceiver.mat','prmQPSKReceiver');
    save('prmPAControl.mat','prmPAControl');

    %% Decoding and save measurements from prototype
    prmPAControl.nthTxBeam = num_codebook_to_train;
    save('prmPAControl.mat','prmPAControl');
    cpr_decoder;
    BMP = permute(BER_MER_Power_Overall_Training,[3 2 1]);
    MeasurementsPrototype.Data = BMP(:,:,7);
    MeasurementsPrototype.Program_ID = prmPAControl.Program_ID;
    maximum = max(max(MeasurementsPrototype.Data));
    [x,y]=find(MeasurementsPrototype.Data==maximum);
    AoD_Estimated_Sweep = Test_Angles(x);
    AoA_Estimated_Sweep = Test_Angles(y);
    MeasurementsPrototype_Train = MeasurementsPrototype;

    %% Run CPR to estimate the AoD and AoA
    prmCode.Shift = 0;
    prmCode.Ind = calibration_method_id;
    save('prmCode.mat','prmCode');
    re_estimate = 0;
    if re_estimate
        BeamformerQuant_W = BeamformerQuant_W(:,1:num_codebook_to_train,:);
        BeamformerQuant_F = BeamformerQuant_F(:,1:num_codebook_to_train,:);
        save BeamformerQuant_W BeamformerQuant_W
        save BeamformerQuant_F BeamformerQuant_F
    end
    [AoD_Estimated, AoA_Estimated] = run_cpr(MeasurementsPrototype_Train,num_codebook_to_train,prmCode,Searching_Area,G)

    %% Generate the best tx and rx beam and write to prmPAControl
    num_codebook_to_train_new = 2;
    nthRxBeam = prmPAControl.nthRxBeam;
    nthTxBeam = prmPAControl.nthTxBeam;
    [prmPAControl, ~] = phase_array_control_init(prmQPSKReceiver, SerialPortTx, SerialPortRx, Program_ID, 'num_codebook_to_train', num_codebook_to_train_new);
    prmPAControl.nthRxBeam = nthRxBeam;
    prmPAControl.nthTxBeam = nthTxBeam;
    %
    Output_tx_new = get_beam_pattern(AoD_Estimated,phase_array_id_tx,true);
    BeamformerQuant_F = [BeamformerQuant_F Output_tx_new.BeamformerQuant_F];
    save('BeamformerQuant_F.mat','BeamformerQuant_F');
    % prmPAControl.CodeBook_Tx(1:num_codebook_to_train,:) = Output_tx.Phvec_command_fine_m1;
    % prmPAControl.CodeBook_Tx(num_codebook_to_train+1,:) = Output_tx_new.Phvec_command_fine_m1;
    % prmPAControl.CodeBook_Tx(num_codebook_to_train+2,:) = Output_tx_new.Command_F;
    Output_rx_new = get_beam_pattern(AoA_Estimated,phase_array_id_rx,true);
    BeamformerQuant_W = [BeamformerQuant_W Output_rx_new.BeamformerQuant_W];
    save('BeamformerQuant_W.mat','BeamformerQuant_W');
    % prmPAControl.CodeBook_Rx(1:num_codebook_to_train,:) = Output_rx.Phvec_command_fine_m1;
    % prmPAControl.CodeBook_Rx(num_codebook_to_train+1,:) = Output_rx_new.Phvec_command_fine_m1;
    % prmPAControl.CodeBook_Rx(num_codebook_to_train+2,:) = Output_rx_new.Command_W;
    switch calibration_method_id
        case 1
            estimated_best_angle_tx = Output_tx_new.Phvec_command_no_calibration;
            estimated_best_beam_tx = Output_tx_new.Command_F1;
            estimated_best_angle_rx = Output_rx_new.Phvec_command_no_calibration;
            estimated_best_beam_rx = Output_rx_new.Command_W1;
        case 2
            estimated_best_angle_tx = Output_tx_new.Phvec_command_coarse;
            estimated_best_beam_tx = Output_tx_new.Command_F2;
            estimated_best_angle_rx = Output_rx_new.Phvec_command_coarse;
            estimated_best_beam_rx = Output_rx_new.Command_W2;
        case 3
            estimated_best_angle_tx = Output_tx_new.Phvec_command_fine_m1;
            estimated_best_beam_tx = Output_tx_new.Command_F3;
            estimated_best_angle_rx = Output_rx_new.Phvec_command_fine_m1;
            estimated_best_beam_rx = Output_rx_new.Command_W3;
        case 4
            estimated_best_angle_tx = Output_tx_new.Phvec_command_fine_m2;
            estimated_best_beam_tx = Output_tx_new.Command_F4;
            estimated_best_angle_rx = Output_rx_new.Phvec_command_fine_m2;
            estimated_best_beam_rx = Output_rx_new.Command_W4;
    end
    prmPAControl.CodeBook_Tx(1,:) = estimated_best_angle_tx;
    prmPAControl.CodeBook_Tx(2,:) = estimated_best_beam_tx;
    prmPAControl.CodeBook_Rx(1,:) = estimated_best_angle_rx;
    prmPAControl.CodeBook_Rx(2,:) = estimated_best_beam_rx;
    
    %% Re do the beam training to compare the result
    Number_Test = 1;
    for nthTest = 1 : Number_Test
        disp('-------------------------------------------------------------')
        fprintf('-----------%d-th collection of received raw signal------------\n',nthTest);
        disp('-------------------------------------------------------------')
        for nthTxBeam = 1:1:prmPAControl.Number_CodeBook_Tx_to_train
            prmPAControl.nthTxBeam = nthTxBeam;
            %clear Control_Phase_Array_External_mex %#ok<UNRCH>
            Control_Phase_Array_External_mex(prmPAControl);
            clear run_usrp_rx_mex %#ok<UNRCH>
            run_usrp_rx_mex(prmQPSKReceiver, prmPAControl);
            fprintf('\n-----------Saving data (%d-th transmitter beaming pattern applied)-----------  \n', nthTxBeam);
            save(['Raw_corruptSignal' num2str(nthTest) '_' num2str(nthTxBeam) '.mat'],'Raw_corruptSignal');
            fprintf('\n-----------Saved data (%d-th transmitter beaming pattern applied)----------- %d \n', nthTxBeam)
        end
    end
    save('prmPAControl.mat');
    cpr_decoder;
    BMP = permute(BER_MER_Power_Overall_Training,[3 2 1]);
    MeasurementsPrototype_Test.Data = BMP(:,:,7);

    %% Show result
    cpr_result.AoA_Estimated = AoA_Estimated;
    cpr_result.AoD_Estimated = AoD_Estimated;
    cpr_result.AoA_Estimated_Sweep = AoA_Estimated_Sweep;
    cpr_result.AoD_Estimated_Sweep = AoD_Estimated_Sweep;
    cpr_result.Training_measurement = MeasurementsPrototype_Train.Data;
    cpr_result.Testing_measurement = MeasurementsPrototype_Test.Data;
    Measurement_train_test = zeros(num_codebook_to_train+2,num_codebook_to_train);
    Measurement_train_test(1:num_codebook_to_train,:) = cpr_result.Training_measurement;
    Measurement_train_test(num_codebook_to_train+1:num_codebook_to_train+2,1:2) = cpr_result.Testing_measurement;
    disp(cpr_result)
    figure;
    bar(Measurement_train_test);
    folder = char(strcat(pwd,'/Experimental_Validation/data/val_cpr_data'));
    mat_name = ['/CPR_experiment_result_' num2str(AoD_True) '_' num2str(AoA_True) '_calM' num2str(calibration_method_id) '_M' num2str(num_codebook_to_train) '.mat'];
    matfile = fullfile(folder,mat_name);
    save(matfile)
end

