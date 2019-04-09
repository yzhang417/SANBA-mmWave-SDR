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
% This function visualize the performance metrics that are stored in the
% structure array BERMER.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input arguments:
% prmPAControl: a structure array that groups the related parameters
% on phased arrays including available codebooks, receiving power gain,
% antenna switching on/off settings and etc. With prmPAControl, the RF 
% configuration is set.
% BERMER: a structure array that groups the related evaluation 
% metrics which includes 'BER' (1st-row), 'MER in dB' (4-th row), 
% '95% MER in dB' (5-th row), 'Received signal strength indication (RSSI)' 
% (7-th row). The average operation is performed to BERMER_Non_Count to
% obtain BERMER.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Result_Visulization(prmPAControl, BERMER)
    %% Post-processing and visulization
    switch prmPAControl.Program_ID 
        case 1
            BMP = zeros(prmPAControl.Number_rfvga_to_test, prmPAControl.Number_rxvga_to_test, length(BERMER(:,1)));
            for i = 1:1:length(BERMER(:,1))
                BMP(:,:,i) = reshape(BERMER(i,:), prmPAControl.Number_rxvga_to_test, prmPAControl.Number_rfvga_to_test)';
            end
        case {2, 3, 2.1, 9, 9.1, 9.2}
            BMP = zeros(prmPAControl.Number_antenna_to_calibrate, prmPAControl.Number_Phase_Shifter_State, length(BERMER(:,1)));
            for i = 1:1:length(BERMER(:,1))
                BMP(:,:,i) = reshape(BERMER(i,:), prmPAControl.Number_Phase_Shifter_State, prmPAControl.Number_antenna_to_calibrate)';
            end        
        case {4, 5, 4.1, 8}
            BMP = zeros(prmPAControl.Number_CodeBook_Tx_to_train, 1, length(BERMER(:,1)));
            for i = 1:1:length(BERMER(:,1))
                BMP(:,:,i) = reshape(BERMER(i,:), 1, prmPAControl.Number_CodeBook_Tx_to_train)';
            end   
        case {6, 8.1}
            BMP = permute(BERMER,[3 2 1]);
            %% Save measured data
            MeasurementsPrototype.Data = BMP(:,:,7);
            MeasurementsPrototype.Program_ID = prmPAControl.Program_ID;
            save MeasurementsPrototype MeasurementsPrototype;
    end
    %Turn power into dB unit
    %BMP(:,:,7) = 10*log10(BMP(:,:,7)*1000);
    if prmPAControl.Program_ID == 2 || prmPAControl.Program_ID == 3 || prmPAControl.Program_ID == 2.1
        Best_Phase = zeros(5,prmPAControl.Number_antenna_to_calibrate);
        [~, Best_Phase(1,:)] = min(BMP(:,:,1)');
        [~, Best_Phase(2,:)] = max(BMP(:,:,4)');
        [~, Best_Phase(3,:)] = max(BMP(:,:,5)');
        [~, Best_Phase(4,:)] = max(BMP(:,:,6)');
        [~, Best_Phase(5,:)] = max(BMP(:,:,7)');

        %% Figure 1
        figure
        subplot(4,1,1);
        bar(Best_Phase(1,:));
        xlabel('Index of Phase Shifter')
        ylabel('Index of Best Shift State')
        title('Bit Error Rate (BER) (%)');
        subplot(4,1,2);
        bar(Best_Phase(2,:));
        title('Modulation Error Rate (MER) (dB)');
        xlabel('Index of Phase Shifter')
        ylabel('Index of Best Shift State')
        subplot(4,1,3);
        bar(Best_Phase(4,:));
        title('95% Modulation Error Rate (MER) (dB)');
        xlabel('Index of Phase Shifter')
        ylabel('Index of Best Shift State')
        subplot(4,1,4);
        bar(Best_Phase(5,:));
        Measured_CodeBook = [1 Best_Phase(5,:)]
        title('Received signal strength indication (RSSI)');  
        xlabel('Index of Phase Shifter')
        ylabel('Index of Best Shift State') 

        %% Figure 2
        figure
        Best_Phase_Plot=[Best_Phase(1,:);Best_Phase(2,:);Best_Phase(4,:);Best_Phase(5,:)];
        bar(Best_Phase_Plot')
        title('Calibration of Phase Shifter', 'FontSize', 16); 
        xlabel('Index of Phase Shifter','FontSize', 16)
        ylabel('Index of Best Shift State','FontSize', 16)
        lgd = legend('BER','MER','95% MER','Received signal strength indication (RSSI)');
        lgd.FontSize = 14;  
        
    end
    
    %% Figure 3
    figure
    %Training_angular_range = linspace(-10,10,CodeBookToTrain)
    xticklabels_var = {};
    beam_label = find(prmPAControl.beam_label>0);
    for i=1:1:length(prmPAControl.beam_label)
        xticklabels_var{i} = strcat(num2str(beam_label-50),'\circ');
    end
    
    subplot(4,1,1);
    bar(BMP(:,:,1));
    title('Bit Error Rate (BER) (%)');
    xticks([1:prmPAControl.Number_CodeBook_Tx_to_train])
    xticklabels(xticklabels_var);
    
    subplot(4,1,2); 
    bar(BMP(:,:,4));
    title('Modulation Error Rate (MER) (dB)');
    xticks([1:prmPAControl.Number_CodeBook_Tx_to_train])
    xticklabels(xticklabels_var);
    
    subplot(4,1,3);
    bar(BMP(:,:,6));
    title('95% Modulation Error Rate (MER) (dB)');
    xticks([1:prmPAControl.Number_CodeBook_Tx_to_train])
    xticklabels(xticklabels_var);
    
    subplot(4,1,4); 
    bar(BMP(:,:,7));
    title('Received signal strength indication (RSSI)');  
    xlabel('Directionality of Receving Beam Pattern (Broadside = 0\circ, anticlockwise)', 'Fontsize', 20);
    xticks([1:prmPAControl.Number_CodeBook_Tx_to_train])
    xticklabels(xticklabels_var);
    if prmPAControl.Program_ID == 6 || prmPAControl.Program_ID == 8.1
        xticklabels(xticklabels_var);
        %title(leg_Tx,{'Directionality of';'Transmitting Beam Pattern';'(Broadside = 90\circ, anticlockwise)'})
        %leg_Tx.Title.Visible = 'on';
    end        
end

