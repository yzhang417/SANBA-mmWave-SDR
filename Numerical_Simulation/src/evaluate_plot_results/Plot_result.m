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
% Label_Name: X-axis label (character array).
% Simulation_result: a structure array that groups all related simulation.
% result obtained. Please refer to main scripts in folder main_programs for
% its fields.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Plot_result(Label_Name, Simulation_result)
    plot_init;
    
    %% Success rate of recovery
    if Simulation_result.L == 1
        Num_Quantization_Error = Simulation_result.Num_Quantization_Error;
        [v, Active_ind] = find(Simulation_result.Method.State>0);
        Max_Error = 20;
        figure
        j = 1;
        for i = Active_ind
            Num_Quantization_Error_i = Num_Quantization_Error(:,:,1,i); 
            Dist_i = hist(Num_Quantization_Error_i,0:Max_Error)/length(Num_Quantization_Error(:,1,1,1));
            subplot(round(sqrt(Simulation_result.Method.Number)),ceil(sqrt(Simulation_result.Method.Number)),j);
            j = j+1;
            %colormap hsv;
            bar3(Dist_i);
            set(gca,'XTickLabel',Simulation_result.Range);
            set(gca,'YTickLabel',0:Max_Error)
            xlabel(Label_Name);
            ylabel('Number of Quantization Error');
            zlabel('Percentage');
            title('Distribution of Recovery In Grid Level '+Simulation_result.Method.Name_of_method(i));
            view(-50,30);
            camlight left;
        end    
    end
    
    
    %% Error of angle in degree
    figure;
    %-----NonCoherent/PerfectPhase/NoisePhase-----
    [v, Active_ind] = find(Simulation_result.Method.State>0);
    for i = Active_ind
        plot(Simulation_result.Range,Simulation_result.Mean_Evaluation(1,:,2,i),Simulation_result.Method.Color_Line(i)); 
        if Simulation_result.Method.Name_of_method(i) == "PLGAMP" && (contains(Label_Name,'Measurements')||contains(Label_Name,'SNR'))
            for j=1:1:length(Simulation_result.Range)
                if contains(Label_Name,'Measurements')
                    arrow = annotation('arrow');  
                    arrow.Parent = gca;           
                    arrow.X = [Simulation_result.Range(j),Simulation_result.Range(j)]; 
                    arrow.Y = [Simulation_result.Mean_Evaluation(1,j,2,i)-0.05,Simulation_result.Mean_Evaluation(1,j,2,i)-0.4];
                    arrow.LineWidth  = 2;      
                    arrow.HeadWidth  = 7;
                    arrow.HeadLength = 7;
                    arrow.Color = 'red';
                    str = ['$G$=' num2str(Simulation_result.G(j))];
                    prop = text(arrow.X(2)-0.5,arrow.Y(2)-0.1,str,'Interpreter','latex');
                else
                    str = ['$G$=' num2str(Simulation_result.G(j))];
                    prop = text(Simulation_result.Range(j)+0.5,Simulation_result.Mean_Evaluation(1,j,2,i)+0.25,str,'Interpreter','latex');
                end
                prop.Color = 'red';
                prop.FontSize = 16;
            end
        end
        hold on;
    end
    %----Beam sweep----
    plot(Simulation_result.Range,Simulation_result.Mean_Evaluation(1,:,12,i),'k^-');
    % label
    xticks(Simulation_result.Range)
    xlabel(Label_Name,'Interpreter','latex');
    ylabel('Mean Angle Estimation Error (MAEE) in Degree','Interpreter','latex');        
    legend_name = Simulation_result.Method.Name_of_method;
    legend_name = legend_name(Active_ind);    
    legend_name = [legend_name, "Beam Sweeping"];      
    legend(legend_name); grid on;
    if contains(Label_Name,'Searching') 
        set(gca, 'XDir','reverse')
    end

    
    %% NMSE of Array Response Vector
    figure;
    %-----NonCoherent/PerfectPhase/NoisePhase-----
    [v, Active_ind] = find(Simulation_result.Method.State>0);
    for i = Active_ind
        plot(Simulation_result.Range,10*log10(Simulation_result.Mean_Evaluation(1,:,3,i)),Simulation_result.Method.Color_Line(i)); 
        hold on;
    end
    %-----Labels and legends-----
    xlabel(Label_Name,'Interpreter','latex');
    ylabel('NMSE of Array Response Vector (dB)','Interpreter','latex'); 
    legend_name = Simulation_result.Method.Name_of_method;
    legend_name = legend_name(Active_ind);        
    legend(legend_name); grid on;
    if contains(Label_Name,'Searching') 
        set(gca, 'XDir','reverse')
    end


    %% NMSE of Channel Matrix
    figure;
    %-----NonCoherent/PerfectPhase/NoisePhase-----
    [v, Active_ind] = find(Simulation_result.Method.State>0);
    for i = Active_ind
        plot(Simulation_result.Range,10*log10(Simulation_result.Mean_Evaluation(1,:,4,i)),Simulation_result.Method.Color_Line(i)); 
        hold on;
    end
    %-----Labels and legends-----
    xlabel(Label_Name,'Interpreter','latex');
    ylabel('NMSE of Channel Matrix H (dB)','Interpreter','latex');  
    legend_name = Simulation_result.Method.Name_of_method;
    legend_name = legend_name(Active_ind);        
    legend(legend_name); grid on;
    if contains(Label_Name,'Searching') 
        set(gca, 'XDir','reverse')
    end

    
    %% Spectrum Efficiency
    figure;
    %-----Perfect CSI-----
    plot(Simulation_result.Range,Simulation_result.Mean_Evaluation(1,:,6,i),'cd-');hold on;
    %-----NonCoherent/PerfectPhase/NoisePhase-----
    [v, Active_ind] = find(Simulation_result.Method.State>0);
    for i = Active_ind
        plot(Simulation_result.Range,Simulation_result.Mean_Evaluation(1,:,7,i),Simulation_result.Method.Color_Line(i));
        if Simulation_result.Method.Name_of_method(i) == "PLGAMP" && (contains(Label_Name,'Measurements')||contains(Label_Name,'SNR'))
            for j=1:1:length(Simulation_result.Range)
                arrow = annotation('arrow');  
                arrow.Parent = gca;           
                arrow.X = [Simulation_result.Range(j),Simulation_result.Range(j)]; 
                arrow.Y = [Simulation_result.Mean_Evaluation(1,j,7,i)-0.04,Simulation_result.Mean_Evaluation(1,j,7,i)-0.14];
                arrow.LineWidth  = 2;      
                arrow.HeadWidth  = 7;
                arrow.HeadLength = 7;
                arrow.Color = 'red';
                str = ['$G$=' num2str(Simulation_result.G(j))];
                prop = text(arrow.X(2)-0.5,arrow.Y(2)-0.05,str,'Interpreter','latex');
                prop.Color = 'red';
                prop.FontSize = 16;
            end
        end
        hold on;
    end
    %-----Beam Sweep-----
    plot(Simulation_result.Range,Simulation_result.Mean_Evaluation(1,:,10,i),'k^-');hold on;
    %-----Labels and legends-----
    xlabel(Label_Name,'Interpreter','latex');
    xticks(Simulation_result.Range)
    ylabel('Spectrum Efficiency (Bits/s/hz)','Interpreter','latex');
    legend_name = Simulation_result.Method.Name_of_method;
    legend_name = legend_name(Active_ind);
    legend_name = ["With Perfect CSI",legend_name,"Beam Sweeping"];            
    legend(legend_name); grid on;
    if contains(Label_Name,'Searching') 
        set(gca, 'XDir','reverse')
    end 
end

