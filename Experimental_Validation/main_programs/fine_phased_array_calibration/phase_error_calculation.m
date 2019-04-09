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
% Script description:
% This script calculates the phase error with two fine calibration methods. 
% These two fine phased array calibration methods result in similar 
% performance and the method 1 is detailed in the above paper.
% Actually, three calibration methods are performed in this script.
% 1. coarse phased array calibration: in this method, the four phase states
%    are considered to be exactly 0, 90, 180, 270. The results are stored
%    in coarse_calibrate_Tx.mat and coarse_calibrate_Rx.mat.
%    For this method, the phase error is stored in coarse_calibrate_Tx.mat 
%    and coarse_calibrate_Rx.mat in folder [Experimental_Validation/data/
%    cal_result].
%    For this method, coarse_calibrate_Tx.mat/coarse_calibrate_Rx.mat only 
%    have four possibilities of values which are four types of error 0, 90, 
%    180 and 270 degrees.
% 2. fine phased array calibration method 1: in this method, the four phase
%    states are considered to be 0+epsilon_n, 90+epsilon_n, 180+epsilon_n,
%    and 270+epsilon_n for the n-th antenna element. Thus, the phase error
%    is antenna dependent but independent of the phase state.
%    For this method, the phase error is stored in Ave_Error_Phase_Tx.mat 
%    and Ave_Error_Phase_Rx.mat in folder [Experimental_Validation/data/
%    cal_result]. For this method, the Ave_Error_Phase_Tx.mat and 
%    Ave_Error_Phase_Rx.mat store the estimated phase shift error epsilon_n 
%    which needs to be added to the nominal phase 0, 90, 180 or 270 to be 
%    the real phase corresponding to the nominal four types 00, 01, 10, 11 
%    phase states.
% 3. fine phased array calibration method 2: in this method, the four phase
%    states are considered to be 0+epsilon_n,1, 90+epsilon_n,2, 
%    180+epsilon_n,3, and 270+epsilon_n,4. This means the phase error is
%    not only antenna dependent but also phase state dependent, which is
%    a more realistic model.
%    For this method, the phase error is stored in Fine_calibrate_Tx.mat 
%    and Fine_calibrate_Rx.mat in folder [Experimental_Validation/data/
%    cal_result]. In particular, the Fine_calibrate_Tx.mat and 
%    Fine_calibrate_Rx.mat store the estimated phase of the nominal four 
%    types 00, 01, 10, 11 phase states.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Initialization
clc; clear all; %close all

%% Add paths
folder_name = ...
[
"/Experimental_Validation";
"/Numerical_Simulation";
];
for i=1:length(folder_name)
    Lib_path = char(strcat(pwd,folder_name(i)));
    addpath(genpath(Lib_path));
end

%% Load data   
Activate_Tx_ID_range = 1:1:12;
Activate_Tx_Phase_ID_range = 1:1:4;
Activate_Rx_ID_range = 1:1:12;
Activate_Rx_Phase_ID_range = 1:1:4;
Element_Gain = zeros(length(Activate_Tx_ID_range),...
                     length(Activate_Tx_Phase_ID_range),...
                     length(Activate_Rx_ID_range),...
                     length(Activate_Rx_Phase_ID_range));   
Ant_To_Be_Calibrated_Ind_range = 2:12;
Ant_To_Be_Calibrated_Phase_range = 1:1:4;
Ant_Ref_Ind_Range = [1];
Ant_Ref_Phase_Range = [1 2 3 4];
Fixed_Side_Antenna_ID_Range = [5];
Fixed_Side_Antenna_Phase_Range = [1];
REV_Gain = zeros(12,4,12,4,12,4);                     
REV_Gain = 1000*REV_Gain;
Element_Gain = 1000*Element_Gain;

prmPAControl.Program_ID = 9.1; 
cal_data_path = char(strcat(pwd,'/Experimental_Validation/data/cal_result/'));
load(char(strcat(cal_data_path,['REV_Gain' num2str(prmPAControl.Program_ID) '.mat'])),'REV_Gain');
load(char(strcat(cal_data_path,'Element_Gain.mat')),'Element_Gain');

%% Plot initializtion
plot_init

%% Fine phase error calibration method 1
Error_Array = zeros(6*2,...
                  length(Ant_To_Be_Calibrated_Ind_range)+1,...
                  length(Ant_Ref_Ind_Range),...
                  length(Ant_Ref_Phase_Range),...
                  length(Fixed_Side_Antenna_ID_Range),...
                  length(Fixed_Side_Antenna_Phase_Range));
Reference = ones(12,length(Ant_Ref_Ind_Range),...
                  length(Ant_Ref_Phase_Range),...
                  length(Fixed_Side_Antenna_ID_Range),...
                  length(Fixed_Side_Antenna_Phase_Range));
            
for Fix_Antenna = 1:length(Fixed_Side_Antenna_ID_Range)
    for Fix_Antenna_Phase = Fixed_Side_Antenna_Phase_Range              
        for Ant_Ref_Ind = 1:length(Ant_Ref_Ind_Range)
            for Ant_Ref_Phase = Ant_Ref_Phase_Range              
                Reference(1,Ant_Ref_Ind,Ant_Ref_Phase,Fix_Antenna,Fix_Antenna_Phase) = Ant_Ref_Phase;
                for Ant_To_Be_Calibrated_Ind = Ant_To_Be_Calibrated_Ind_range
                    n = Ant_To_Be_Calibrated_Ind;
                    m = REV_Gain(n,:,Ant_Ref_Ind,Ant_Ref_Phase,Fixed_Side_Antenna_ID_Range(Fix_Antenna),Fix_Antenna_Phase);
                    [v, closest] = max(m);
                    Reference(n,Ant_Ref_Ind,Ant_Ref_Phase,Fix_Antenna,Fix_Antenna_Phase) = closest;
                    Delta = closest - 1;
                    error_min = -pi/4-Delta*pi/2;
                    error_max = pi/4-Delta*pi/2;
                    initial_guess = (error_min+error_max)/2;      
                    Error_Array(:,n,Ant_Ref_Ind,Ant_Ref_Phase,Fix_Antenna,Fix_Antenna_Phase) = inf;
                    P = [1 1 1 2 2 3];
                    Q = [2 3 4 3 4 4];
                    Remove_Ambiguity.CheckPoint = [mod(closest+1,4) mod(closest+3,4)];
                    Remove_Ambiguity.CheckPoint(find(Remove_Ambiguity.CheckPoint==0)) = 4;
                    if prmPAControl.Program_ID == 9.1
                        Remove_Ambiguity.a_1_1 = Element_Gain(Ant_Ref_Ind_Range(Ant_Ref_Ind),Ant_Ref_Phase,...
                                                              Fixed_Side_Antenna_ID_Range(Fix_Antenna),Fix_Antenna_Phase);
                        Remove_Ambiguity.a_n_p = Element_Gain(n,Remove_Ambiguity.CheckPoint(1),...
                                                              Fixed_Side_Antenna_ID_Range(Fix_Antenna),Fix_Antenna_Phase);
                        Remove_Ambiguity.a_n_q = Element_Gain(n,Remove_Ambiguity.CheckPoint(2),...
                                                              Fixed_Side_Antenna_ID_Range(Fix_Antenna),Fix_Antenna_Phase);  
                    else
                        Remove_Ambiguity.a_1_1 = Element_Gain(Fixed_Side_Antenna_ID_Range(Fix_Antenna),Fix_Antenna_Phase,...
                                                              Ant_Ref_Ind_Range(Ant_Ref_Ind),Ant_Ref_Phase);
                        Remove_Ambiguity.a_n_p = Element_Gain(Fixed_Side_Antenna_ID_Range(Fix_Antenna),Fix_Antenna_Phase,...
                                                              n,Remove_Ambiguity.CheckPoint(1));
                        Remove_Ambiguity.a_n_q = Element_Gain(Fixed_Side_Antenna_ID_Range(Fix_Antenna),Fix_Antenna_Phase,...
                                                              n,Remove_Ambiguity.CheckPoint(2));             
                    end
                    Remove_Ambiguity.mp = m(Remove_Ambiguity.CheckPoint(1));
                    Remove_Ambiguity.mq = m(Remove_Ambiguity.CheckPoint(2));
                    error_ind = 1;
                    for i = 1:length(P)
                        p = P(i);
                        q = Q(i);
                        if prmPAControl.Program_ID == 9.1
                            a_n_minus_1_1 = Element_Gain(Ant_Ref_Ind_Range(Ant_Ref_Ind),Ant_Ref_Phase,...
                                                         Fixed_Side_Antenna_ID_Range(Fix_Antenna),Fix_Antenna_Phase);
                            a_n_p = Element_Gain(n,p,Fixed_Side_Antenna_ID_Range(Fix_Antenna),Fix_Antenna_Phase);
                            a_n_q = Element_Gain(n,q,Fixed_Side_Antenna_ID_Range(Fix_Antenna),Fix_Antenna_Phase);
                            if a_n_minus_1_1 == 0 || a_n_p==0 || a_n_q ==0
                                fprintf('Issue')
                            end
                        else
                            a_n_minus_1_1 = Element_Gain(Fixed_Side_Antenna_ID_Range(Fix_Antenna),Fix_Antenna_Phase,...
                                                         Ant_Ref_Ind_Range(Ant_Ref_Ind),Ant_Ref_Phase);
                            a_n_p = Element_Gain(Fixed_Side_Antenna_ID_Range(Fix_Antenna),Fix_Antenna_Phase,n,p);
                            a_n_q = Element_Gain(Fixed_Side_Antenna_ID_Range(Fix_Antenna),Fix_Antenna_Phase,n,q);           
                            if a_n_minus_1_1 == 0 || a_n_p==0 || a_n_q ==0
                                fprintf('Issue')
                            end
                        end
                        Cpq = m(p)^2/m(q)^2 * (a_n_minus_1_1^2+a_n_q^2) - (a_n_minus_1_1^2+a_n_p^2);
                        Apq = 2*a_n_minus_1_1*a_n_p*cos((p-1)*pi/2) - m(p)^2/m(q)^2*(2*a_n_minus_1_1*a_n_q)*cos((q-1)*pi/2);
                        Bpq = 2*a_n_minus_1_1*a_n_p*sin((p-1)*pi/2) - m(p)^2/m(q)^2*(2*a_n_minus_1_1*a_n_q)*sin((q-1)*pi/2);
                        a = Apq^2 + Bpq^2;
                        b = 2*Bpq*Cpq;
                        c = Cpq^2 - Apq^2;
                        D = b^2-4*a*c;
                        if D<=1e-10
                            if D<=0
                                %fprintf("\nOne solutions\n");
                                sol = -b/2*a;
                                if abs(sol) <= 1
                                    [error_array, L] = Error_Transformation(sol, error_min, error_max, Remove_Ambiguity);
                                    for l=1:L
                                        Error_Array(error_ind,n,Ant_Ref_Ind,Ant_Ref_Phase,Fix_Antenna,Fix_Antenna_Phase)...
                                            = error_array(l)+pi/2*(Ant_Ref_Phase-1);
                                        error_ind = error_ind + 1;
                                    end
                                end
                            end
                        else
                            %fprintf("\nTwo solutions\n");
                            sol1 = (-b+sqrt(D))/(2*a);
                            sol2 = (-b-sqrt(D))/(2*a);
                            if abs(sol1) <= 1
                                [error_array, L] = Error_Transformation(sol1, error_min, error_max, Remove_Ambiguity);
                                for l=1:L
                                    Error_Array(error_ind,n,Ant_Ref_Ind,Ant_Ref_Phase,Fix_Antenna,Fix_Antenna_Phase) ...
                                        = error_array(l)+pi/2*(Ant_Ref_Phase-1);
                                    error_ind = error_ind + 1;
                                end
                            end
                            if abs(sol2) <= 1
                                [error_array, L] = Error_Transformation(sol2, error_min, error_max, Remove_Ambiguity);
                                for l=1:L
                                    Error_Array(error_ind,n,Ant_Ref_Ind,Ant_Ref_Phase,Fix_Antenna,Fix_Antenna_Phase) ...
                                        = error_array(l)+pi/2*(Ant_Ref_Phase-1);
                                    error_ind = error_ind + 1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
Error_Array_Deg = rad2deg(Error_Array);

%% Calculate and plot reference (coarse phase error calibration) with method 1
for Fix_Antenna = 1:length(Fixed_Side_Antenna_ID_Range)
    for Fix_Antenna_Phase = Fixed_Side_Antenna_Phase_Range              
        for Ant_Ref_Ind = 1:length(Ant_Ref_Ind_Range)
            for Ant_Ref_Phase = Ant_Ref_Phase_Range 
                Reference_Temp = Reference(:,Ant_Ref_Ind,Ant_Ref_Phase,Fix_Antenna,Fix_Antenna_Phase);
                Delta = Reference_Temp(1) - 1;
                for n=1:12
                    Reference_Temp(n) = Reference_Temp(n) - Delta;
                    if Reference_Temp(n) <= 0
                        Reference_Temp(n) = Reference_Temp(n) + 4;
                    end
                end
                if Ant_Ref_Phase == 1
                    coarse_calibrate = Reference_Temp;
                    if prmPAControl.Program_ID == 9.1
                        save(char(strcat(cal_data_path,['coarse_calibrate_Tx.mat'])),'coarse_calibrate');
                    else
                        save(char(strcat(cal_data_path,['coarse_calibrate_Rx.mat'])),'coarse_calibrate');
                    end
                    figure
                    bar(Reference_Temp)
                end
            end
        end
    end
end


%% Plot fine phase error calibration result method 1
Error_Phase = zeros(4,12);
figure_ID = randi([1 10000],1,1);
for Fix_Antenna = 1:length(Fixed_Side_Antenna_ID_Range)
    for Fix_Antenna_Phase = Fixed_Side_Antenna_Phase_Range              
        for Ant_Ref_Ind = 1:length(Ant_Ref_Ind_Range)
            for Ant_Ref_Phase = Ant_Ref_Phase_Range  
                for n=1:1:12
                    Error_Array_Deg_Temp = Error_Array_Deg(:,n,Ant_Ref_Ind,Ant_Ref_Phase,Fix_Antenna,Fix_Antenna_Phase);
                    Error_Array_Deg_Temp = Error_Array_Deg_Temp(find(Error_Array_Deg_Temp<inf));
                    Error_Array_Deg_Temp = Error_Array_Deg_Temp(find(Error_Array_Deg_Temp>-inf));
                    figure(figure_ID)
                    subplot(3,5,n)
                    if prmPAControl.Program_ID == 9.1 
                        title(['Tx antenna index: ' num2str(n)])
                    else
                        title(['Rx antenna index: ' num2str(n)])
                    end
                    for j=1:length(Error_Array_Deg_Temp)
                        if Ant_Ref_Phase == 1
                            set(0,'DefaultLineMarkerSize',8);
                            set(0,'DefaultLineLineWidth', 2);
                            h0 = polar(deg2rad(Error_Array_Deg_Temp(j))*ones(1,6),0:0.2:1,'kx-'); hold on;
                        end
                    end
                    switch Ant_Ref_Phase 
                        case 1
                            set(0,'DefaultLineLineWidth', 6);
                            h1 = polar(mean(deg2rad(Error_Array_Deg_Temp))*ones(1,11),0:0.1:1,'r-'); hold on;
                            h2 = polar((mean(deg2rad(Error_Array_Deg_Temp))+pi/2)*ones(1,11),0:0.1:1,'y'); hold on;
                            h3 = polar((mean(deg2rad(Error_Array_Deg_Temp))+2*pi/2)*ones(1,11),0:0.1:1,'b'); hold on;
                            h4 = polar((mean(deg2rad(Error_Array_Deg_Temp))+3*pi/2)*ones(1,11),0:0.1:1,'g'); hold on;
                        case 2
                            %polar(mean(deg2rad(Error_Array_Deg_Temp))*ones(1,100),1:100,'y*'); hold on;
                        case 3
                            %polar(mean(deg2rad(Error_Array_Deg_Temp))*ones(1,100),1:100,'b*'); hold on;
                        case 4
                            %polar(mean(deg2rad(Error_Array_Deg_Temp))*ones(1,100),1:100,'g*'); hold on;
                    end
                    Error_Phase(Ant_Ref_Phase,n) = deg2rad(mean(Error_Array_Deg_Temp));
                end
            end
        end
    end
end
%Ave_Error_Phase = wrapToPi(mean(wrapTo2Pi(Error_Phase)));
Ave_Error_Phase = wrapToPi((wrapTo2Pi(Error_Phase(1,:))));
if prmPAControl.Program_ID == 9.1
    save(char(strcat(cal_data_path,['Ave_Error_Phase_Tx.mat'])),'Ave_Error_Phase');
else
    save(char(strcat(cal_data_path,['Ave_Error_Phase_Rx.mat'])),'Ave_Error_Phase');
end
figure(figure_ID)
legend([h0 h1 h2 h3 h4],'Possible values of phase state 1',...
                        'Averaged value of phase state 1',...
                        'Averaged value of phase state 2',...
                        'Averaged value of phase state 3',...
                        'Averaged value of phase state 4');


%% Fine phase error calibration method 2
Error_Array2 = zeros(1,...
                  length(Ant_To_Be_Calibrated_Ind_range)+1,...
                  length(Ant_To_Be_Calibrated_Phase_range),...
                  length(Ant_Ref_Ind_Range),...
                  length(Ant_Ref_Phase_Range),...
                  length(Fixed_Side_Antenna_ID_Range),...
                  length(Fixed_Side_Antenna_Phase_Range));
              
for Fix_Antenna = 1:length(Fixed_Side_Antenna_ID_Range)
    for Fix_Antenna_Phase = Fixed_Side_Antenna_Phase_Range              
        for Ant_Ref_Ind = 1:length(Ant_Ref_Ind_Range)
            for Ant_Ref_Phase = Ant_Ref_Phase_Range              
                for Ant_To_Be_Calibrated_Ind = Ant_To_Be_Calibrated_Ind_range              
                    n = Ant_To_Be_Calibrated_Ind;
                    m = REV_Gain(n,:,Ant_Ref_Ind,Ant_Ref_Phase,Fixed_Side_Antenna_ID_Range(Fix_Antenna),Fix_Antenna_Phase);
                    [v, closest] = max(m);              
                    for j = Ant_To_Be_Calibrated_Phase_range
                        error_ind = 1;
                        if prmPAControl.Program_ID == 9.1
                            a_n_minus_1_1 = Element_Gain(Ant_Ref_Ind_Range(Ant_Ref_Ind),Ant_Ref_Phase,...
                                                         Fixed_Side_Antenna_ID_Range(Fix_Antenna),Fix_Antenna_Phase); 
                            a_n_p = Element_Gain(n,j,Fixed_Side_Antenna_ID_Range(Fix_Antenna),Fix_Antenna_Phase); 
                            if a_n_minus_1_1 == 0 || a_n_p==0
                                fprintf('Issue')
                            end
                        else
                            a_n_minus_1_1 = Element_Gain(Fixed_Side_Antenna_ID_Range(Fix_Antenna),Fix_Antenna_Phase,...
                                                         Ant_Ref_Ind_Range(Ant_Ref_Ind),Ant_Ref_Phase); 
                            a_n_p = Element_Gain(Fixed_Side_Antenna_ID_Range(Fix_Antenna),Fix_Antenna_Phase,n,j);     
                            if a_n_minus_1_1 == 0 || a_n_p==0
                                fprintf('Issue')
                            end
                        end
                        sol = (m(j)^2-(a_n_minus_1_1.^2+a_n_p.^2)) / (2*a_n_minus_1_1*a_n_p);
                        if abs(sol) > 1
                            sol = sign(sol);
                        end
                        [error_array2, L] = Error_Transformation2(sol, closest, j);
                        if isempty(error_array2)
                            fprintf('\nEmpty error_array2\n')
                        end
                        for l=1:L
                            Error_Array2(error_ind,n,j,Ant_Ref_Ind,Ant_Ref_Phase,Fix_Antenna,Fix_Antenna_Phase) = error_array2(l);
                            error_ind = error_ind + 1;
                        end
                    end
                end
            end
        end
    end
end
Error_Array_Deg2 = rad2deg(Error_Array2);

%% Calculate and plot reference (coarse phase error calibration) with method 2
Error_Array_Deg2_Ave = zeros(12,4);
for n=2:12
    Error_Array_Deg_Temp = squeeze(Error_Array_Deg2(1,n,:,1,:))';
    for i=0:1:3
        if i==0
            temp = diag(Error_Array_Deg_Temp,i);
            Error_Array_Deg2_Ave(n,i+1) = mean(temp);
            Error_Array_Deg2_Ave(n,i+1) = temp(1);
        else
            temp = [diag(Error_Array_Deg_Temp,i);diag(Error_Array_Deg_Temp,i-4)];
            Error_Array_Deg2_Ave(n,i+1) = mean(temp);
            Error_Array_Deg2_Ave(n,i+1) = temp(1);
        end
    end
end

[v p] = min(Error_Array_Deg2_Ave');
figure
bar(p);

%% Remove ambiguity for fine phase error calibration method 2
A = Error_Array_Deg2_Ave';
Sign_Matrix_Overall = zeros(4,12,11);
for ref_antenna = 2:12
    %Identify sign of the first two state
    j = 2;
    ind = 0;
    Delta = zeros(11,4);
    for sign1 = [+1 -1]
        for sign2 = [+1 -1]
            ind = ind + 1;
            Delta(ref_antenna-1,ind) = wrapTo360(A(j,ref_antenna)*sign2 - A(j-1,ref_antenna)*sign1);
            for n=2:12
                Delta1 = A(j,n)*1 - A(j-1,n)*1;
                Delta2 = A(j,n)*1 - A(j-1,n)*(-1);
                Delta3 = A(j,n)*(-1) - A(j-1,n)*1;
                Delta4 = A(j,n)*(-1) - A(j-1,n)*(-1);
                DeltaSet = wrapTo360([Delta1 Delta2 Delta3 Delta4]);
                [v p] = min(abs(DeltaSet-Delta(ref_antenna-1,ind)));
                Delta(n-1,ind) = DeltaSet(p);
            end
        end
    end
    Var_Delta = var(Delta);
    [v p] = min(Var_Delta);
    if (abs(Var_Delta(1)-Var_Delta(4))<1e-5) || (abs(Var_Delta(2)==Var_Delta(3))<1e-5)
        if p == 1 || p == 4 %[1 1] %[-1 -1]
            p = inf;
            if wrapTo360(A(j,ref_antenna)*1 - A(j-1,ref_antenna)*1) > 0 && ...
                    wrapTo360(A(j,ref_antenna)*1 - A(j-1,ref_antenna)*1)<180
                p = 1;
            end
            if wrapTo360(A(j,ref_antenna)*(-1) - A(j-1,ref_antenna)*(-1)) > 0 && ...
                    wrapTo360(A(j,ref_antenna)*(-1) - A(j-1,ref_antenna)*(-1))<180 
                p = 4;
            end
            if p == inf
                fprintf('Issue1');
            end           
        else %[1 -1] %[-1 1]
            p = inf;
            if wrapTo360(A(j,ref_antenna)*(-1) - A(j-1,ref_antenna)*1) > 0 && ...
                    wrapTo360(A(j,ref_antenna)*(-1) - A(j-1,ref_antenna)*1)<180
                p = 2;
            end
            if wrapTo360(A(j,ref_antenna)*(1) - A(j-1,ref_antenna)*(-1)) > 0 && ...
                    wrapTo360(A(j,ref_antenna)*(1) - A(j-1,ref_antenna)*(-1))<180
                p = 3;
            end
            if p == inf
                fprintf('Issue2');
            end        
        end
    else
        fprintf('Issue3');
    end

    %Identify sign of the first two states
    SIGN = zeros(1,4);
    switch p
        case 1
            SIGN(1) = 1; SIGN(2) = 1;
        case 2
            SIGN(1) = 1; SIGN(2) = -1;
        case 3
            SIGN(1) = -1; SIGN(2) = 1;
        case 4
            SIGN(1) = -1; SIGN(2) = -1;
    end

    %Identify sign of the last two states
    for j = 3:4
        ind = 0;
        Delta = zeros(11,4);
        for sign1 = [+1 -1]
            for sign2 = [+1 -1]
                ind = ind + 1;
                Delta(ref_antenna-1,ind) = wrapTo360(A(j,ref_antenna)*sign2 - A(j-1,ref_antenna)*sign1);
                for n=2:12
                    Delta1 = A(j,n)*1 - A(j-1,n)*1;
                    Delta2 = A(j,n)*1 - A(j-1,n)*(-1);
                    Delta3 = A(j,n)*(-1) - A(j-1,n)*1;
                    Delta4 = A(j,n)*(-1) - A(j-1,n)*(-1);
                    DeltaSet = wrapTo360([Delta1 Delta2 Delta3 Delta4]);
                    [v p] = min(abs(DeltaSet-Delta(ref_antenna-1,ind)));
                    Delta(n-1,ind) = DeltaSet(p);
                end
            end
        end
        Var_Delta = var(Delta);
        [v p] = min(Var_Delta);
        if (abs(Var_Delta(1)-Var_Delta(4))<1e-5) || (abs(Var_Delta(2)==Var_Delta(3))<1e-5)
            ;
        else
            fprintf('Issue4')
        end
        if (p == 1 || p == 4) 
            SIGN(j) = SIGN(j-1);
        else
            SIGN(j) = -SIGN(j-1);
        end
    end
    SIGN = SIGN;

    %
    % Calulate All Sign
    Ref_SIGN = SIGN;
    Sign_Matrix = zeros(4,12);
    % First two states
    j = 2;
    Delta = wrapTo360(A(j,ref_antenna)*Ref_SIGN(j) - A(j-1,ref_antenna)*Ref_SIGN(j-1));
    for n=2:12
        Delta1 = A(j,n)*1 - A(j-1,n)*1;
        Delta2 = A(j,n)*1 - A(j-1,n)*(-1);
        Delta3 = A(j,n)*(-1) - A(j-1,n)*1;
        Delta4 = A(j,n)*(-1) - A(j-1,n)*(-1);
        DeltaSet = wrapTo360([Delta1 Delta2 Delta3 Delta4]);
        [v p] = min(abs(DeltaSet-Delta));   
        switch p
            case 1
                if wrapTo360(A(j,n)*(1) - A(j-1,n)*1) > 0 && ...
                    wrapTo360(A(j,n)*(1) - A(j-1,n)*1) < 180 
                    Sign_Matrix(j,n) = 1;
                    Sign_Matrix(j-1,n) = 1;
                else
                    fprintf('Issue')
                end
            case 2
                if wrapTo360(A(j,n)*(1) - A(j-1,n)*(-1)) > 0 && ...
                    wrapTo360(A(j,n)*(1) - A(j-1,n)*(-1)) < 180 
                    Sign_Matrix(j,n) = 1;
                    Sign_Matrix(j-1,n) = -1;
                else
                    fprintf('Issue')
                end            
            case 3
                if wrapTo360(A(j,n)*(-1) - A(j-1,n)*1) > 0 && ...
                    wrapTo360(A(j,n)*(-1) - A(j-1,n)*1) < 180 
                    Sign_Matrix(j,n) = -1;
                    Sign_Matrix(j-1,n) = 1;
                else
                    fprintf('Issue')
                end            
            case 4
                if wrapTo360(A(j,n)*(-1) - A(j-1,n)*(-1)) > 0 && ...
                    wrapTo360(A(j,n)*(-1) - A(j-1,n)*(-1)) < 180 
                    Sign_Matrix(j,n) = -1;
                    Sign_Matrix(j-1,n) = -1;
                else
                    fprintf('Issue')
                end
        end
    end
    % Last two states
    for j=3:4
        Delta = wrapTo360(A(j,ref_antenna)*Ref_SIGN(j) - A(j-1,ref_antenna)*Ref_SIGN(j-1));
        for n=2:12
            Delta1 = A(j,n)*1 - A(j-1,n)*Sign_Matrix(j-1,n);
            Delta2 = A(j,n)*(-1) - A(j-1,n)*Sign_Matrix(j-1,n);
            DeltaSet = wrapTo360([Delta1 Delta2]);
            [v p] = min(abs(DeltaSet-Delta));   
            switch p
                case 1
                    if wrapTo360(A(j,n)*(1) - A(j-1,n)*Sign_Matrix(j-1,n)) > 0 && ...
                        wrapTo360(A(j,n)*(1) - A(j-1,n)*Sign_Matrix(j-1,n)) < 180 
                        Sign_Matrix(j,n) = 1;
                    else
                        fprintf('Issue5 %d\n',n);
                        Sign_Matrix(j,n) = 0;
                    end
                case 2
                    if wrapTo360(A(j,n)*(-1) - A(j-1,n)*Sign_Matrix(j-1,n)) > 0 && ...
                        wrapTo360(A(j,n)*(-1) - A(j-1,n)*Sign_Matrix(j-1,n)) < 180 
                        Sign_Matrix(j,n) = -1;
                    else
                        fprintf('Issue5 %d\n',n);
                        Sign_Matrix(j,n) = 0;
                    end            
            end
        end
    end
    Sign_Matrix_Overall(:,:,ref_antenna-1) = Sign_Matrix;  
end

Sign_Matrix_Ave = Sign_Matrix;
Sign_Matrix_Ave(:,:) = 0;
for i = 1:4
    for j=1:12
        Sign_Matrix_Ave(i,j) = mode(squeeze(Sign_Matrix_Overall(i,j,:)));
    end
end

%% Plot fine phase error calibration result method 2
Error_Deg = A.*Sign_Matrix_Ave;
Error_Rad = wrapTo2Pi(deg2rad(Error_Deg));
% figure
% for n=1:12
%     subplot(3,4,n)
%     polar(Error_Rad(1,n)*ones(1,10),1:10,'r*'); hold on;
%     polar(Error_Rad(2,n)*ones(1,10),1:10,'y*'); hold on;
%     polar(Error_Rad(3,n)*ones(1,10),1:10,'b*'); hold on;
%     polar(Error_Rad(4,n)*ones(1,10),1:10,'g*'); hold on;
%     title(['Antenna Index: ' num2str(n)])
% end

Error_Rad2 = wrapTo2Pi(Error_Rad - repmat(Error_Rad(1,:),4,1));
% figure
% for n=1:12
%     subplot(3,4,n)
%     polar(Error_Rad2(1,n)*ones(1,10),1:10,'r*'); hold on;
%     polar(Error_Rad2(2,n)*ones(1,10),1:10,'y*'); hold on;
%     polar(Error_Rad2(3,n)*ones(1,10),1:10,'b*'); hold on;
%     polar(Error_Rad2(4,n)*ones(1,10),1:10,'g*'); hold on;
%     title(['Antenna Index: ' num2str(n)]);
% end
Fine_calibrate = Error_Rad;
Fine_calibrate(:,1) = mean(Error_Rad2(:,2:end)');
Fine_calibrate = wrapToPi(Fine_calibrate);
if prmPAControl.Program_ID == 9.1
    save(char(strcat(cal_data_path,['Fine_calibrate_Tx.mat'])),'Fine_calibrate');
else
    save(char(strcat(cal_data_path,['Fine_calibrate_Rx.mat'])),'Fine_calibrate');
end

figure_ID2 = randi([1 10000],1,1);
figure(figure_ID2);
for n=1:12
    subplot(3,5,n)
    h1 = polar(Fine_calibrate(1,n)*ones(1,11),0:0.1:1,'r-'); hold on;
    h2 = polar(Fine_calibrate(2,n)*ones(1,11),0:0.1:1,'y-'); hold on;
    h3 = polar(Fine_calibrate(3,n)*ones(1,11),0:0.1:1,'b-'); hold on;
    h4 = polar(Fine_calibrate(4,n)*ones(1,11),0:0.1:1,'g-'); hold on;   
    if prmPAControl.Program_ID == 9.1 
        title(['Tx antenna index: ' num2str(n)])
    else
        title(['Rx antenna index: ' num2str(n)])
    end
end
figure(figure_ID2)
legend([h1 h2 h3 h4],'Averaged value of phase state 1',...
                     'Averaged value of phase state 2',...
                     'Averaged value of phase state 3',...
                     'Averaged value of phase state 4');
                 
                 