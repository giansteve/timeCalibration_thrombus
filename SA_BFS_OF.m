clear all
close all
clc

uqlab
set(0,'DefaultTextInterpreter','latex')
set(0,'defaultAxesFontSize',11)
% 861108077
degree = 3; % degree PCE
Ns = round(2.05*factorial(9+degree)/(factorial(9)*factorial(degree))); % number of simulations
%% Input
var_names = {'Dc','Kc','K_{BP}','Ct','BPt','RTt','K_{Cwall}','C_{AP}','BP_{bt}'};
input.Marginals(1).Name = 'Dc.'; 
input.Marginals(1).Type = 'Uniform';
input.Marginals(1).Parameters = [1.00e-10 1.00e-6];

input.Marginals(2).Name = 'Kc.'; 
input.Marginals(2).Type = 'Uniform';
input.Marginals(2).Parameters = [200e2 200e4];

input.Marginals(3).Name = 'kBP.'; 
input.Marginals(3).Type = 'Uniform';
input.Marginals(3).Parameters = [8.000e-11 8.000e-9];

input.Marginals(4).Name = 'Ct.'; 
input.Marginals(4).Type = 'Uniform';
input.Marginals(4).Parameters = [10.00e2 10.00e4];

input.Marginals(5).Name = 'BPt.'; 
input.Marginals(5).Type = 'Uniform';
input.Marginals(5).Parameters = [20.00e2 20.00e4];

input.Marginals(6).Name = 'RTt.'; 
input.Marginals(6).Type = 'Uniform';
input.Marginals(6).Parameters = [1.00e-01 3.00e+00];

input.Marginals(7).Name = 'kCwall.'; 
input.Marginals(7).Type = 'Uniform';
input.Marginals(7).Parameters = [100 10.00e4];

input.Marginals(8).Name = 'cAP.'; 
input.Marginals(8).Type = 'Uniform';
input.Marginals(8).Parameters = [1.5e14/20 4.5e14/20];

input.Marginals(9).Name = 'BPbt.'; 
input.Marginals(9).Type = 'Uniform';
input.Marginals(9).Parameters = [100 250e3];

M = size(input.Marginals,2);
INPUT = uq_createInput(input);
exp_design = uq_getSample(INPUT,Ns,'lhs');

%% Model
model.mFile = 'model_funct';
MODEL = uq_createModel(model);

EVAL_model = uq_evalModel(MODEL,exp_design);
save('SA_prePCE.mat')

OUTPUT.signal = EVAL_model;

% growth_rate
OUTPUT.growth_rate = diff(EVAL_model,1,2);
time_rate = size(OUTPUT.signal,2)-1;

% max output
OUTPUT.max_value = max(OUTPUT.signal,[],2);
figure('Name','Frequency max value Volume thrombus')
histogram(OUTPUT.max_value,round(sqrt(length(OUTPUT.max_value))))
grid on
xlabel('$max(\bar{V}_t)$')
ylabel('Frequency')
GM_printEPS(360/2,360/2,'frqMaxVol');

% mean output
OUTPUT.mean_value = mean(OUTPUT.signal,2);
figure('Name','Frequency mean value Volume thrombus')
histogram(OUTPUT.mean_value,round(sqrt(length(OUTPUT.mean_value))))
grid on
xlabel('$mean(\bar{V}_t)$')
ylabel('Frequency')
GM_printEPS(360/2,360/2,'frqMeanVol');


% max value growth rate
[OUTPUT.max_growthRate,OUTPUT.time_maxGR] = max(OUTPUT.growth_rate,[],2);
figure('Name','Frequency max growth rate thrombus')
histogram(OUTPUT.max_growthRate,round(sqrt(length(OUTPUT.max_growthRate))))
grid on
xlabel('$max(\dot{\bar{V}}_t)$')
ylabel('Frequency')
GM_printEPS(360/2,360/2,'frqMaxGR');


figure('Name','Frequency time instant max growth rate thrombus')
histogram(OUTPUT.time_maxGR,round(sqrt(length(OUTPUT.time_maxGR))))
grid on
xlabel('$time max(\dot{\bar{V}}_t) [s] $')
ylabel('Frequency')
GM_printEPS(360/2,360/2,'frqTimeMaxGR');


% evaluating variances
var_sign = var(OUTPUT.signal,1,1);
var_growth = var(OUTPUT.growth_rate,1,1);
figure('Name','Thrombus volume')
for tt = 1:Ns
    plot( (1:time_rate+1)/2 , OUTPUT.signal(tt,:) )
    hold on
end
xlabel('Time [s]')
ylabel('$\bar{V} [1/s]$')
grid on
xlim([0 50])
GM_printEPS(2*360/3,120,'ThrombusVol')


figure('Name','Variance Volumetric average of thrombus')
plot((1:100)/2,var_sign)
grid on
% title()
xlabel('Time [s]')
ylabel('$\sigma^2[\bar{V}]$ [-]','Interpreter','latex')
GM_printEPS(2*360/3,120,'VarThrVol')


figure('Name','Thrombus Growth Rate')
for tt = 1:Ns
    plot( (1:time_rate)/2 , OUTPUT.growth_rate(tt,:) )
    hold on
%     plot(OUTPUT.time_maxGR(tt,1),OUTPUT.max_growthRate(tt,1),'ro')
%     ylim([0 0.008])
%     pause(1)
%     hold off
end
xlabel('Time [s]')
ylabel('growth rate [1/s]')
grid on
GM_printEPS(2*360/3,120,'GrowthRate')

figure('Name','Variance Growth rate of thrombus')
plot((1:99)/2,var_growth)
grid on
% title()
xlabel('Time [s]')
ylabel('$\sigma^2[\dot{\bar{V}}]$ [-]','Interpreter','latex')
GM_printEPS(2*360/3,120,'VarGR')

%% Metamodel
metamodel.Type = 'Metamodel';
metamodel.MetaType = 'PCE';
metamodel.Display = 'quiet';
metamodel.Degree = degree;
metamodel.Input = INPUT;
metamodel.FullModel = MODEL;
% metamodel.ExpDesign.Sampling = 'lhs';
metamodel.ExpDesign.NSamples = Ns;
metamodel.ExpDesign.X = exp_design;

% time signal
for tt = 1:time_rate+1
    metamodel.ExpDesign.Y = OUTPUT.signal(:,tt);
    PCE_signal.time(tt,1) = uq_createModel(metamodel);
    error_signal(tt,1) = PCE_signal.time(tt,1).Error.LOO;
end
figure('Name','LOO error Y_{PCE}:= Volume thrombus')
plot((1:time_rate+1)/2,error_signal)
grid on
% title()
xlabel('Tim [s]')
ylabel('LOO error [-]')

% growth rate
for tt = 1:time_rate
    metamodel.ExpDesign.Y = OUTPUT.growth_rate(:,tt);
    PCE_growth.time(tt,1) = uq_createModel(metamodel);
    error_growth(tt,1) = PCE_growth.time(tt,1).Error.LOO;
end
figure('Name','LOO error Y_{PCE}:= Thrombus growth rate')
plot((1:time_rate)/2,error_growth)
grid on
% title()
xlabel('Time [s]')
ylabel('LOO error [-]')

% maximum value
metamodel.ExpDesign.Y = OUTPUT.max_value;
PCE_max = uq_createModel(metamodel);
% mean value
metamodel.ExpDesign.Y = OUTPUT.mean_value;
PCE_mean = uq_createModel(metamodel);
% maximum value growth rate
metamodel.ExpDesign.Y = OUTPUT.max_growthRate;
PCE_max_growthR = uq_createModel(metamodel);
% time max GR
metamodel.ExpDesign.Y = OUTPUT.time_maxGR;
PCE_timeGR = uq_createModel(metamodel);

save('SA_postPCE')

%% Sensitivity Analysis

% mean value
[SA_1_mean,SA_TOT_mean,SA_names_mean] = SA_Sobols(PCE_mean,M,var_names);
[Sobol_int2_mean,int_names_mean] = SA_int2(PCE_mean,M,var_names);

% max value
[SA_1_max,SA_TOT_max,SA_names_max] = SA_Sobols(PCE_max,M,var_names);
[Sobol_int2_max,int_names_max] = SA_int2(PCE_max,M,var_names);

% max value growth rate
[SA_1_maxGR,SA_TOT_maxGR,SA_names_maxGR] = SA_Sobols(PCE_max_growthR,M,var_names);
[Sobol_int2_maxGR,int_names_maxGR] = SA_int2(PCE_max_growthR,M,var_names);

% time peak growth rate
[SA_1_timeGR,SA_TOT_timeGR,SA_names_timeGR] = SA_Sobols(PCE_timeGR,M,var_names);
[Sobol_int2_timeGR,int_names_timeGR] = SA_int2(PCE_timeGR,M,var_names);

% Time signal
[S_1_sig,S_TOT_sig] = SA_time(PCE_signal.time,M);
S_1_sig(isnan(S_1_sig))=0;
S_TOT_sig(isnan(S_TOT_sig))=0;

% Growth rate
[S_1_growth,S_TOT_growht] = SA_time(PCE_growth.time,M);
S_1_growth(isnan(S_1_growth))=0;
S_TOT_growht(isnan(S_TOT_growht))=0;

%% plot
figure('Name','Mean thrombus Volume')
subplot(1,2,1)
bar(1:M,[SA_1_mean;SA_TOT_mean]',1)
% hold on
% bar(1:M,SA_1_mean,0.75)
ylim([0 1])
hold off
grid on
ylabel('$S_i, S^t_i$')
legend({'Main Order','Total Order'},'Location','ne')
bb = gca;
bb.XTickLabel = SA_names_mean;
bb.XTickLabelRotation = 90;
subplot(1,2,2)
bar(1:length(Sobol_int2_mean),Sobol_int2_mean,0.75)
grid on
ylabel('$S^2_i$')
bb = gca;
bb.XTickLabel = int_names_mean;
bb.XTickLabelRotation = 90;
GM_printEPS(370,200,'SA_meanVol')

figure('Name','Max thrombus Volume')
subplot(1,2,1)
bar(1:M,[SA_1_max;SA_TOT_max]',1)
% hold on
% bar(1:M,SA_1_max,0.75)
ylim([0 1])
hold off
grid on
ylabel('$S_i, S^t_i$')
legend({'Main Order','Total Order'},'Location','ne')
bb = gca;
bb.XTickLabel = SA_names_max;
bb.XTickLabelRotation = 90;
subplot(1,2,2)
bar(1:length(Sobol_int2_max),Sobol_int2_max,0.75)
grid on
ylabel('$S^2_i$')
bb = gca;
bb.XTickLabel = int_names_max;
bb.XTickLabelRotation = 90;
GM_printEPS(370,200,'SA_maxVol')


figure('Name','Max Growth Rate')
subplot(1,2,1)
bar(1:M,[SA_1_maxGR;SA_TOT_maxGR]',1)
% hold on
% bar(1:M,SA_1_maxGR,0.75)
ylim([0 1])
hold off
grid on
ylabel('$S_i, S^t_i$')
legend({'Main Order','Total Order'},'Location','ne')
bb = gca;
bb.XTickLabel = SA_names_mean;
bb.XTickLabelRotation = 90;
subplot(1,2,2)
bar(1:length(Sobol_int2_maxGR),Sobol_int2_maxGR,0.75)
grid on
ylabel('$S^2_i$')
bb = gca;
bb.XTickLabel = int_names_maxGR;
bb.XTickLabelRotation = 90;
GM_printEPS(370,200,'SA_maxGR')

figure('Name','Time Max Growth Rate')
subplot(1,2,1)
bar(1:M,[SA_1_timeGR;SA_TOT_timeGR]',1)
% hold on
% bar(1:M,SA_1_timeGR,0.75)
ylim([0 1])
hold off
grid on
ylabel('$S_i, S^t_i$')
legend({'Main Order','Total Order'},'Location','ne')
bb = gca;
bb.XTickLabel = SA_names_timeGR;
bb.XTickLabelRotation = 90;
subplot(1,2,2)
bar(1:length(Sobol_int2_timeGR),Sobol_int2_timeGR,0.75)
grid on
ylabel('$S^2_i$')
bb = gca;
bb.XTickLabel = int_names_timeGR;
bb.XTickLabelRotation = 90;
GM_printEPS(370,200,'SA_timeMaxGR')

figure('Name','Main_SA_Vol')
plot((1:time_rate+1)/2,S_1_sig)
% title('Main sensitivity indices - Volume of thrombus')
ylim([0 1])
grid on
xlabel('Time [s]')
ylabel('$S_i$ [-]')
legend('Dc','Kc','K_{BP}','Ct','BPt','RTt','K_{Cwall}','C_{AP}','BP_{bt}','Location','bestoutside')
GM_printEPS(360,170,'SA1_timeVol')


figure('Name','TOT_SA_Vol')
plot((1:time_rate+1)/2,S_TOT_sig)
% title('Total sensitivity indices - Volume of thrombus')
ylim([0 1])
grid on
xlabel('Time [s]')
ylabel('$S^T_i$ [-]')
legend('Dc','Kc','K_{BP}','Ct','BPt','RTt','K_{Cwall}','C_{AP}','BP_{bt}','Location','bestoutside')
GM_printEPS(360,170,'SAT_timeVol')


figure('Name','Main_SA_growth')
plot((1:time_rate)/2,S_1_growth)
% title('Main sensitivity indices - Growth rate of thrombus')
ylim([0 1])
grid on
xlabel('Time [s]')
ylabel('$S_i$ [-]')
legend('Dc','Kc','K_{BP}','Ct','BPt','RTt','K_{Cwall}','C_{AP}','BP_{bt}','Location','bestoutside')
GM_printEPS(360,170,'SA1_timeGR')

figure('Name','TOT_SA_growth')
plot((1:time_rate)/2,S_TOT_growht)
% title('Total sensitivity indices - Growth rate of thrombus')
ylim([0 1])
grid on
xlabel('Time [s]')
ylabel('$S^T_i$ [-]')
legend('Dc','Kc','K_{BP}','Ct','BPt','RTt','K_{Cwall}','C_{AP}','BP_{bt}','Location','bestoutside')
GM_printEPS(360,170,'SAT_timeGR')

save('SA_finalPostProc')






