function [diff_signal] = timeCal_statsPlot(stats_prjct,output_matrix,fittedData)
% plot the stats and output analysis of the moprhoical aorta project

% M = size(exp_design,2);

time_vec = linspace(min(min(output_matrix.Time_model)),max(max(output_matrix.Time_model)),size(output_matrix.Time_model,1));
% Normalise model time
time_vec_norm = (time_vec - min(time_vec))./(max(time_vec) - min(time_vec));

%% Time variation of model
% H/S
timeCal_signalStatsPlot(time_vec_norm,output_matrix.H_S_model,stats_prjct.H_S)
% xlim([min(min(output_matrix.time)) max(max(output_matrix.time))])
% xlim([0 1])
% ylim([0 1.5])
grid on
xlabel('$t^*_\mathrm{mod}$ [-]')
ylabel('$H/S$ [-]')
GM_printBMP(300,300,'stats_H_S')
GM_printEPS(300,300,'stats_H_S')

% L/S
timeCal_signalStatsPlot(time_vec_norm,output_matrix.L_S_model,stats_prjct.L_S)
% xlim([min(min(output_matrix.time)) max(max(output_matrix.time))])
% ylim([0 30])
grid on
xlabel('$t^*_\mathrm{mod}$ [-]')
ylabel('$L/S$ [-]')
GM_printBMP(300,300,'stats_L_S')
GM_printEPS(300,300,'stats_L_S')

% Surface area
% timeCal_signalStatsPlot(time_vec_norm,output_matrix.SA_model,stats_prjct.SA)
% % xlim([min(min(output_matrix.time)) max(max(output_matrix.time))])
% % ylim([0 1.5])
% grid on
% xlabel('$t^*_\mathrm{mod}$ [-]')
% ylabel('$A$ [$m^2$]')
% GM_printBMP(300,300,'stats_SA')
% GM_printEPS(300,300,'stats_SA')

% Volume
% timeCal_signalStatsPlot(time_vec_norm,output_matrix.VOL_model,stats_prjct.Vol)
% % xlim([min(min(output_matrix.time)) max(max(output_matrix.time))])
% % ylim([0 10e-5])
% grid on
% xlabel('$t^*_\mathrm{mod}$ [-]')
% ylabel('$V$ [$m^3$]')
% GM_printBMP(300,300,'stats_Vol')
% GM_printEPS(300,300,'stats_Vol')

%% Test to delete
% figure
% histogram(output_matrix.H_S_model(end,:))
% hold on
% histogram(stats_prjct.H_S.signal(end,:))
% 
% figure
% histogram(output_matrix.L_S_model(end,:))
% hold on
% histogram(stats_prjct.L_S.signal(end,:))
%% Difference with MRI data
N = 10;
diff_signal.H_S = timeCal_modelDataDifference(time_vec,output_matrix.H_S_model,fittedData(:,2),'H/S [-]',N);
GM_printBMP(300,500,'stats_DifferenceSignals_H_S')
GM_printEPS(300,500,'stats_DifferenceSignals_H_S')
diff_signal.L_S = timeCal_modelDataDifference(time_vec,output_matrix.L_S_model,fittedData(:,3),'L/S [-]',N);
GM_printBMP(300,500,'stats_DifferenceSignals_L_S')
GM_printEPS(300,500,'stats_DifferenceSignals_L_S')

% weights_likelihood_all = [2.5 2.5 2.5 2.5]./10;
weights_likelihood_HS_LS = [1 9]./10;
diff_signal.likelihood_HS_LS = (weights_likelihood_HS_LS(1)./diff_signal.H_S.likelihood.^(-1/N) + weights_likelihood_HS_LS(2)./diff_signal.L_S.likelihood.^(-1/N) ).^N;
% diff_signal.likelihood_all = (weights_likelihood_all(1)./diff_signal.H_S.likelihood.^(-1/N) + weights_likelihood_all(2)./diff_signal.L_S.likelihood.^(-1/N) + weights_likelihood_all(3)./diff_signal.SA.likelihood.^(-1/N) + weights_likelihood_all(4)./diff_signal.Vol.likelihood.^(-1/N)).^N;

% re-scaling likelihood so that sum Likelihood = 1
diff_signal.likelihood_HS_LS = (diff_signal.likelihood_HS_LS)./(sum(diff_signal.likelihood_HS_LS));
% diff_signal.likelihood_all = (diff_signal.likelihood_all)./(sum(diff_signal.likelihood_all));
% log transform
% diff_signal.likelihood_HS_LS = log(diff_signal.likelihood_HS_LS);
% diff_signal.likelihood_all = log(diff_signal.likelihood_all);

% plot likelihood histogram
figure('Visible','off')
subplot(221)
histogram(log(diff_signal.H_S.likelihood),'FaceColor','r','Normalization','probability')
% set(gca,'XScale','log')
title('H/S')
subplot(222)
histogram(log(diff_signal.L_S.likelihood),'FaceColor','b','Normalization','probability')
% set(gca,'XScale','log')
title('L/S')
% subplot(423)
% histogram((diff_signal.SA.likelihood),'FaceColor','g','Normalization','probability')
% set(gca,'XScale','log')
% title('SA')
% subplot(424)
% histogram((diff_signal.Vol.likelihood),'FaceColor','m','Normalization','probability')
% set(gca,'XScale','log')
% title('Vol')
subplot(2,2,[3 4])
histogram(log(diff_signal.likelihood_HS_LS),'FaceColor',[.7 .7 .7],'Normalization','probability')
% set(gca,'XScale','log')
title('H/S + L/S')
% subplot(4,2,[7 8])
% histogram((diff_signal.likelihood_all),'FaceColor',[.7 .7 .7],'Normalization','probability')
% set(gca,'XScale','log')
% title('All likelihood')
GM_printBMP(300,500,'stats_LikelihoodFunc')
GM_printEPS(300,500,'stats_LikelihoodFunc')

% % plot SSE error
% figure()
% subplot(421)
% histogram(diff_signal.H_S.indicatorSSE,'FaceColor','r','Normalization','probability')
% title('H/S')
% subplot(422)
% histogram(diff_signal.L_S.indicatorSSE,'FaceColor','b','Normalization','probability')
% title('L/S')
% subplot(423)
% histogram(diff_signal.SA.indicatorSSE,'FaceColor','g','Normalization','probability')
% title('SA')
% subplot(424)
% histogram(diff_signal.Vol.indicatorSSE,'FaceColor','m','Normalization','probability')
% title('Vol')
% subplot(4,2,[5 6])
% histogram(diff_signal.totalSSE,'FaceColor',[.7 .7 .7],'Normalization','probability')
% title('H/S + L/S')
% subplot(4,2,[7 8])
% histogram(diff_signal.totalSSE_all,'FaceColor',[.7 .7 .7],'Normalization','probability')
% title('All SSE')
% GM_printBMP(300,500,'stats_SSEerror')
% GM_printEPS(300,500,'stats_SSEerror')

close all

end
