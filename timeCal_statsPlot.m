function [] = timeCal_statsPlot(OUT_HS,OUT_LS,fittedData)
% plot the stats and output analysis of the moprhoical aorta project


% Normalise model time
% OUT_time_norm = (OUT_time - min(OUT_time))./(max(OUT_time) - min(OUT_time));
OUT_time_norm = linspace(0,1,size(OUT_HS,1));
%% Time variation of model
% H/S
figure('Visible','off')
timeCal_signalStatsPlot(OUT_time_norm,OUT_HS)
plot(OUT_time_norm,fittedData.H_S,'r','LineWidth',1)
legend('median','mean','+ 1 std','- 1 std','data','Location','best')

% xlim([min(min(output_matrix.time)) max(max(output_matrix.time))])
% xlim([0 1])
% ylim([0 1.5])
grid on
xlabel('$t^*$ [-]')
ylabel('$H/S$ [-]')
GM_printBMP(300,300,'stats_H_S')
GM_printEPS(300,300,'stats_H_S')

% L/S
figure('Visible','off')
timeCal_signalStatsPlot(OUT_time_norm,OUT_LS)
plot(OUT_time_norm,fittedData.L_S,'r','LineWidth',1)
legend('median','mean','+ 1 std','- 1 std','data','Location','best')
% xlim([min(min(output_matrix.time)) max(max(output_matrix.time))])
% ylim([0 30])
grid on
xlabel('$t^*$ [-]')
ylabel('$L/S$ [-]')
GM_printBMP(300,300,'stats_L_S')
GM_printEPS(300,300,'stats_L_S')

%% ================== COmmented on 08/12/21
%% Difference with MRI data
% N = 10;
% diff_signal.H_S = timeCal_modelDataDifference(OUT_time,OUT_HS,fittedData(:,2),'H/S [-]',N);
% GM_printBMP(300,500,'stats_DifferenceSignals_H_S')
% GM_printEPS(300,500,'stats_DifferenceSignals_H_S')
% diff_signal.L_S = timeCal_modelDataDifference(OUT_time,OUT_LS,fittedData(:,3),'L/S [-]',N);
% GM_printBMP(300,500,'stats_DifferenceSignals_L_S')
% GM_printEPS(300,500,'stats_DifferenceSignals_L_S')
% 
% % weights_likelihood_all = [2.5 2.5 2.5 2.5]./10;
% weights_likelihood_HS_LS = [1 9]./10;
% diff_signal.likelihood_HS_LS = (weights_likelihood_HS_LS(1)./diff_signal.H_S.likelihood.^(-1/N) + weights_likelihood_HS_LS(2)./diff_signal.L_S.likelihood.^(-1/N) ).^N;
% % diff_signal.likelihood_all = (weights_likelihood_all(1)./diff_signal.H_S.likelihood.^(-1/N) + weights_likelihood_all(2)./diff_signal.L_S.likelihood.^(-1/N) + weights_likelihood_all(3)./diff_signal.SA.likelihood.^(-1/N) + weights_likelihood_all(4)./diff_signal.Vol.likelihood.^(-1/N)).^N;
% 
% % re-scaling likelihood so that sum Likelihood = 1
% diff_signal.likelihood_HS_LS = (diff_signal.likelihood_HS_LS)./(sum(diff_signal.likelihood_HS_LS));
% % diff_signal.likelihood_all = (diff_signal.likelihood_all)./(sum(diff_signal.likelihood_all));
% % log transform
% % diff_signal.likelihood_HS_LS = log(diff_signal.likelihood_HS_LS);
% % diff_signal.likelihood_all = log(diff_signal.likelihood_all);
% 
% % plot likelihood histogram
% figure('Visible','off')
% subplot(221)
% histogram(log(diff_signal.H_S.likelihood),'FaceColor','r','Normalization','probability')
% % set(gca,'XScale','log')
% title('H/S')
% subplot(222)
% histogram(log(diff_signal.L_S.likelihood),'FaceColor','b','Normalization','probability')
% % set(gca,'XScale','log')
% title('L/S')
% % subplot(423)
% % histogram((diff_signal.SA.likelihood),'FaceColor','g','Normalization','probability')
% % set(gca,'XScale','log')
% % title('SA')
% % subplot(424)
% % histogram((diff_signal.Vol.likelihood),'FaceColor','m','Normalization','probability')
% % set(gca,'XScale','log')
% % title('Vol')
% subplot(2,2,[3 4])
% histogram(log(diff_signal.likelihood_HS_LS),'FaceColor',[.7 .7 .7],'Normalization','probability')
% % set(gca,'XScale','log')
% title('H/S + L/S')
% % subplot(4,2,[7 8])
% % histogram((diff_signal.likelihood_all),'FaceColor',[.7 .7 .7],'Normalization','probability')
% % set(gca,'XScale','log')
% % title('All likelihood')
% GM_printBMP(300,500,'stats_LikelihoodFunc')
% GM_printEPS(300,500,'stats_LikelihoodFunc')
% ===================================================


close all

end
