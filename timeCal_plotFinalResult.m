function [] = timeCal_plotFinalResult(OUT_HS,OUT_LS,fittedData)
% plot the probability distribution of the results from the model
% simulations at the end of the simulation and compare them with the MRI
% measurements

figure('Visible','off')
subplot(121) % H/S
histogram(OUT_HS(end,:),'FaceColor','r','Normalization','probability')
hold on
xline(fittedData.H_S(end),'k','LineWidth',3);
% xlim([0.5 1]);
grid on
legend('$H/S$','MRI','Location','nw')

subplot(122) % L/S
histogram(OUT_LS(end,:),'FaceColor','b','Normalization','probability')
hold on
xline(fittedData.L_S(end),'k','LineWidth',3);
grid on
legend('$L/S$','MRI','Location','ne')

% subplot(223) % SA
% histogram(output_matrix.SA_model(end,:),'FaceColor','g','Normalization','probability')
% hold on
% xline(fittedData(end,4),'k');
% grid on
% legend('SA','MRI','Location','nw')
% 
% subplot(224) % Vol
% histogram(output_matrix.VOL_model(end,:),'FaceColor','m','Normalization','probability')
% hold on
% xline(fittedData(end,5),'k');
% grid on
% legend('Vol','MRI','Location','nw')
GM_printBMP(300,300,'stats_finalTime')
GM_printEPS(300,300,'stats_finalTime')


end