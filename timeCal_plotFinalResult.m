function [] = timeCal_plotFinalResult(output_matrix,fittedData)
% plot the probability distribution of the results from the model
% simulations at the end of the simulation and compare them with the MRI
% measurements


figure('Visible','off')
subplot(121) % H/S
histogram(output_matrix.H_S_model(end,:),'FaceColor','r','Normalization','probability')
hold on
xline(fittedData(end,2),'k');
xlim([0.5 1]);
grid on
legend('H/S','MRI','Location','nw')

subplot(122) % l/S
histogram(output_matrix.L_S_model(end,:),'FaceColor','b','Normalization','probability')
hold on
xline(fittedData(end,3),'k');
grid on
legend('L/S','MRI','Location','nw')

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