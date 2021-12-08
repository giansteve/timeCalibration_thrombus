function [] = timeCal_plotSA_likelihood(SA,time,M,var_names)
% plot the SA results

time_vec = (time - min(time))./(max(time) - min(time));

%% H/S
SA_plot_time(SA.H_S,time_vec,'H_S','$t^*_\mathrm{mod}$')
subplot(211); legend(var_names,'Interpreter','latex','Location','bestoutside')
GM_printBMP(300,300,'SA_H_S')
GM_printEPS(300,300,'SA_H_S')

%% L/S
SA_plot_time(SA.L_S,time_vec,'L_S','$t^*_\mathrm{mod}$')
subplot(211); legend(var_names,'Interpreter','latex','Location','bestoutside')
GM_printBMP(300,300,'SA_L_S')
GM_printEPS(300,300,'SA_L_S')

% SA
% SA_plot_time(SA.SA,time_vec,'SA','$t^*_\mathrm{mod}$')
% subplot(211); legend(var_names,'Interpreter','latex','Location','bestoutside')
% GM_printBMP(300,300,'SA_SA')
% GM_printEPS(300,300,'SA_SA')

% Vol
% SA_plot_time(SA.Vol,time_vec,'Vol','$t^*_\mathrm{mod}$')
% subplot(211); legend(var_names,'Interpreter','latex','Location','bestoutside')
% GM_printBMP(300,300,'SA_Vol')
% GM_printEPS(300,300,'SA_Vol')

%% LOG Likelihood
% figure('Visible','off')
% b = bar([SA.LOGlikeL.main;SA.LOGlikeL.total]');
% b(1).FaceColor = [.75 .75 .75];
% b(1).BarWidth = 1;
% b(2).BarWidth = 1;
% b(2).FaceColor = [0 0 0];
% hold on
% set(gca,'XTickLabel',var_names)
% grid on
% ylim([0 1])
% GM_printBMP(204*1.3,204,'SA_LOGLikelihood')
% GM_printEPS(204*1.3,204,'SA_LOGLikelihood')

%% LS LOG Likelihood
figure('Visible','off')
b = bar([SA.LOGlikeLS.main;SA.LOGlikeLS.total]');
b(1).FaceColor = [.75 .75 .75];
b(1).BarWidth = 1;
b(2).BarWidth = 1;
b(2).FaceColor = [0 0 0];
hold on
set(gca,'XTickLabel',var_names)
grid on
ylim([0 1])
GM_printBMP(204*1.3,204,'SA_LS_LOGLikelihood')
GM_printEPS(204*1.3,204,'SA_LS_LOGLikelihood')

%%
close all

end