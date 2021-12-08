function [] = timeCal_plotSA(SA,time,M,var_names)
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

%% SA
SA_plot_time(SA.SA,time_vec,'SA','$t^*_\mathrm{mod}$')
subplot(211); legend(var_names,'Interpreter','latex','Location','bestoutside')
GM_printBMP(300,300,'SA_SA')
GM_printEPS(300,300,'SA_SA')

%% Vol
SA_plot_time(SA.Vol,time_vec,'Vol','$t^*_\mathrm{mod}$')
subplot(211); legend(var_names,'Interpreter','latex','Location','bestoutside')
GM_printBMP(300,300,'SA_Vol')
GM_printEPS(300,300,'SA_Vol')

% %% Diff signal H/S
% SA_plot(SA.diff_indicat.H_S,'diffSignal',0)
% GM_printBMP(300,300,'SA_diffSignal_HS')
% GM_printEPS(300,300,'SA_diffSignal_HS')
% 
% %% Diff signal L/S
% SA_plot(SA.diff_indicat.L_S,'diffSignal',1)
% GM_printBMP(300,300,'SA_diffSignal_LS')
% GM_printEPS(300,300,'SA_diffSignal_LS')
% 
% %% Diff signal SA
% SA_plot(SA.diff_indicat.SA,'diffSignal',1)
% GM_printBMP(300,300,'SA_diffSignal_SA')
% GM_printEPS(300,300,'SA_diffSignal_SA')
% 
% %% Diff signal Vol
% SA_plot(SA.diff_indicat.Vol,'diffSignal',1)
% GM_printBMP(300,300,'SA_diffSignal_Vol')
% GM_printEPS(300,300,'SA_diffSignal_Vol')

%% H/S diff indicator
figure('Visible','off')
b = bar([SA.diff_indicat.H_S.main;SA.diff_indicat.H_S.total]');
b(1).FaceColor = [.75 .75 .75];
b(1).BarWidth = 1;
b(2).BarWidth = 1;
b(2).FaceColor = [0 0 0];
hold on
set(gca,'XTickLabel',var_names)
grid on
ylim([0 1])
GM_printBMP(204*1.3,204,'SA_diffIndic_HS')
GM_printEPS(204*1.3,204,'SA_diffIndic_HS')

%% L/S diff indicator
figure('Visible','off')
b = bar([SA.diff_indicat.L_S.main;SA.diff_indicat.L_S.total]');
b(1).FaceColor = [.75 .75 .75];
b(1).BarWidth = 1;
b(2).BarWidth = 1;
b(2).FaceColor = [0 0 0];
hold on
set(gca,'XTickLabel',var_names)
grid on
ylim([0 1])
GM_printBMP(204*1.3,204,'SA_diffIndic_LS')
GM_printEPS(204*1.3,204,'SA_diffIndic_LS')

%% SA diff indicator
figure('Visible','off')
b = bar([SA.diff_indicat.SA.main;SA.diff_indicat.SA.total]');
b(1).FaceColor = [.75 .75 .75];
b(1).BarWidth = 1;
b(2).BarWidth = 1;
b(2).FaceColor = [0 0 0];
hold on
set(gca,'XTickLabel',var_names)
grid on
ylim([0 1])
GM_printBMP(204*1.3,204,'SA_diffIndic_SA')
GM_printEPS(204*1.3,204,'SA_diffIndic_SA')

%% Vol diff indicator
figure('Visible','off')
b = bar([SA.diff_indicat.Vol.main;SA.diff_indicat.Vol.total]');
b(1).FaceColor = [.75 .75 .75];
b(1).BarWidth = 1;
b(2).BarWidth = 1;
b(2).FaceColor = [0 0 0];
hold on
set(gca,'XTickLabel',var_names)
grid on
ylim([0 1])
GM_printBMP(204*1.3,204,'SA_diffIndic_Vol')
GM_printEPS(204*1.3,204,'SA_diffIndic_Vol')


%%
close all

end