function [] = timeCal_plotFittingMRI(fittedData,human_thr,picName)
% 


figure('Visible','off')
subplot(211)
plot(fittedData.time,fittedData.H_S,'r')
hold on
plot(human_thr.time,human_thr.H_S(:,1),'r.','MarkerSize',10)
errorbar(human_thr.time,human_thr.H_S(:,1),human_thr.H_S(:,2),human_thr.H_S(:,2),'LineStyle','none','Marker','none','Color','r')
ylabel('$H/S$ [-]')
legend('fit','MRI','$\sigma_{MRI}$','Location','se')
grid on
subplot(212)
hold on
plot(fittedData.time,fittedData.L_S,'b')
plot(human_thr.time,human_thr.L_S(:,1),'b.','MarkerSize',10)
errorbar(human_thr.time,human_thr.L_S(:,1),human_thr.L_S(:,2),human_thr.L_S(:,2),'LineStyle','none','Marker','none','Color','b')
xlabel('Time [min]')
ylabel('$L/S$ [-]')
legend('fit','MRI','$\sigma_{MRI}$','Location','se')
grid on

% subplot(2,1,2)
% plot(fittedData(:,1),fittedData(:,4),'g')
% hold on
% plot(fittedData.time,human_thr(:,6),'g.','MarkerSize',10)
% errorbar(fittedData.time,human_thr(:,6),human_thr(:,7),human_thr(:,7),'LineStyle','none','Marker','none','Color','g')
% plot(fittedData(:,1),fittedData(:,5),'m')
% plot(fittedData.time,human_thr(:,8),'m.','MarkerSize',10)
% errorbar(fittedData.time,human_thr(:,8),human_thr(:,9),human_thr(:,9),'LineStyle','none','Marker','none','Color','m')
% set(gca, 'YScale', 'log')
% xlabel('Time [min]')
% ylabel('SA [$m^2$]; Vol [$m^3$]')
% % legend('SA','Vol')
% grid on

GM_printBMP(300,300,picName)
GM_printEPS(300,300,picName)


figure('Visible','off')
subplot(2,1,1)
plot(fittedData.time(2:end),human_thr.growthRate.H_S,'r')
hold on
ylabel('$d/dt(H/S)$ [1/min]')
grid on
subplot(2,1,2)
hold on
plot(fittedData.time(2:end),human_thr.growthRate.L_S,'b')
xlabel('Time [min]')
ylabel('$d/dt(L/S)$ [1/min]')
grid on


% subplot(2,1,2)
% plot(fittedData.time(2:end),human_thr.growthRate.exposedArea,'g')
% hold on
% plot(fittedData.time(2:end),human_thr.growthRate.volume,'m')
% set(gca, 'YScale', 'log')
% xlabel('Time [min]')
% ylabel('$\frac{d}{dt}$(SA),$\frac{d}{dt}$(Vol) [$\frac{m^2}{min}$,$\frac{m^3}{min}$]')
% legend('SA','Vol')
% grid on

GM_printBMP(300,300,sprintf('%s_GR',picName))
GM_printEPS(300,300,sprintf('%s_GR',picName))


close all

end