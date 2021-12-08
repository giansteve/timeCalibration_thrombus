function [] = timeCal_plotFittingMRI(fittedData,growthRate,human_thr,picName)

figure('Visible','off')
subplot(2,1,1)
plot(fittedData(:,1),fittedData(:,2),'r')
hold on
plot(human_thr(:,1),human_thr(:,2),'r.','MarkerSize',10)
errorbar(human_thr(:,1),human_thr(:,2),human_thr(:,3)./sqrt(6),human_thr(:,3)./sqrt(6),'LineStyle','none','Marker','none','Color','r')
plot(fittedData(:,1),fittedData(:,3),'b')
plot(human_thr(:,1),human_thr(:,4),'b.','MarkerSize',10)
errorbar(human_thr(:,1),human_thr(:,4),human_thr(:,5)./sqrt(6),human_thr(:,5)./sqrt(6),'LineStyle','none','Marker','none','Color','b')
xlabel('Time [min]')
ylabel('H/S; L/S [-]')
% legend('H/S','L/S')
grid on

subplot(2,1,2)
plot(fittedData(:,1),fittedData(:,4),'g')
hold on
plot(human_thr(:,1),human_thr(:,6),'g.','MarkerSize',10)
errorbar(human_thr(:,1),human_thr(:,6),human_thr(:,7)./sqrt(6),human_thr(:,7)./sqrt(6),'LineStyle','none','Marker','none','Color','g')
plot(fittedData(:,1),fittedData(:,5),'m')
plot(human_thr(:,1),human_thr(:,8),'m.','MarkerSize',10)
errorbar(human_thr(:,1),human_thr(:,8),human_thr(:,9)./sqrt(6),human_thr(:,9)./sqrt(6),'LineStyle','none','Marker','none','Color','m')
set(gca, 'YScale', 'log')
xlabel('Time [min]')
ylabel('SA [$m^2$]; Vol [$m^3$]')
% legend('SA','Vol')
grid on

GM_printBMP(300,300,picName)
GM_printEPS(300,300,picName)


figure('Visible','off')
subplot(2,1,1)
plot(fittedData(2:end,1),growthRate(:,1),'r')
hold on
plot(fittedData(2:end,1),growthRate(:,2),'b')
xlabel('Time [min]')
ylabel('$d/dt(H/S,L/S)$ [1/min]')
legend('H/S','L/S')
grid on


subplot(2,1,2)
plot(fittedData(2:end,1),growthRate(:,3),'g')
hold on
plot(fittedData(2:end,1),growthRate(:,4),'m')
set(gca, 'YScale', 'log')
xlabel('Time [min]')
ylabel('$\frac{d}{dt}$(SA),$\frac{d}{dt}$(Vol) [$\frac{m^2}{min}$,$\frac{m^3}{min}$]')
legend('SA','Vol')
grid on

GM_printBMP(300,300,sprintf('%s_GR',picName))
GM_printEPS(300,300,sprintf('%s_GR',picName))


close all

end