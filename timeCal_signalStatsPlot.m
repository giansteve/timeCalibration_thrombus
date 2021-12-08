function [] = timeCal_signalStatsPlot(time,signal_data,stats_sign)
% Plot the statistical variation of a given signal



figure('Visible','off')
plot(time,signal_data,'Color',[.7 .7 .7])
hold on
plot(time,stats_sign.mean,'k-','LineWidth',2)
% hold on
plot(time,stats_sign.CI(:,1),'k:','LineWidth',2)
plot(time,stats_sign.CI(:,2),'k:','LineWidth',2)




end