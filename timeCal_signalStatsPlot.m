function [] = timeCal_signalStatsPlot(time,signal_data)
% Plot the statistical variation of a given signal



figure('Visible','off')
plot(time,signal_data,'Color',[.7 .7 .7])
hold on
plot(time,mean(signal_data,2),'k-','LineWidth',2)
plot(time,mean(signal_data,2) + 1.*std(signal_data,[],2),'k:','LineWidth',2)
plot(time,mean(signal_data,2) - 1.*std(signal_data,[],2),'k:','LineWidth',2)
ylim([0 inf])



end