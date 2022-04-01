function [] = timeCal_signalStatsPlot(time,signal_data)
% Plot the statistical variation of a given signal

hold on

% function to fill between 2 areas
% fill_between_lines = @(X,Y1,Y2,C) fill( [X fliplr(X)],  [Y1 fliplr(Y2)], C );
% fill_between_lines(time',min(signal_data,[],2),max(signal_data,[],2),'b');

% plot all simulations
% plot(time,signal_data,'Color',[.7 .7 .7])

% plot specific moments
plot(time,mean(signal_data,2),'k--','LineWidth',1.5)              % mean
plot(time,smooth(median(signal_data,2)),'k-','LineWidth',1.5)     % median
plot(time,mean(signal_data,2) + 1.*std(signal_data,[],2),'k:','LineWidth',1)    % confidence interval
plot(time,mean(signal_data,2) - 1.*std(signal_data,[],2),'k:','LineWidth',1)    % confidence interval
ylim([0 inf])

% plot only min and max simulations
% plot(time',min(signal_data,[],2),'Color',[.7 .7 .7])
% plot(time',max(signal_data,[],2),'Color',[.7 .7 .7])


end