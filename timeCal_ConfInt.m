function [statistics_struct] = timeCal_ConfInt(signal,Ns)
%MOPRH_CONFINT Compute confidence interval for a given output in time
%series
% INPUT:
%       - signal: output to analyse
%       - Ns: number of simulations performed

signal_test = signal;
%%
% figure
% hold on
% plot(mean(signal,2),'r','LineWidth',1)

[signal_test,tf] = rmoutliers(signal');
signal_test = signal_test';
% plot(mean(signal_test,2),'Color',[0 0 0])

[signal_testMedian,tfMedian] = rmoutliers(signal','median');
signal_testMedian = signal_testMedian';
% plot(mean(signal_testMedian,2))

[signal_testMean,tfMean] = rmoutliers(signal','mean');
signal_testMean = signal_testMean';
% plot(mean(signal_testMean,2))

[signal_testQua,Qua] = rmoutliers(signal','quartiles');
signal_testQua = signal_testQua';
% plot(mean(signal_testQua,2))

[signal_testGru,tfGru] = rmoutliers(signal','grubbs');
signal_testGru = signal_testGru';
% plot(mean(signal_testGru,2))

[signal_testGesd,tfGesd] = rmoutliers(signal','gesd');
signal_testGesd = signal_testGesd';
% plot(mean(signal_testGesd,2))
% 
% lgd = {'Original',...
%     sprintf('Norm, Ns:%.0f',size(signal_test,2)),...
%     sprintf('Median, Ns:%.0f',size(signal_testMedian,2)),...
%     sprintf('Mean, Ns:%.0f',size(signal_testMean,2)),...
%     sprintf('Quartile, Ns:%.0f',size(signal_testQua,2)),...
%     sprintf('Grubbs, Ns:%.0f',size(signal_testGru,2)),...
%     sprintf('GESD, Ns:%.0f',size(signal_testGesd,2))};
% legend(lgd);


%%
statistics_struct.tf = tfGesd;
statistics_struct.signal = signal_testGesd;
statistics_struct.median = median(signal_testGesd,2);
statistics_struct.mean = mean(signal_testGesd,2);
statistics_struct.SE = std(signal_testGesd,0,2)./sqrt(Ns);
statistics_struct.t_score = tinv([0.025 0.975],Ns-1);
statistics_struct.CI = statistics_struct.mean + statistics_struct.t_score.*statistics_struct.SE;

% figure
% plot(signal_testMean,'Color',[.7 .7 .7])
% hold on
% plot(statistics_struct.mean,'k')
% plot(statistics_struct.median,'r')
end

