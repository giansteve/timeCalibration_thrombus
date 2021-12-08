function [] = timeCal_plotPCEvsMRI(PCE,outToPCE,human_thr,figName)
% plot the MRI data with error bars and the computed PCE with errorbars
time_array = [1 11 21 31 41 51 61]; % indices for model instances


figure('Visible','off')
%% H/S
subplot(221)
% MRI
plot(linspace(0,1,size(human_thr,1)),human_thr(:,2),'r.','MarkerSize',10)
hold on
% model output
plot(linspace(0,1,size(human_thr,1)),mean(outToPCE.H_S_model(time_array,:),2),'b.','MarkerSize',10)
% PCE output
plot(linspace(0,1,size(human_thr,1)),mean(PCE.H_S),'k.','MarkerSize',10)
% errorbars
errorbar(linspace(0,1,size(human_thr,1)),human_thr(:,2),human_thr(:,3)./sqrt(6),human_thr(:,3)./sqrt(6),'LineStyle','none','Marker','none','Color','r')
errorbar(linspace(0,1,size(human_thr,1)),mean(outToPCE.H_S_model(time_array,:),2),...
    std(outToPCE.H_S_model(time_array,:),[],2)./sqrt(size(outToPCE.H_S_model,2)),std(outToPCE.H_S_model(time_array,:),[],2)./sqrt(size(outToPCE.H_S_model,2)),'LineStyle','none','Marker','none','Color','b')
errorbar(linspace(0,1,size(human_thr,1)),mean(PCE.H_S),std(PCE.H_S)./sqrt(size(PCE.H_S,1)),std(PCE.H_S)./sqrt(size(PCE.H_S,1)),'LineStyle','none','Marker','none','Color','k')
% quantiles
plot(linspace(0,1,size(human_thr,1)),quantile(outToPCE.H_S_model(time_array,:),0.05,2),'--b')
plot(linspace(0,1,size(human_thr,1)),quantile(outToPCE.H_S_model(time_array,:),0.95,2),'--b')
plot(linspace(0,1,size(human_thr,1)),quantile(PCE.H_S(time_array,:),0.05,1),'--k')
plot(linspace(0,1,size(human_thr,1)),quantile(PCE.H_S(time_array,:),0.95,1),'--k')
% ylim([0.85 1])
legend('MRI','Model','PCE','Location','best')
ylabel('H/S')
xlabel('t/T')
grid on

%% L/S
subplot(222)
% MRI
plot(linspace(0,1,size(human_thr,1)),human_thr(:,4),'r.','MarkerSize',10)
hold on
% model output
plot(linspace(0,1,size(human_thr,1)),mean(outToPCE.L_S_model(time_array,:),2),'b.','MarkerSize',10)
% PCE output
plot(linspace(0,1,size(human_thr,1)),mean(PCE.L_S),'k.','MarkerSize',10)
errorbar(linspace(0,1,size(human_thr,1)),human_thr(:,4),human_thr(:,4)./sqrt(6),human_thr(:,4)./sqrt(6),'LineStyle','none','Marker','none','Color','r')
errorbar(linspace(0,1,size(human_thr,1)),mean(outToPCE.L_S_model(time_array,:),2),...
    std(outToPCE.L_S_model(time_array,:),[],2)./sqrt(size(outToPCE.L_S_model,2)),std(outToPCE.L_S_model(time_array,:),[],2)./sqrt(size(outToPCE.L_S_model,2)),'LineStyle','none','Marker','none','Color','b')
errorbar(linspace(0,1,size(human_thr,1)),mean(PCE.L_S),std(PCE.L_S)./sqrt(size(PCE.L_S,1)),std(PCE.L_S)./sqrt(size(PCE.L_S,1)),'LineStyle','none','Marker','none','Color','k')
% quantiles
plot(linspace(0,1,size(human_thr,1)),quantile(outToPCE.L_S_model(time_array,:),0.05,2),'--b')
plot(linspace(0,1,size(human_thr,1)),quantile(outToPCE.L_S_model(time_array,:),0.95,2),'--b')
plot(linspace(0,1,size(human_thr,1)),quantile(PCE.L_S(time_array,:),0.05,1),'--k')
plot(linspace(0,1,size(human_thr,1)),quantile(PCE.L_S(time_array,:),0.95,1),'--k')
legend('MRI','Model','PCE','Location','best')
ylabel('L/S')
xlabel('t/T')
grid on

%% Likelihood
subplot(2,2,3)
hold on
plot(linspace(0,1,size(PCE.LOGlikelihood,1)),PCE.LOGlikelihood,'bo','MarkerSize',2)
plot(linspace(0,1,size(outToPCE.diff_signal.L_S.likelihood,2)),log(outToPCE.diff_signal.L_S.likelihood),'ko','MarkerSize',2)
ylabel('$log(1/\sigma^2)$')
xlabel('IN')
grid on
legend('PCE','Model','Location','best')

subplot(2,2,4)
hold on
[h,nbins] = histcounts(PCE.LOGlikelihood,'Normalization','probability');
plot(nbins(1:end-1),(h),'b')
[h,nbins] = histcounts(log(outToPCE.diff_signal.L_S.likelihood),nbins,'Normalization','probability');
plot(nbins(1:end-1),smooth(h),'k')
legend('PCE','Model','Location','best')
ylabel('$f(log(1/\sigma^2))$')
xlabel('$1/\sigma^2$')
grid on


% %% SA
% subplot(223)
% % MRI
% plot(linspace(0,1,size(human_thr,1)),human_thr(:,6),'r.','MarkerSize',10)
% hold on
% % model output
% plot(linspace(0,1,size(human_thr,1)),mean(outToPCE.SA_model(time_array,:),2),'b.','MarkerSize',10)
% % PCE output
% plot(linspace(0,1,size(human_thr,1)),mean(PCE.SA),'k.','MarkerSize',10)
% errorbar(linspace(0,1,size(human_thr,1)),human_thr(:,6),human_thr(:,7)./sqrt(6),human_thr(:,7)./sqrt(6),'LineStyle','none','Marker','none','Color','r')
% errorbar(linspace(0,1,size(human_thr,1)),mean(outToPCE.SA_model(time_array,:),2),...
%     std(outToPCE.SA_model(time_array,:),[],2)./sqrt(size(outToPCE.SA_model,2)),std(outToPCE.SA_model(time_array,:),[],2)./sqrt(size(outToPCE.SA_model,2)),'LineStyle','none','Marker','none','Color','b')
% errorbar(linspace(0,1,size(human_thr,1)),mean(PCE.SA),std(PCE.SA)./sqrt(size(PCE.SA,1)),std(PCE.SA)./sqrt(size(PCE.SA,1)),'LineStyle','none','Marker','none','Color','k')
% % quantiles
% plot(linspace(0,1,size(human_thr,1)),quantile(outToPCE.SA_model(time_array,:),0.05,2),'--b')
% plot(linspace(0,1,size(human_thr,1)),quantile(outToPCE.SA_model(time_array,:),0.95,2),'--b')
% plot(linspace(0,1,size(human_thr,1)),quantile(PCE.SA(time_array,:),0.05,1),'--k')
% plot(linspace(0,1,size(human_thr,1)),quantile(PCE.SA(time_array,:),0.95,1),'--k')
% legend('MRI','Model','PCE','Location','best')
% ylabel('SA')
% xlabel('t/T')
% grid on
% 
% %% VOL
% subplot(224)
% % MRI
% plot(linspace(0,1,size(human_thr,1)),human_thr(:,8),'r.','MarkerSize',10)
% hold on
% % model output
% plot(linspace(0,1,size(human_thr,1)),mean(outToPCE.VOL_model(time_array,:),2),'b.','MarkerSize',10)
% % PCE output
% plot(linspace(0,1,size(human_thr,1)),mean(PCE.VOL),'k.','MarkerSize',10)
% errorbar(linspace(0,1,size(human_thr,1)),human_thr(:,8),human_thr(:,9)./sqrt(6),human_thr(:,9)./sqrt(6),'LineStyle','none','Marker','none','Color','r')
% errorbar(linspace(0,1,size(human_thr,1)),mean(outToPCE.VOL_model(time_array,:),2),...
%     std(outToPCE.VOL_model(time_array,:),[],2)./sqrt(size(outToPCE.VOL_model,2)),std(outToPCE.VOL_model(time_array,:),[],2)./sqrt(size(outToPCE.VOL_model,2)),'LineStyle','none','Marker','none','Color','b')
% errorbar(linspace(0,1,size(human_thr,1)),mean(PCE.VOL),std(PCE.VOL)./sqrt(size(PCE.VOL,1)),std(PCE.VOL)./sqrt(size(PCE.VOL,1)),'LineStyle','none','Marker','none','Color','k')
% % quantiles
% plot(linspace(0,1,size(human_thr,1)),quantile(outToPCE.VOL_model(time_array,:),0.05,2),'--b')
% plot(linspace(0,1,size(human_thr,1)),quantile(outToPCE.VOL_model(time_array,:),0.95,2),'--b')
% plot(linspace(0,1,size(human_thr,1)),quantile(PCE.VOL(time_array,:),0.05,1),'--k')
% plot(linspace(0,1,size(human_thr,1)),quantile(PCE.VOL(time_array,:),0.95,1),'--k')
% legend('MRI','Model','PCE','Location','best')
% ylabel('Vol')
% xlabel('t/T')
% grid on

GM_printBMP(500,500,figName)
GM_printEPS(500,500,figName)



end