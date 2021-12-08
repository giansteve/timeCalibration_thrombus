function [Pi,T_m,T_lb,T_ub] = PAWNsens(Y_uncond,Y_cond,alfa,funct_handle_KS,Nboot,var_names,Eplot,Ereg_handle)
% Compute the PAWN sensirtivity indices for the morphological aorta
% project. The code is inspired by the SAFE Toolbox provided by F.Pianosi

int_num = size(Y_cond,2);
NC = size(Y_cond{1},1);
NU = size(Y_uncond,1);
M = size(var_names,2);

% compute CDFs
[ YF, Fu, Fc  ] = pawn_cdfs(Y_uncond,Y_cond);

% Compute KS statistics:
KS = pawn_ks(YF,Fu,Fc) ;
if     alfa == 0.10 ; c = 1.22;
elseif alfa == 0.05 ; c = 1.36;
elseif alfa == 0.025; c = 1.48;
elseif alfa == 0.01 ; c = 1.63;
elseif alfa == 0.005; c = 1.73;
elseif alfa == 0.001; c = 1.95;
else
    error('''alfa'' must take values in {0.1,0.05,0.025,0.01,0.005,0.001}');
end

% Compute PAWN index by taking a statistic of KSs (e.g. max):
Pi = feval(funct_handle_KS,KS);

% Use bootstrapping to assess robustness of PAWN indices:
if Ereg_handle.switch == 1
    if strcmpi(Ereg_handle.CompOp,'between')
        [ T_m, T_lb, T_ub, ks ] = pawn_indices(Y_uncond,Y_cond,funct_handle_KS,[],Nboot,alfa,Ereg_handle.condition,Ereg_handle.Threshold);
    else
        % for regional usage
        [ T_m, T_lb, T_ub, ks ] = pawn_indices(Y_uncond,Y_cond,funct_handle_KS,[],Nboot,alfa,Ereg_handle.CompOp,Ereg_handle.Threshold);
    end
else
    [ T_m, T_lb, T_ub ] = pawn_indices(Y_uncond,Y_cond,funct_handle_KS,[],Nboot);
end
% Convergence analysis:
% NCb = [ NC/10 NC/2 NC ] ;
% NUb = [ NU/10 NU/2 NU ] ;
% NCb = linspace(NC/10,NC,10);
% NUb = linspace(NU/10,NU,10);
% [ T_m_n, T_lb_n, T_ub_n ] = pawn_convergence(outToPCE.VFT_Tot_Dev.VFT_D(end,:)',pceCELL_PAWN_CDF, funct_handle_KS, NUb, NCb,[],Nboot );
% NN = NUb+int_num*NCb ;
% figure; plot_convergence(T_m_n,NN,T_lb_n,T_ub_n,[],'no of evals',[],var_names)
% figure()
% plot()


if Eplot == 1
    ks_crit = zeros(M,1);
    % Plot CDFs:
    figure()
    for m = 1:M
        if mod(M,2) == 0
            subplot(2,M/2,m)
        else
            subplot(1,M,m)
        end
        ecdf(Y_uncond)
        %         [f_orig,x_orig] = ecdf(pceNonKNST.VFT_Tot_Dev.eval(:,20));
        hold on
        for ntr = 1:int_num
            ecdf(Y_cond{m,ntr})
            %             [f_pawn(ntr,m).cdf,x_pawn(ntr,m).cdf] = ecdf(pceKNST(m).VFT_Tot_Dev(ntr).eval(:,20));
        end
        title(var_names{m})
        xlabel('Y')
        ylabel('F(Y)')
    end
    
    % Plot KS stat for each RV
    figure()
    for m = 1:M
        ks_crit(m) = c*sqrt((NC+NU)/(NC*NU)) ;
        if mod(M,2) == 0
            subplot(2,M/2,m)
        else
            subplot(1,M,m)
        end
        plot(KS(:,m),'k','Marker','x','LineWidth',1.5)
        hold on
        %     yline(median(KS(:,m)),'r--');
        yline(feval(funct_handle_KS,KS(:,m)),'b-');
        yline(ks_crit(m),'r--','LineWidth',1.5);
        ylim([0 1])
        xlim([1 int_num])
        title(var_names{m})
        ylabel(sprintf('%s(KS)',funct_handle_KS))
    end
    
    % Plot sensitivity indices with booting limits
    figure()
    bar(Pi,'BarWidth',0.5,'FaceColor',[.7 .7 .7])
    %     hold on
    %     bar(Pi,'BarWidth',0.5,'FaceColor',[1 1 1])
    set(gca,'XTickLabel',var_names)
    % ylim([0 round(max(Pi).*1.2,2)])
    ylim([0 1])
    hold on
    ylabel('$\hat{\mathrm{T}}_\mathrm{i}$')
    er = errorbar(1:M,Pi,T_m - T_lb,T_m - T_ub);
    er.Color = [0 0 0];
    er.LineStyle = 'none';
end

%% TEST on implementing my own PAWN SA
% % the code that is implemented seems to have some small errors in the
% % computation of the indices, which lead to the use of the SAFE toolbox
% % get the VFT tot dev CDF
% [VFTTotDev_cdf,x_VFTTotDev_cdf] = ecdf(outToPCE.VFT_Tot_Dev.VFT_D(20,:));
% output_to_const = PCE.VFT_Tot_Dev;
% Ns_const = 1000;
% int_num = 10;
% pceNonKNST.VFT_Tot_Dev = morph_PCE_constRV_eval(99,Ns_const,output_to_const);
% for m = 1:M
%     k_NonConst = m;
%     pceKNST(m).VFT_Tot_Dev = morph_PCE_PAWN_eval(k_NonConst,Ns_const,output_to_const,int_num);
% end
% % plot test on stats output_PCE
% xlimit = linspace(0,100,41);
% time_edges = 0:.05:1.0; % define normal. time interval's edges
% normTime_mat = outToPCE.VFT_Tot_Dev.normTime_mat;
% normTime_VFT = outToPCE.VFT_Tot_Dev.normTime_VFT;
%
% figure()
% subplot(121)
% idx_col = 1;
% idx_intensity = 1:length(time_edges);
% interval_pdf = linspace(0,100,41);
% for lmt_idx = 2:length(time_edges)
%     % get normalised time interval
%     for dev_th = 1:size(normTime_mat,2)
%         limit1 = time_edges(lmt_idx - 1);
%         limit2 = time_edges(lmt_idx);
%         time_intervalIdx = find(normTime_mat(:,dev_th) >= limit1 & normTime_mat(:,dev_th) <= limit2);
%         discr_VFTD(lmt_idx - 1,dev_th) = mean(normTime_VFT(time_intervalIdx,dev_th))*100;
%         if isnan(discr_VFTD(lmt_idx - 1,dev_th))
%             discr_VFTD(lmt_idx - 1,dev_th) = discr_VFTD(lmt_idx - 1 - 1,dev_th);
%         end
%     end
%     prob_VFTD(lmt_idx - 1,:) = histcounts(discr_VFTD(lmt_idx - 1,:),'BinEdges',interval_pdf,'Normalization','probability');
%     pp(lmt_idx - 1) = plot(interval_pdf(2:end)-interval_pdf(2)/2,smooth(smooth(prob_VFTD(lmt_idx - 1,:))),'Color',[1-idx_intensity(idx_col)/length(time_edges) 1-idx_intensity(idx_col)/length(time_edges) 1-idx_intensity(idx_col)/length(time_edges)],'LineWidth',1);
%     if lmt_idx == length(time_edges)
%         %         pp(lmt_idx - 1).Marker = '.';
%         pp(lmt_idx - 1).Color = 'k';
%         pp(lmt_idx - 1).LineWidth = 1.5;
%         pp(lmt_idx - 1).LineStyle = '--';
%     elseif round(time_edges(lmt_idx - 1),3) == 1
%         pp(lmt_idx - 1).Marker = '.';
%         pp(lmt_idx - 1).Color = 'k';
%         pp(lmt_idx - 1).LineWidth = 1.5;
%     end
%     hold on
%     idx_col = idx_col + 1;
%     lgd_text{idx_col - 1} = sprintf('$t/T_D$ = %.2f',time_edges(idx_col - 1));
% end
% title('Output')
% ylim([0 .2])
% subplot(122)
% for tt = 1:size(pceNonKNST.VFT_Tot_Dev.eval,2)
% % h2 = histcounts(pceNonKNST.VFT_Tot_Dev.eval(:,tt),xlimit,'Normalization','probability');
% h2 = histcounts(pceNonKNST.VFT_Tot_Dev.eval(:,tt),'Normalization','probability');
% plot(((h2)),'Color',[1 1 1] - tt/size(outToPCE.VFT_Tot_Dev.VFT_D,1))
% hold on
% end
% ylim([0 .2])
% title('PCE')
% % PAWN sensitivity analysis
% % plot test CDF
% figure()
% for m = 1:M
%     subplot(2,4,m)
%     ecdf(pceNonKNST.VFT_Tot_Dev.eval(:,20))
%     [f_orig,x_orig] = ecdf(pceNonKNST.VFT_Tot_Dev.eval(:,20));
%     hold on
%     for ntr = 1:int_num
%         ecdf(pceKNST(m).VFT_Tot_Dev(ntr).eval(:,20))
%         [f_pawn(ntr,m).cdf,x_pawn(ntr,m).cdf] = ecdf(pceKNST(m).VFT_Tot_Dev(ntr).eval(:,20));
%     end
%     xlabel(var_names{m})
% end
% % Kolmogorov-Smirnov test
% for m = 1:M
%     for ntr = 1:int_num
%         [h(ntr,m),t(ntr,m),ks2stat(ntr,m)] = kstest2(pceNonKNST.VFT_Tot_Dev.eval(:,20),...
%             pceKNST(m).VFT_Tot_Dev(ntr).eval(:,20),...
%             'Alpha',0.01);
%     end
%     crKSstat(:,m) = sum(h(:,m))/int_num;
% end
%
% figure()
% for m = 1:M
%     subplot(2,4,m)
%     plot(ks2stat(:,m))
%     hold on
%     yline(median(ks2stat(:,m)),'r--');
%     yline(mean(ks2stat(:,m)),'k-');
%     yline(crKSstat(m),'g','LineWidth',1);
%     ylim([0 1])
%     xlabel(var_names{m})
% end
% % plot sensitivity on PAWN
% figure()
% b = bar(1:M,[mean(ks2stat); crKSstat]');
% b(1).FaceColor = [.75 .75 .75];
% b(1).BarWidth = 1;
% b(2).BarWidth = 1;
% b(2).FaceColor = [0 0 0];
% set(gca,'XTickLabel',var_names)
% grid on
% ylim([0 1])
% legend('$T_\mathrm{i}$','$T_\mathrm{i}^\mathrm{TOT}$','Location','nw')


end