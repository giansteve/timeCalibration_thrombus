function [Pi,T_m,T_lb,T_ub] = PAWNmain(Y_uncond,output_to_const,Ns_const,int_num,var_names,Ereg_handle)
% main function for setting up the PAWN sensitivity analysis


% morph_PAWNmain(Y_uncond,output_to_const,Ns_const,int_num)
M = size(output_to_const.ExpDesign.X,2);
alfa = 0.01;
funct_handle_KS = 'max';
Nboot = 100  ; % number of boostrap resamples
Eplot_PAWN = 0;

% to morph_PAWNmain.m
% Create the conditional output
% pceNonKNST.output = morph_PCE_constRV_eval(99,Ns_const,output_to_const);
for m = 1:M
    k_NonConst = m;
    pceKNST(m).output = PCE_PAWN_eval(k_NonConst,Ns_const,output_to_const,int_num);
    for cache = 1:int_num
        XX_condInput{m,cache} = pceKNST(m).output(cache).ED;
        YY_condInput{m,cache} = pceKNST(m).output(cache).eval;
        xc_temp(1,cache) = XX_condInput{m,cache}(1,1);
    end
    xc{m,1} = xc_temp;
end
Xu = output_to_const.ExpDesign.X;
Yu = Y_uncond';
NU = size(Yu,1);
distrpar=cell(M,1); for i=1:M; distrpar{i}=[min(Xu(:,i)) max((Xu(:,i)))]; end
% Estimate unconditional and conditional CDFs:
[ YF, Fu, Fc  ] = pawn_cdfs(Yu,YY_condInput) ;
% Plot CDFs:
if Eplot_PAWN == 1
    figure
else
    figure('Visible','off')
end
if mod(M,2)==0
    for i=1:M
        subplot(2,M/2,i)
        pawn_plot_cdf(YF, Fu, Fc(i,:),[],var_names{i})
    end
else
    for i=1:M
        subplot(1,M,i)
        pawn_plot_cdf(YF, Fu, Fc(i,:),[],var_names{i})
    end
    
end
% Compute KS statistics:
KS = pawn_ks(YF,Fu,Fc) ;
% Plot KS statistics:
if Eplot_PAWN == 1
    figure
else
    figure('Visible','off')
end
if mod(M,2)==0
    for i=1:M
        subplot(2,M/2,i)
        pawn_plot_kstest(KS(:,i),Ns_const,size(Yu,1),alfa,xc{i}',var_names{i})
    end
else
    for i=1:M
        subplot(1,M,i)
        pawn_plot_kstest(KS(:,i),Ns_const,size(Yu,1),alfa,xc{i}',var_names{i})
    end
end
% Compute PAWN index by taking a statistic of KSs (e.g. max):
Pi = max(KS);
% Plot:
if Eplot_PAWN == 1
    figure
else
    figure('Visible','off')
end
boxplot1(Pi,var_names)
% Use bootstrapping to assess robustness of PAWN indices:
stat = funct_handle_KS ; % statistic to be applied to KSs
% Nboot = 100  ; % number of boostrap resamples
[ T_m, T_lb, T_ub ] = pawn_indices(Yu,YY_condInput,stat,[],Nboot);
% Plot:
% if Eplot_PAWN == 1
figure
% else
%     figure('Visible','off')
% end
boxplot1(T_m,var_names,[],T_lb,T_ub)
% Convergence analysis:
stat = funct_handle_KS ; % statistic to be applied to KSs
NCb = [ Ns_const/10 Ns_const/2 Ns_const ] ;
NUb = [ NU/10 NU/2 NU ] ;
[ T_m_n, T_lb_n, T_ub_n ] = pawn_convergence(Yu,YY_condInput,stat,NUb,NCb,[],Nboot );
n = int_num;
NN = NUb+n*NCb ;
if Eplot_PAWN == 1
    figure
else
    figure('Visible','off')
end
plot_convergence(T_m_n,NN,T_lb_n,T_ub_n,[],'no of evals',[],var_names)


%% My version - some problems
% Pi = zeros(size(pceKNST(1).output(1).eval,2),M);
% T_m = Pi;
% T_lb = Pi;
% T_ub = Pi;
% pceCELL_PAWN_CDF = cell(M,int_num);
% for tt = 1:size(pceKNST(1).output(1).eval,2)
%     for m = 1:M
%         for ntr = 1:int_num
%             pceCELL_PAWN(tt).CDF{m,ntr} = pceKNST(m).output(ntr).eval(:,tt);
%         end
%     end
% end
%
% if Ereg_handle.switchMain == 1
%     tt = Ereg_handle.timeInstance;
%
%     Pi = zeros(1,M);
%     T_m = Pi;
%     T_lb = Pi;
%     T_ub = Pi;
%     Ereg_handle.switch = 1;
%     Y_cond = pceCELL_PAWN(tt).CDF;
%     Y_uncond_PAWN = Y_uncond(:,tt);
%     [Pi,T_m,T_lb,T_ub] = PAWNsens...
%         (Y_uncond_PAWN,Y_cond,alfa,funct_handle_KS,Nboot,var_names,Eplot_PAWN,Ereg_handle);
% else
%     if size(Y_uncond,2) == 1 %% not working! check it
%         Ereg_handle.switch = 0;
%         Y_cond = pceCELL_PAWN.CDF;
%         Y_uncond_PAWN = Y_uncond;
%         [Pi,T_m,T_lb,T_ub] = PAWNsens...
%             (Y_uncond_PAWN,Y_cond,alfa,funct_handle_KS,Nboot,var_names,Eplot_PAWN,Ereg_handle);
%
%     else
%         for tt = 1:size(pceKNST(1).output(1).eval,2)
%             %     if tt == Ereg_handle.timeInstance && Ereg_handle.switchMain == 1
%             %         Ereg_handle.switch = 1;
%             %         Y_cond = pceCELL_PAWN(tt).CDF;
%             %         Y_uncond_PAWN = Y_uncond(:,tt);
%             %         [Pi(tt,:),T_m(tt,:),T_lb(tt,:),T_ub(tt,:)] = morph_PAWNsens...
%             %             (Y_uncond_PAWN,Y_cond,alfa,funct_handle_KS,Nboot,var_names,Eplot_PAWN,Ereg_handle);
%             %     else
%             Ereg_handle.switch = 0;
%             Y_cond = pceCELL_PAWN(tt).CDF;
%             Y_uncond_PAWN = Y_uncond(:,tt);
%             [Pi(tt,:),T_m(tt,:),T_lb(tt,:),T_ub(tt,:)] = PAWNsens...
%                 (Y_uncond_PAWN,Y_cond,alfa,funct_handle_KS,Nboot,var_names,Eplot_PAWN,Ereg_handle);
%             %     end
%         end
%     end
% end
% % Y_uncond = outToPCE.VFT_Tot_Dev.VFT_D(end,1:end-2)';



end