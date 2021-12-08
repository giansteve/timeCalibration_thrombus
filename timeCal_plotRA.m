function timeCal_plotRA(RA,SA,var_names)
% function to plot the RA results


% H/S
% Define RVs to analyse
[~,idx_cache] = sort(SA.diff_indicat.H_S.total);
nonConst = idx_cache(end-3:end);
% plot
figure('Visible','off')
subplot(121)
scatter(RA.H_S.Results.History.X(:,nonConst(end-1)),RA.H_S.Results.History.X(:,nonConst(end)),1,RA.H_S.Results.History.G,'filled')
cc = colorbar;
colormap(winter)
ylabel(cc,'indicator H/S','Interpreter','latex')
xlabel(var_names{nonConst(end-1)})
ylabel(var_names{nonConst(end)})
subplot(122)
scatter(RA.H_S.Results.History.X(:,nonConst(end-3)),RA.H_S.Results.History.X(:,nonConst(end-2)),1,RA.H_S.Results.History.G,'filled')
cc = colorbar;
colormap(winter)
ylabel(cc,'indicator H/S','Interpreter','latex')
xlabel(var_names{nonConst(end-3)})
ylabel(var_names{nonConst(end-2)})
GM_printBMP(400,250,'RA_indicHS')
GM_printEPS(400,250,'RA_indicHS')
close all

% L/S
% Define RVs to analyse
[~,idx_cache] = sort(SA.diff_indicat.L_S.total);
nonConst = idx_cache(end-3:end);
% plot
figure('Visible','off')
subplot(121)
scatter(RA.L_S.Results.History.X(:,nonConst(end-1)),RA.L_S.Results.History.X(:,nonConst(end)),1,RA.L_S.Results.History.G,'filled')
cc = colorbar;
colormap(winter)
ylabel(cc,'indicator L/S','Interpreter','latex')
xlabel(var_names{nonConst(end-1)})
ylabel(var_names{nonConst(end)})
axis tight
subplot(122)
scatter(RA.L_S.Results.History.X(:,nonConst(end-3)),RA.L_S.Results.History.X(:,nonConst(end-2)),1,RA.L_S.Results.History.G,'filled')
cc = colorbar;
colormap(winter)
ylabel(cc,'indicator L/S','Interpreter','latex')
xlabel(var_names{nonConst(end-3)})
ylabel(var_names{nonConst(end-2)})
axis tight
GM_printBMP(400,250,'RA_indicLS')
GM_printEPS(400,250,'RA_indicLS')
close all

% SA
% Define RVs to analyse
[~,idx_cache] = sort(SA.diff_indicat.SA.total);
nonConst = idx_cache(end-3:end);
% plot
figure('Visible','off')
subplot(121)
scatter(RA.SA.Results.History.X(:,nonConst(end-1)),RA.SA.Results.History.X(:,nonConst(end)),1,RA.SA.Results.History.G,'filled')
cc = colorbar;
colormap(winter)
ylabel(cc,'indicator SA','Interpreter','latex')
xlabel(var_names{nonConst(end-1)})
ylabel(var_names{nonConst(end)})
subplot(122)
scatter(RA.SA.Results.History.X(:,nonConst(end-3)),RA.SA.Results.History.X(:,nonConst(end-2)),1,RA.SA.Results.History.G,'filled')
cc = colorbar;
colormap(winter)
ylabel(cc,'indicator SA','Interpreter','latex')
xlabel(var_names{nonConst(end-3)})
ylabel(var_names{nonConst(end-2)})
GM_printBMP(400,250,'RA_indicSA')
GM_printEPS(400,250,'RA_indicSA')
close all

% Vol
% Define RVs to analyse
[~,idx_cache] = sort(SA.diff_indicat.Vol.total);
nonConst = idx_cache(end-3:end);
% plot
figure('Visible','off')
subplot(121)
scatter(RA.Vol.Results.History.X(:,nonConst(end-1)),RA.Vol.Results.History.X(:,nonConst(end)),1,RA.Vol.Results.History.G,'filled')
cc = colorbar;
colormap(winter)
ylabel(cc,'indicator Vol','Interpreter','latex')
xlabel(var_names{nonConst(end-1)})
ylabel(var_names{nonConst(end)})
subplot(122)
scatter(RA.Vol.Results.History.X(:,nonConst(end-3)),RA.Vol.Results.History.X(:,nonConst(end-2)),1,RA.Vol.Results.History.G,'filled')
cc = colorbar;
colormap(winter)
ylabel(cc,'indicator Vol','Interpreter','latex')
xlabel(var_names{nonConst(end-3)})
ylabel(var_names{nonConst(end-2)})
GM_printBMP(400,250,'RA_indicVol')
GM_printEPS(400,250,'RA_indicVol')

close all
end % function