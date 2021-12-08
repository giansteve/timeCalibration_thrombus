function [] = GM_printEPS(width,height,figName)
% print a figure in .eps with exact specifications and save it in new
% folder whose name is "Plots"

% try
%     cd Plots
% catch
%     mkdir Plots
%     cd Plots
% end

name_fig = sprintf('%s.eps',figName);

% set(gca,'FontName','Times New Roman')
set(gcf,'Units', 'points');
set(gcf,'Position', [100 100 width height]);
set(gcf,'PaperUnits','points');
set(gcf,'PaperSize',[width height]);
set(gcf,'PaperPosition',[0 0 width height]);
set(gcf, 'PaperPositionMode', 'auto');
print(name_fig,'-depsc');

% cd ..\

end