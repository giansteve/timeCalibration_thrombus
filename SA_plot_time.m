function [] = SA_plot_time(strc,time_vec,pic_name,x_label_text)
%SA_plot It displays the results of the generalized sensitivity analysis in
%time
% INPUT:
% 		- strc: structure with the sensitivity indeces
%       - pic_name: name of the picture if needed


if nargin > 2
    figure('Name',pic_name,'Visible','off')
else
    figure('Visible','off')
end
subplot(2,1,1)
if size(strc.main,2) > 4 && size(strc.main,2) <= 10
    p = plot(time_vec,strc.main,'LineWidth',1);
    for m = 6:size(strc.main,2)
        p(m).LineStyle = '-.';
    end
else
    p = plot(time_vec,strc.main,'LineWidth',1);
end
ylim([0 1])
if isempty(round(time_vec(max(find(isnan(strc.main(:,1)))) + 1)))
    xlim([0 round(max(time_vec))])
elseif strcmpi(pic_name,'VFT Tot Dev') || strcmpi(pic_name,'VFT Tot Dev I') || strcmpi(pic_name,'VFT Tot Dev II') || strcmpi(pic_name,'VFT Tot Dev III')
    xlim([(time_vec(max(find(isnan(strc.main(:,1)))) + 1)) round(max(time_vec))])
else
    xlim([round(time_vec(max(find(isnan(strc.main(:,1)))) + 1)) round(max(time_vec))])
end
grid on
xline(6);
xlabel(x_label_text)
ylabel('$S_{\mathrm{{i}}}$ [-]','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')

subplot(2,1,2)
if size(strc.main,2) > 4 && size(strc.main,2) <= 10
    p = plot(time_vec,strc.total,'LineWidth',1);
    for m = 6:size(strc.main,2)
        p(m).LineStyle = '-.';
    end
else
    p = plot(time_vec,strc.main,'LineWidth',1);
end
ylim([0 1])
if isempty(round(time_vec(max(find(isnan(strc.main(:,1)))) + 1)))
    xlim([0 round(max(time_vec))])
elseif strcmpi(pic_name,'VFT Tot Dev') || strcmpi(pic_name,'VFT Tot Dev I') || strcmpi(pic_name,'VFT Tot Dev II') || strcmpi(pic_name,'VFT Tot Dev III')
    xlim([(time_vec(max(find(isnan(strc.main(:,1)))) + 1)) round(max(time_vec))])
else
    xlim([round(time_vec(max(find(isnan(strc.main(:,1)))) + 1)) round(max(time_vec))])
end
grid on
xlabel(x_label_text)
xline(6);
ylabel('$S^{\mathrm{{T}}}_{\mathrm{{i}}}$ [-]','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')


end