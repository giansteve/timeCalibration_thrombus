function [] = GM_pdf_matrix(mat)
% generate probability distribution function (pdf) of a matrix. The second
% 1ension is intended as time, therefore the analysis of the PDF is
% dynamic
% the "time" variable has to be in position 1, e.g. mat[t,any]

%% define colors
colors = zeros(size(mat,1),3);
for i = 1:size(mat,1)
    colors(i,:) = [1-i/size(mat,1) 1-i/size(mat,1) 1-i/size(mat,1)];
%     colors(i,:) = [1-i/size(mat,1) 0.5 0.5];
end
colors(colors==0)=.05;
colormap gray
c = colorbar;
% set( c, 'YDir', 'reverse' );

%% define max and min of variable
N_edges = 100;
if min(min(mat)) == -inf
    edges_hist = linspace(-1,max(max(mat)),N_edges);
else
    edges_hist = linspace(min(min(mat)),max(max(mat)),N_edges);
end
centres = edges_hist(1:end-1)+ diff(edges_hist)/2;

% figure
hold on
for i = 1:size(mat,1)
    h = histcounts(mat(i,:),edges_hist,'Normalization','probability');
    plot(centres,smooth(h),'Color',colors(i,:))
end
xlim([min(edges_hist) max(edges_hist)])
grid on










end