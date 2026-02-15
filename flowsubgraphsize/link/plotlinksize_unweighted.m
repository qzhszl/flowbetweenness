clear,clc
% plot the relation between the flow subgraph node size with the real ave degree

% % for small N:-------------------------------------------------------------
fig = figure; 
fig.Position = [100 100 600 300]; 
hold on;

h_bar  = plot(nan, nan, 'k-', 'LineWidth', 2);                        % 黑线
h_star = plot(nan, nan, 'kp', 'MarkerSize', 6, 'LineStyle', 'none','MarkerFaceColor', "k");  % 黑星


% colors = ["#D08082", "#C89FBF", "#62ABC7", "#7A7DB1", "#6FB494", "#D9B382"];
% colors = ["#D08082", "#C89FBF", "#62ABC7", "#7A7DB1", ...
%           "#6FB494", "#D9B382", "#B3C47A", "#A989C2", ...
%           "#E2C572", "#A0A0A0"];

colors = ["#D08082","#6FB494","#D9B382","#7A7DB1","#62ABC7","#A0A0A0"];
% colors = ["#C89FBF","#62ABC7","#B3C47A","#E2C572",];

% count = 1;
% for N = [10,20,30,50,80,100]
%     avg = 0:0.2:8;
%     avg = [avg,N-1];
%     p_vals = avg/(N-1);
% 
%     s_vals = zeros(size(p_vals));
%     for i = 1:length(p_vals)
%         s_vals(i) = compute_S_link_from_ER(N, p_vals(i));
%     end
%     plot(avg, s_vals, 'LineWidth', 2, Color=colors(count))
%     count = count+1;
%     hold on
% end


count = 1;
for N = [10,20,30,50,80,100]
    avg = 0:0.2:N-1;
    p_vals = avg/(N-1);

    s_vals = 1 - (p_vals.^2 + (1-p_vals).^2).^(N-2);
    plot(avg, s_vals, 'LineWidth', 2, Color=colors(count))
    count = count+1;
    hold on
end



count = 1;
for N = [10,20,30,50,80,100]
    filefolder_name = "D:\\data\\flow betweenness\\sizeofflowsubgraph\\new\\unweighted";
    outname = fullfile(filefolder_name, sprintf('%dnode_results_summary.csv', N));
    result_table = readtable(outname);
    x = result_table.RealAveDegree;
    y = result_table.LinkSizeFSG./result_table.LinkNum;
    if N>=60
        idx = 1:2:length(x);
    else
        idx = 1:3:length(x);
    end
    plot(x(idx),y(idx),'Marker', 'p','LineStyle','none', 'MarkerSize', 6,'MarkerEdgeColor',colors(count), 'MarkerFaceColor', colors(count))
    hold on
    count = count+1;
end
ylim([0,1.05]);
xlim([0.075,160]);
xticks([1,10,100])

% lgd = legend({sprintf('simultion:$\\Lambda_l$, $N=%d$',n), s}, 'interpreter','latex','Location', 'best',FontSize=30);
xlabel('$E[D]$',Interpreter='latex',FontSize=16);
ylabel('$E[\rho_L]$','interpreter','latex',FontSize=16)
% set(legend, 'Position', [0.446, 0.73, 0.2, 0.1]);
box on

lgd = legend({'Analytical value', 'Simulation result', ...
        '$N=10$', '$N=20$', '$N=30$', '$N=50$', '$N=80$', '$N=100$'}, ...
        'Interpreter', 'latex', ...
        'Position', [0.225 0.555 0.1 0.1],...
        'FontSize', 14, ...
        'Box', 'on');
lgd.ItemTokenSize = [10, 8];
% set(findobj(lgd, 'type', 'line'), 'LineWidth', 4);
ax = gca;  % Get current axis
ax.FontSize = 14;  % Set font size for tick label
set(gca,"xscale","log")
picname = sprintf("D:\\data\\flow betweenness\\sizeofflowsubgraph\\new\\linksize_flow_subgraph_unweighted.pdf");
exportgraphics(fig, picname,'BackgroundColor', 'none','Resolution', 600);
