clear,clc
% plot the relation between the flow subgraph node size with the real ave degree

% % for small N:-------------------------------------------------------------
% fig = figure; 
% fig.Position = [100 100 450 300]; 
% hold on;
% 
% h_bar  = plot(nan, nan, 'k-', 'LineWidth', 2);                        % 黑线
% h_star = plot(nan, nan, 'kp', 'MarkerSize', 6, 'LineStyle', 'none','MarkerFaceColor', "k");  % 黑星
% 
% 
% % colors = ["#D08082", "#C89FBF", "#62ABC7", "#7A7DB1", "#6FB494", "#D9B382"];
% % colors = ["#D08082", "#C89FBF", "#62ABC7", "#7A7DB1", ...
% %           "#6FB494", "#D9B382", "#B3C47A", "#A989C2", ...
% %           "#E2C572", "#A0A0A0"];
% 
% colors = ["#D08082","#6FB494","#D9B382","#7A7DB1",];
% % colors = ["#C89FBF","#62ABC7","#B3C47A","#E2C572",];
% 
% count = 1;
% for N = [10,20,30,50]
%     avg = 0:0.1:10;
%     p_vals = avg/(N-1);
%     s_vals = zeros(size(p_vals));
%     for i = 1:length(p_vals)
%         s_vals(i) = compute_S_link_from_ER(N, p_vals(i));
%     end
%     plot(avg, s_vals, 'LineWidth', 2, Color=colors(count))
%     count = count+1;
%     hold on
% end
% 
% count = 1;
% for N = [10,20,30,50]
%     filefolder_name = "D:\\data\\flow betweenness\\sizeofflowsubgraph\\new";
%     outname = fullfile(filefolder_name, sprintf('%dnode_results_summary.csv', N));
%     result_table = readtable(outname);
%     plot(result_table.RealAveDegree,result_table.LinkSizeFSG./result_table.LinkNum,'Marker', 'p','LineStyle','none', 'MarkerSize', 6,'MarkerEdgeColor',colors(count), 'MarkerFaceColor', colors(count))
%     hold on
%     count = count+1;
% end
% 
% % lgd = legend({sprintf('simultion:$\\Lambda_l$, $N=%d$',n), s}, 'interpreter','latex','Location', 'best',FontSize=30);
% xlabel('$E[D]$',Interpreter='latex',FontSize=16);
% ylabel('$E\left[\rho_L \right]$','interpreter','latex',FontSize=16)
% ylim = ([0,1.05]);
% xlim = ([0,10]);
% % set(legend, 'Position', [0.446, 0.73, 0.2, 0.1]);
% box on
% 
% lgd = legend({'Analytical value', 'Simulation result', ...
%         '$N=10$', '$N=20$', '$N=30$', '$N=50$'}, ...
%         'Interpreter', 'latex', ...
%         'FontSize', 14, ...
%         'Location', 'southeast', ...
%         'Box', 'on');
% lgd.ItemTokenSize = [12, 10];
% % set(findobj(lgd, 'type', 'line'), 'LineWidth', 4);
% ax = gca;  % Get current axis
% ax.FontSize = 12;  % Set font size for tick label
% picname = sprintf("D:\\data\\flow betweenness\\sizeofflowsubgraph\\new\\linksize_flow_subgraph_small.pdf");
% exportgraphics(fig, picname,'BackgroundColor', 'none','Resolution', 600);



% % for large N:-------------------------------------------------------------
fig = figure; 
fig.Position = [100 100 450 300]; 
hold on;

h_bar  = plot(nan, nan, 'k-', 'LineWidth', 2);                        % 黑线
h_star = plot(nan, nan, 'kp', 'MarkerSize', 6, 'LineStyle', 'none','MarkerFaceColor', "k");  % 黑星


% colors = ["#D08082", "#C89FBF", "#62ABC7", "#7A7DB1", "#6FB494", "#D9B382"];
% colors = ["#D08082", "#C89FBF", "#62ABC7", "#7A7DB1", ...
%           "#6FB494", "#D9B382", "#B3C47A", "#A989C2", ...
%           "#E2C572", "#A0A0A0"];

% colors = ["#D08082","#6FB494","#D9B382","#7A7DB1",];
colors = ["#C89FBF","#62ABC7","#E2C572","#A0A0A0"];

count = 1;
for N = [80,100,1000,10000]
    avg = 0:0.1:10;
    p_vals = avg/(N-1);
    s_vals = zeros(size(p_vals));
    for i = 1:length(p_vals)
        s_vals(i) = compute_S_link_from_ER(N, p_vals(i));
    end
    % s_vals = s_vals.^1.5
    plot(avg, s_vals, 'LineWidth', 2, Color=colors(count))
    count = count+1;
    hold on
end

count = 1;
for N = [80,100,1000,10000]
    filefolder_name = "D:\\data\\flow betweenness\\sizeofflowsubgraph\\new";
    outname = fullfile(filefolder_name, sprintf('%dnode_results_summary.csv', N));
    result_table = readtable(outname);
    plot(result_table.RealAveDegree,result_table.LinkSizeFSG./result_table.LinkNum,'Marker', 'p','LineStyle','none', 'MarkerSize', 6,'MarkerEdgeColor',colors(count), 'MarkerFaceColor', colors(count))
    hold on
    count = count+1;
end

% lgd = legend({sprintf('simultion:$\\Lambda_l$, $N=%d$',n), s}, 'interpreter','latex','Location', 'best',FontSize=30);
xlabel('$E[D]$',Interpreter='latex',FontSize=16);
ylabel('$E\left[\rho_L \right]$','interpreter','latex',FontSize=16)
ylim = ([0,1.05]);
xlim = ([0,10]);
% set(legend, 'Position', [0.446, 0.73, 0.2, 0.1]);
box on

lgd = legend({'Analytical value', 'Simulation result', ...
        '$N=80$', '$N=100$', '$N=1000$', '$N=10000$'}, ...
        'Interpreter', 'latex', ...
        'FontSize', 14, ...
        'Location', 'southeast', ...
        'Box', 'on');
lgd.ItemTokenSize = [12, 10];
% set(findobj(lgd, 'type', 'line'), 'LineWidth', 4);
ax = gca;  % Get current axis
ax.FontSize = 12;  % Set font size for tick label
picname = sprintf("D:\\data\\flow betweenness\\sizeofflowsubgraph\\new\\linksize_flow_subgraph_large.pdf");
exportgraphics(fig, picname,'BackgroundColor', 'none','Resolution', 600);