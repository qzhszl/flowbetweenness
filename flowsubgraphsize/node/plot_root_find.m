clear,clc

% plot where is the root

p = linspace(0, 1, 100);     % p范围
u_values = [0.5, 1, 2, 5];   % 不同u值

fig = figure; 
fig.Position = [100 100 400 300]; 
hold on;


% colors = ["#D08082", "#C89FBF", "#62ABC7", "#7A7DB1", "#6FB494", "#D9B382"];
% colors = ["#D08082", "#C89FBF", "#62ABC7", "#7A7DB1", ...
%           "#6FB494", "#D9B382", "#B3C47A", "#A989C2", ...
%           "#E2C572", "#A0A0A0"];

h1 = plot(p, p, 'k-', 'LineWidth', 2);

colors = ["#D08082","#6FB494","#D9B382","#7A7DB1","#62ABC7","#A0A0A0"];
% colors = ["#C89FBF","#62ABC7","#B3C47A","#E2C572",];
h_curves = gobjects(1, length(u_values)); % 存储句柄
for k = 1:length(u_values)
    u = u_values(k);
    y = 1 - exp(-u * p);
    h_curves(k) = plot(p, y, 'LineWidth', 2, 'Color', colors(k));
end

for p_star = [0,0.797,0.993]
    y_star = p_star;         % 因为 y = x 在交点处相等
    
    % ---- 标记交点 ----
   
    % ---- 在交点上方添加文字（黑色）----
    if p_star ==0
        plot(p_star, y_star, 'ko', 'MarkerFaceColor', colors(1),'MarkerSize', 6);
        text(p_star+0.2, y_star + 0.05, sprintf('$p^*=0$'), ...
         'Interpreter', 'latex', 'FontSize', 16, ...
         'HorizontalAlignment', 'center', 'Color', 'k');
    elseif p_star==0.797
        plot(p_star, y_star, 'ko', 'MarkerFaceColor', colors(3),'MarkerSize', 6);
        text(p_star-0.15, y_star + 0.05, sprintf('$p^*=%.2f$', p_star), ...
         'Interpreter', 'latex', 'FontSize', 16, ...
         'HorizontalAlignment', 'center', 'Color', 'k');
    else
        plot(p_star, y_star, 'ko', 'MarkerFaceColor', colors(4),'MarkerSize', 6);
        text(p_star-0.15, y_star + 0.08, sprintf('$p^*=%.2f$', p_star), ...
         'Interpreter', 'latex', 'FontSize', 16, ...
         'HorizontalAlignment', 'center', 'Color', 'k');
    end
    
end


leg_text = ["$y=p$", "$E[D]=0.5$", "$E[D]=1$", "$E[D]=2$", "$E[D]=5$"];
lgd =legend([h1 h_curves], leg_text, ...
       'Interpreter','latex', 'FontSize',14, ...
       'Box','on');
lgd.ItemTokenSize = [12, 10];
set(lgd, 'Units', 'normalized');   % 使用相对坐标
ylim([0 1.2])
pos = lgd.Position;
pos(1) = 0.14;                    % 距左边 0.02（越小越贴边）
pos(2) = 1 - pos(4) - 0.08;       % 紧贴顶部
lgd.Position = pos;

% lgd = legend({sprintf('simultion:$\\Lambda_l$, $N=%d$',n), s}, 'interpreter','latex','Location', 'best',FontSize=30);
xlabel('$p$',Interpreter='latex',FontSize=16);
ylabel('$y$','interpreter','latex',FontSize=16)
% set(legend, 'Position', [0.446, 0.73, 0.2, 0.1]);
box on

% lgd = legend({'Analytical value', 'Simulation result', ...
%         '$N=10$', '$N=20$', '$N=30$', '$N=50$', '$N=80$', '$N=100$'}, ...
%         'Interpreter', 'latex', ...
%         'FontSize', 14, ...
%         'Location', 'east', ...
%         'Box', 'on');
% lgd.ItemTokenSize = [12, 10];
% set(findobj(lgd, 'type', 'line'), 'LineWidth', 4);
ax = gca;  % Get current axis
ax.FontSize = 12;  % Set font size for tick label
picname = sprintf("D:\\data\\flow betweenness\\sizeofflowsubgraph\\findroot.pdf");
exportgraphics(fig, picname,'BackgroundColor', 'none','Resolution', 600);
