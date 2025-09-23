clear,clc
filefolder_name = "D:\\data\\flow betweenness\\";

n = 20;
p = 0.66;
resname  = sprintf('power_dissipation_N%dp%.2fER_unweighted.mat',n,p);
filename = filefolder_name+resname;

% load_and_plot_powerdissipation.m
S = load(filename);
results = S.results;

% % 2. 拼接所有实验的 energy
total_energy_path = [];
for r = 1:numel(results)
    total_energy_path = [total_energy_path; results(r).total_SP(:)];
end


total_energy_flow = [];
for r = 1:numel(results)
    total_energy_flow = [total_energy_flow; results(r).total_Flow(:)];
end


% % % 3. Figure. 1 plot the distribution of total energy 
% figure;
% histogram(total_energy_path, 30, 'Normalization', 'pdf'); % 50 bins
% % set(gca,"YScale", "log")
% hold on
% 
% histogram(total_energy_flow, 30, 'Normalization', 'pdf'); % 50 bins
% 
% ylabel('$f_e(x)$','interpreter','latex','FontSize',30)
% xlabel('$x$','interpreter','latex','FontSize',30);



% % 2. load power dissipation for each link
link_energy_path = [];
for r = 1:numel(results)
    link_energy_path = [link_energy_path; results(r).linkP_SP(:)];
end
link_energy_path = link_energy_path./nchoosek(n,2);

link_energy_flow = [];
for r = 1:numel(results)
    link_energy_flow = [link_energy_flow; results(r).linkP_Flow(:)];
end
link_energy_flow = link_energy_flow./nchoosek(n,2);

% % 3. 画直方图 (分布图)
fig = figure; 
fig.Position = [100 100 900 600]; 
colors = ["#D08082", "#C89FBF", "#62ABC7", "#7A7DB1", "#6FB494", "#D9B382"];

% histogram(link_energy_path, 60, 'Normalization', 'pdf', FaceColor='#D08082'); % 50 bins
% hold on

h = histogram(link_energy_flow, 60, 'Normalization', 'pdf', FaceColor='#7A7DB1'); % 50 bins
hold on

% set(gca,"XScale", "log")
% set(gca,"YScale", "log")

ylabel('$f_{E_l}(x)$','interpreter','latex','FontSize',30)
xlabel('$x$','interpreter','latex','FontSize',30);


ax = gca;  % Get current axis
ax.FontSize = 20;  % Set font size for tick label
% xlim([0.01 0.55])
% ylim([0.05 0.25])
% xticks([1 2 3 4])
% xticklabels({'10','20','50','100'})
% lgd = legend({'$N = 20$', '$N = 50$', '$N = 100$', '$N = 200$'}, 'interpreter','latex','Location', 'northwest',FontSize=30);
% lgd.NumColumns = 2;
% set(legend, 'Position', [0.446, 0.73, 0.2, 0.1]);
box on
% set(fig, 'Color', 'none');              % figure 背景透明
% set(gca,  'Color', 'none');             % 坐标轴区域背景透明
hold off

picname = sprintf("D:\\data\\flow betweenness\\power_dissipation\\link_power_distribution_N%d_p%.2f.pdf",n,p);
% exportgraphics(fig, picname,'BackgroundColor', 'none','Resolution', 600);
print(fig, picname, '-dpdf', '-r600', '-bestfit');

function gamma_curve_fit
x  =h.BinEdges
y = h.Values
y  = y(1:40)
x = x(2:41)

plot(x, y, 'o',"Color","black"); 
xdata = x(:);
ydata = y(:);
plot(x, y, 'o',"Color","black"); 
hold on
% plot(xdata, ydata, 'o'); hold on
% 定义 Gamma PDF 模型: b(1)=shape, b(2)=scale
gammaFun = @(b,x) gampdf(x, b(1), b(2));

% 初始猜测
b0 = [1, 1];

% 约束边界，确保参数合法
lb = [1e-6, 1e-6];   % 下界：避免为0或负数
ub = [Inf, Inf];     % 上界：无穷大

% 拟合
opts = optimset('Display','off');  % 不输出迭代信息
beta = lsqcurvefit(gammaFun, b0, xdata, ydata, lb, ub, opts);

shape = beta(1);
scale = beta(2);

% 绘图看结果
xplot = linspace(min(xdata), max(xdata), 300);
yfit = gampdf(xplot, shape, scale);


plot(xplot, yfit, '-r', 'LineWidth', 2);
title(sprintf('shape=%.3f, scale=%.3f', shape, scale));
end


