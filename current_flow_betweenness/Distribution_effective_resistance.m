clear,clc
% clear,clc

% load BA netwok
G = BAgraph(1000, 4, 3);
data = readmatrix('D:\\data\\flow betweenness\\BAnetworks\\BAnetworkN5000m3.txt');  % 每行: u v w

% 提取三列
s = data(:,1)+1;    % 起点
t = data(:,2)+1;    % 终点
w = data(:,3);    % 权重

% 构造无向图（如果是有向图，就用 digraph）
G = graph(s, t, w);
A = adjacency(G,"weighted");
n = size(A,1);

% [row,col,v]=find(G.Edges.Weight<=0)


% deg = degree(G)
% % h = histogram(deg,60,Normalization="pdf");
% [counts, edges] = histcounts(deg, 'BinMethod','integers');
% k = edges(1:end-1);
% P = counts / sum(counts);
% nz = P>0;
% k = k(nz);
% P = P(nz);
% 
% plot(k,P,"o")
% set(gca,"YScale", "log")
% set(gca,"XScale", "log")

% % Generate an ER network
% avg = 6;
% 
% target_mean = 2/avg
% target_std = sqrt(2/(avg^3))
% 
% n = 1000;
% p = avg/(1000-1)
% A = GenerateERfast(n,p,0);

% 度矩阵与拉普拉斯
deg = sum(A,2);
D = diag(deg);
L = D - A;                % combinatorial Laplacian

% 特征分解（L 对称，使用 eig 更稳妥）
[V, Lambda] = eig(full(L));   % Lambda 对角矩阵
eigvals = diag(Lambda);

% 按特征值从小到大排序（一般第一个约为0）
[evals_sorted, idx] = sort(eigvals, 'ascend');
V = V(:, idx);

% 阈值去掉零（或接近零）特征值，构造伪逆
tol = 1e-12;                 % 数值容忍度，可酌情调整
pos_idx = find(evals_sorted > tol);   % 非零特征值索引

% 若图不连通，L 会有多个零特征值，pos_idx 对应非零那些
Lplus = zeros(n);  % 初始化
for k = 1:length(pos_idx)
    i = pos_idx(k);
    lambda = evals_sorted(i);
    v = V(:, i);
    Lplus = Lplus + (1/lambda) * (v * v.');
end
% Lplus 现在是 Moore-Penrose 伪逆

% 检查对称性与近似性质（可选）
% disp(norm(Lplus - Lplus', 'fro'));  % 应该接近 0
% disp(norm(L * Lplus * L - L, 'fro')); % Moore-Penrose 条件之一, 应接近 0

Reff = zeros(n);
for i = 1:n
    for j = i+1:n
        Reff(i,j) = Lplus(i,i) + Lplus(j,j) - 2*Lplus(i,j);
        Reff(j,i) = Reff(i,j); % 对称
    end
end

resistances = Reff(triu(true(n),1));
writematrix(resistances,'D:\\data\\flow betweenness\\BAnetworks\\BAnetworkN5000m3_effective_resistance.txt')


% 1. 拟合正态分布
mu = mean(resistances);
sigma = std(resistances);

% 或者用 fitdist 更正规：
pd = fitdist(resistances, 'Normal');
mu = pd.mu;
sigma = pd.sigma;

% 2. 画直方图
fig = figure; 
fig.Position = [100 100 600 600];
colors = ["#D08082", "#C89FBF", "#62ABC7", "#7A7DB1", "#6FB494", "#D9B382"];
% for count  =1:6
%     plot(linspace(1,10,10),count*ones(10,1),Color=colors(count))
%     hold on
% end


histogram(resistances, 50, 'Normalization', 'pdf', FaceColor='#7A7DB1'); % 归一化成概率密度
hold on;

% 3. 画拟合的正态分布曲线
x = linspace(min(resistances), max(resistances), 200);
y = normpdf(x, mu, sigma);
plot(x, y, "Color",'#D08082', 'LineWidth', 5);

ax = gca;  % Get current axis
ax.FontSize = 30;  % Set font size for tick label
box on


xlim_ = xlim; ylim_ = ylim;
dx = 0.02*(xlim_(2)-xlim_(1));   % 横向内缩 2%
dy = 0.03*(ylim_(2)-ylim_(1));   % 纵向下移 5%
str = {'$N = 10^{3}$', '$E[D] = 6$'};  % cell 数组自动分两行
text(xlim_(2)-dx, ylim_(2)-dy, str, ...
     'Interpreter','latex', ...
     'VerticalAlignment','top', ...
     'HorizontalAlignment','right', ...
     'FontSize',30);

dx = 0.02*(xlim_(2)-xlim_(1));
dy = 0.08*(ylim_(2)-ylim_(1));
% LaTeX 格式的字符串
str = sprintf('Gaussian:\n$\\mu = %.4f$\n $\\sigma = %.4f$', mu, sigma);
% 在右下角加文字
text(xlim_(2)-dx, ylim_(1)+dy, str, ...
     'Interpreter','latex', ...
     'HorizontalAlignment','right', ...
     'VerticalAlignment','bottom', ...
     'FontSize',29);

% ylim([0,105])
xlabel('$x$', interpreter = "latex",FontSize=40);
ylabel('$f_{\Lambda_G}(x)$',interpreter = "latex",FontSize=40);
% title(sprintf('Normal Fit: \\mu = %.3f, \\sigma = %.3f', mu, sigma));
fit_legend_name = sprintf('Gaussian Fit:\n \\mu = %.4f, \\sigma = %.4f', mu, sigma);
% legend('Empirical Distribution', fit_legend_name,FontSize=28);
hold off;
% picname = sprintf("D:\\data\\flow betweenness\\power_dissipation\\graphdissipation_distribution_N%d_p%.2f.pdf",n,p);
picname = sprintf("D:\\data\\flow betweenness\\power_dissipation\\graphdissipation_distribution_N%BA.pdf",n);

% exportgraphics(fig, picname,'BackgroundColor', 'none','Resolution', 600);
print(fig, picname, '-dpdf', '-r600', '-bestfit');

