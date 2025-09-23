clear,clc
% clear,clc
% % G = BAgraph(1000, 4, 3);
% data = readmatrix('D:\\data\\flow betweenness\\BA_network\\BAnetworkN5000m3.txt');  % 每行: u v w
% 
% % 提取三列
% s = data(:,1)+1;    % 起点
% t = data(:,2)+1;    % 终点
% w = data(:,3);    % 权重
% 
% % 构造无向图（如果是有向图，就用 digraph）
% G = graph(s, t, w);
% A = adjacency(G,"weighted")
% 
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

A = GenerateERfast(1000,0.05,0);
n = 1000;
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

% 展平所有元素并画 histogram
vals = Lplus(:);

figure('Color','w','Units','normalized','Position',[0.2 0.2 0.5 0.5]);
histogram(vals, 100, 'Normalization', 'pdf'); % 80 个 bin
xlabel('Elements of L^+');
ylabel('Probability density');
title('Histogram of all elements in Laplacian pseudoinverse (L^+)');

% 在图上标注均值和标准差
mu = mean(vals);
sigma = std(vals);
yl = ylim;
hold on;
plot([mu mu], yl, '--k', 'LineWidth',1.5);
text(mu, yl(2)*0.9, sprintf('\\leftarrow mean=%.3e', mu), 'FontSize',10);
hold off;

% 输出一些数值信息
fprintf('Number of zero (<= tol) eigenvalues: %d\n', sum(eigvals <= tol));
fprintf('Mean of L^+ elements: %.6e\n', mu);
fprintf('Std  of L^+ elements: %.6e\n', sigma);

% 如果想看对角/非对角元素分布，可以另外绘制：
figure('Color','w','Units','normalized','Position',[0.2 0.1 0.5 0.35]);
histogram(triu(Lplus), 50, 'Normalization','pdf');
% xscale("log")
yscale("log")
xlabel('upper triangle elements of L^+');
ylabel('Density');
title('upper triangle elements of L^+');

