clear,clc
avg = 200;

target_mean = 2/avg
target_std = sqrt(2/(avg^3))

n = 10000;
p = avg/(1000-1)
A = GenerateERfast(n,p,0);

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
% resistances = resistances+resistances.';

% 或者用 fitdist 更正规：
pd = fitdist(resistances, 'Normal');
mu = pd.mu;
sigma = pd.sigma;

% 2. Plot histogram for resistance
figure;
histogram(resistances, 30, 'Normalization', 'pdf'); % 归一化成概率密度
hold on;
% 3. 画拟合的正态分布曲线
x = linspace(min(resistances), max(resistances), 200);
y = normpdf(x, mu, sigma);
plot(x, y, 'r-', 'LineWidth', 2);

xlabel('Effective Resistance');
ylabel('Probability Density');
title(sprintf('Normal Fit: \\mu = %.3f, \\sigma = %.3f', mu, sigma));
legend('Empirical Distribution', 'Normal Fit');
hold off;



% 2. Plot histogram for ωis−ωjs+ωjt−ωit
s = randi(n)
t = randi(n)

col_s = Reff(:, s);   % n x 1
col_t = Reff(:, t);   % n x 1

% 构造 f(i,j) 矩阵
F = col_s - col_t;     % n x 1 向量
G = -col_s + col_t;    % n x 1 向量
% 组合成 n×n 矩阵：F(i) + G(j)
X_result = F + G.';      % 每个元素就是 ωis−ωjs+ωjt−ωit

X_his = X_result(triu(true(n),1));

% numBins = 100;
% [counts, edges] = histcounts(X_his, numBins);
% validIdx = counts > 1000;
% validEdges = edges([find(validIdx), find(validIdx)+1])
X_his = X_his(X_his >= -0.002 & X_his <= 0.003);



% 用 fitdist：
pd = fitdist(X_his, 'Normal');
mu = pd.mu;
sigma = pd.sigma;

figure;

% h = histogram(X_his, 100, 'Normalization', 'pdf'); % 归一化成概率密度
h = histogram(X_his, 30, 'Normalization', 'count');
h.Values
h.BinEdges
hold on;

% 3. 画拟合的正态分布曲线
x = linspace(min(X_his), max(X_his), 200);
y = normpdf(x, mu, sigma);
plot(x, y, 'r-', 'LineWidth', 2);

xlabel('X');
ylabel('Probability Density');
% set(gca,"YScale",'log')
title(sprintf('Normal Fit: \\mu = %.3f, \\sigma = %.3f', mu, sigma));
legend('Empirical Distribution', 'Normal Fit');
hold off;



