clear,clc
compare_symmetry_vs_p



function compare_symmetry_vs_p
% 比较 ER 网络中电势相等节点对数量与结构对称节点对数量随 p 的变化
%
% - 电势对称: 求解电流注入方程，统计 v_i = v_j 的节点对
% - 结构对称: 邻接矩阵行完全相同的节点对
% - 结果: 画出两条曲线

    % ---- 参数设置 ----
    N = 100;                     % 节点数
    p_values = linspace(0.05, 1, 15);
    num_trials = 100;            % 随机源汇对数量

    mean_equal_pairs = zeros(size(p_values));    % 电势相等节点对均值
    struct_sym_pairs = zeros(size(p_values));    % 结构对称节点对数

    % ---- 遍历 p ----
    for k = 1:length(p_values)
        p = p_values(k);
        fprintf('Processing p = %.2f\n', p);

        % 生成随机 ER 图
        A = rand(N) < p; 
        A = triu(A,1); 
        A = A + A';          % 无向图
        A = double(A);       % 邻接矩阵
        
         A = [0	0	1	1	1	1
            0	0	1	1	1	1
            1	1	0	1	1	1
            1	1	1	0	1	1
            1	1	1	1	0	1
            1	1	1	1	1	0];
        nodei = 1
        nodej = 5
        visualize_current_flow(A, nodei, nodej,1)
        % ---- 计算结构对称节点对数量 ----
        [pair_count, groups] = count_structural_symmetric_pairs(A);
        struct_sym_pairs(k) = pair_count;
        fprintf('共有 %d 对结构对称节点。\n', pair_count);
        disp(groups)

        % ---- 电势对称节点对 ----
        sym_counts = zeros(num_trials,1);
        for trial = 1:num_trials
            nodes = randperm(N,2);
            s = nodes(1); t = nodes(2);
            sym_counts(trial) = count_equal_potential_pairs(A, s, t);
        end
        mean_equal_pairs(k) = mean(sym_counts);
    end

    % ---- 绘图 ----
    figure; hold on; box on;
    plot(p_values, mean_equal_pairs, '-o', 'LineWidth', 1.8, 'MarkerSize',6);
    plot(p_values, struct_sym_pairs, '-s', 'LineWidth', 1.8, 'MarkerSize',6);
    xlabel('$p$', 'Interpreter','latex', 'FontSize',13);
    ylabel('Number of symmetric node pairs', 'Interpreter','latex', 'FontSize',13);
    legend({'Equal-potential pairs (avg)', 'Structurally symmetric pairs'}, ...
           'Interpreter','latex', 'Location','northwest');
    title(sprintf('N = %d, averaged over %d source-destination pairs', N, num_trials));
    grid on;
end

% ==============================================================
% 辅助函数：电势相等节点对
% ==============================================================
function sym_pair_count = count_equal_potential_pairs(A, s, t)
    n = size(A,1);
    L = diag(sum(A,2)) - A;
    b = zeros(n,1); b(s)=1; b(t)=-1;

    ref = n;
    L_red = L(1:end-1,1:end-1);
    b_red = b(1:end-1);
    v = zeros(n,1);
    v(1:end-1) = L_red \ b_red;

    tol = 1e-12;
    visited = false(n,1);
    sym_pair_count = 0;
    for i = 1:n
        if visited(i), continue; end
        same = abs(v - v(i)) < tol;
        idx = find(same);
        if numel(idx) > 1
            sym_pair_count = sym_pair_count + nchoosek(numel(idx),2);
        end
        visited(idx) = true;
    end
end

% ==============================================================
% 辅助函数：结构对称节点对
% ==============================================================
function [pair_count, sym_groups] = count_structural_symmetric_pairs(A)
%COUNT_STRUCTURAL_SYMMETRIC_PAIRS  计算结构对称节点对数量
%
% 若两个节点的邻居集合（不含自身）完全相同，则它们结构对称。
%
% 输入：
%   A - n×n 邻接矩阵（可为 0/1 或权重，只要>0表示连边即可）
%
% 输出：
%   pair_count - 对称节点对数量
%   sym_groups - 对称节点分组（元胞数组，每组内节点邻居相同）

    % === 预处理 ===
    n = size(A,1);
    A = double(A > 0);     % 转换为0/1矩阵
    A(1:n+1:end) = 0;      % 去掉自环

    visited = false(n,1);
    sym_groups = {};
    pair_count = 0;

    % === 主循环：逐行比较邻居向量 ===
    for i = 1:n
        if visited(i), continue; end
        % 找出邻接行完全相同的节点
        same = all(A == A(i,:), 2);
        idx = find(same);
        if numel(idx) > 1
            sym_groups{end+1} = idx; %#ok<AGROW>
            pair_count = pair_count + nchoosek(numel(idx),2);
        end
        visited(idx) = true;
    end
end



function visualize_current_flow(A, i, j,lsg, xy)
% visualize_current_flow(A, i, j, xy)
% xy 可选参数，用于固定节点位置
%
%  - 若 xy 未提供，则自动用 force 布局生成一次；
%  - 若提供 xy（n×2 矩阵），则使用固定坐标。

    n = size(A,1);

    % ---- Laplacian ----
    L = diag(sum(A,2)) - A;
    b = zeros(n,1); b(i)=1; b(j)=-1;
    ref = n;
    L_red = L(1:end-1, 1:end-1);
    b_red = b(1:end-1);
    v = zeros(n,1);
    v(1:end-1) = L_red \ b_red;
    I_edges = A .* (v - v');
    mask = (A > 0) & (abs(I_edges) > 1e-12);

    % ---- 构建图 ----
    G = graph(A);

    % ---- 若未提供坐标，则第一次运行生成并保存 ----
    if nargin < 5 || isempty(xy)
        figure;
        htemp = plot(G, 'Layout', 'force');
        xy = [htemp.XData', htemp.YData'];
        close; % 关闭临时图
    end

    % ---- 绘图 ----
    figure;
    h = plot(G, 'XData', xy(:,1), 'YData', xy(:,2), ...
        'NodeColor', 'k', 'MarkerSize', 6, 'EdgeColor', [0.7 0.7 0.7]);
    hold on;

    [rows, cols] = find(triu(mask));
    highlight(h, rows, cols, 'EdgeColor', 'r', 'LineWidth', 2);
    highlight(h, i, 'NodeColor', 'g', 'MarkerSize', 8);
    highlight(h, j, 'NodeColor', 'm', 'MarkerSize', 8);

    title(sprintf('Current flow from node %d → node %d, link number %d', i, j,lsg), 'FontSize', 12);
    hold off;
end