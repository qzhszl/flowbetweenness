clear,clc
% this.m inverstigate the size of the flow subgraph(nodes) in ER graph
% We need 

N = 6;
pc= log(N)/N;
ave_degree = 0.2:0.2:N-1;
% ave_degree_2 = 5:10;

% ave_degree = 1.2:0.2:4;
% ave_degree_2 = 5:10;
% ave_degree = [ave_degree,ave_degree_2];

ave_degree = [5]
p_list = ave_degree/(N-1);
siumtimes = 1;
xy = [0.802385749219495	1.59761712128446
1.76031580715563	0.102069992210448
0.984548013597186	-1.49701944025127
-1.78749905874878	-0.106709857894361
-0.787506726261186	-1.56875592122334
-0.972243784962346	1.47279810587407];


for p=p_list(1:length(p_list))
    disp(p)
    nodep=0;

    for i = 1:siumtimes
%         A = GenerateERfast(N,p,0);

        A = [0	1	1	1	1	1
            1	0	1	1	1	0
            1	1	0	1	1	1
            1	1	1	0	1	1
            1	1	1	1	0	1
            1	0	1	1	1	0];

        real_ave_degree = mean(sum(A));
        real_ave_degree_list(i) = real_ave_degree;
        
        G = graph(A);

        % h = plot(G,'NodeColor',[0.8500 0.3250 0.0980], ...
        % 'EdgeAlpha',0.5,'LineWidth',1,'MarkerSize',7,'EdgeLabelColor',[0 0.4470 0.7410],'NodeFontSize',10);
        % hold on

        Linknumlist(i) = numedges(G);
        % [max_size, second_max_size] = getLargestComponentsSize(G);
        % Lcc_list(i)  = max_size;
        % sLcc_list(i) = second_max_size;
        nodei = 1;
        nodej = 5;

        [flowsubgraphlink,lsg] = flowsubgraph(G,nodei,nodej);

%         [flowsubgraphlink, lsg] = compute_edge_currents(A, nodei, nodej)
        
        flow_subgraph_linksize_list(i) = lsg

        visualize_current_flow(A, nodei, nodej,lsg,xy)
        h = findobj(gca,'Type','GraphPlot');
        xy = [h.XData' h.YData'];

        
        
        if lsg ~= 0
            FB = tril(A);
            FB(find(FB)) = flowsubgraphlink;
            FB = FB+FB.';
            flowsubgraphnode = sum(abs(FB)).';
            flowsubgraphnode = find(abs(flowsubgraphnode)>0.00000001);
            nodesize = size(flowsubgraphnode,1);
            flow_subgraph_size_list(i) = nodesize
            
        end

               
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



function [edge_currents, num_active_edges] = compute_edge_currents(A, i, j)
% 计算从节点 i 注入 1A 电流、节点 j 流出 1A 时的边电流
% 并统计存在电流的边数
%
% 输入:
%   A - 邻接矩阵 (n×n)，表示导纳 (1/电阻)，A_ij>0表示有边
%   i, j - 源节点与汇节点编号 (1-based)
%
% 输出:
%   edge_currents - 每条边上的电流矩阵（与 A 相同大小, 对称）
%   num_active_edges - 电流不为 0 的边数量（无向边计一次）

    n = size(A, 1);

    % ---- 1. 图拉普拉斯矩阵 ----
    L = diag(sum(A, 2)) - A;

    % ---- 2. 电流注入向量 ----
    b = zeros(n, 1);
    b(i) = 1;
    b(j) = -1;

    % ---- 3. 去掉一个参考节点（防止奇异） ----
    ref = n;
    L_red = L(1:end-1, 1:end-1);
    b_red = b(1:end-1);

    % ---- 4. 求解节点电位 ----
    v = zeros(n, 1);
    v(1:end-1) = L_red \ b_red;   % v(ref)=0

    % ---- 5. 欧姆定律求电流 ----
    edge_currents = A .* (v - v');  % I_ij = A_ij*(v_i - v_j)

    % ---- 6. 统计存在电流的边数（无向边计一次）----
    abs_currents = abs(edge_currents);
    tol = 1e-12;  % 数值容差，小于此视为0
    mask = (A > 0) & (abs_currents > tol);
    
    % 只统计上三角部分（避免重复）
    num_active_edges = nnz(triu(mask));

end

 



function [max_size, second_max_size] = getLargestComponentsSize(G)
% 输入：
%   G - graph 或 digraph 对象
%
% 输出：
%   max_size - 最大联通子图的节点数
%   second_max_size - 第二大联通子图的节点数（如果只有一个，则为 0）

    % 获取每个节点所属的联通分量编号
    comp_ids = conncomp(G);

    % 统计每个联通分量的节点数量
    comp_sizes = histcounts(comp_ids, 1:(max(comp_ids)+1));

    % 如果只存在一个分量，返回第二个为0
    if isscalar(comp_sizes)
        max_size = comp_sizes(1);
        second_max_size = 0;
    else
        % 排序，找到最大和第二大
        sorted_sizes = sort(comp_sizes, 'descend');
        max_size = sorted_sizes(1);
        second_max_size = sorted_sizes(2);
    end
end
