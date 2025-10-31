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

        A = [0	0	1	1	1	1
            0	0	1	1	1	1
            1	1	0	1	1	1
            1	1	1	0	1	0
            1	1	1	1	0	1
            1	1	1	0	1	0];

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

        % [Ac, Gc, p] = complement_graph(A);
        % [L, GL, p, edgeTable] = line_graph_from_adj(A);
        
        Q = sum(A)-A
        [V,D] = eigs(Q)

        
        
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


function sym_groups = find_symmetric_nodes(A)
% 返回所有在结构上对称的节点组
% 输入：A - 无向邻接矩阵 (n×n)
% 输出：sym_groups - 元胞数组，每个元素是一组对称节点

    n = size(A,1);
    A = A + A';        % 若是有向图，可注释此行
    A(A>0) = 1;        % 无权化
    
    sym_groups = {};
    visited = false(n,1);
    
    for i = 1:n
        if visited(i), continue; end
        % 找所有邻接向量完全相同的节点
        same = all(A == A(i,:), 2);
        idx = find(same);
        if numel(idx) > 1
            sym_groups{end+1} = idx; %#ok<AGROW>
        end
        visited(idx) = true;
    end
end



function [L, GL, p, edgeTable] = line_graph_from_adj(A, nodeNames)
% LINE_GRAPH_FROM_ADJ  由邻接矩阵 A 构造线图的邻接矩阵 L，并绘图。
% 假设 A 为无向简单图；函数会做合理化处理（同上）。
%
% 输入：
%   A          (n×n) 邻接矩阵
%   nodeNames  (可选) 1×n 的字符串/字符元胞数组——用于原图的结点命名，
%                     线图结点将按“u-v”命名（对应原图的一条边）
%
% 输出：
%   L         (m×m) 线图的邻接矩阵（逻辑），m 为原图边数
%   GL        MATLAB graph 对象（线图）
%   p         plot 句柄
%   edgeTable 原图的边表（table），含 EndNodes，便于溯源

    arguments
        A {mustBeNumericOrLogical, mustBeSquare(A)}
        nodeNames = []
    end

    n = size(A,1);

    % —— 规范化 A —— %
    A = double(A ~= 0);
    A(1:n+1:end) = 0;
    A = triu(A,1);
    A = A + A.';

    % 原图（用于拿边表）
    if isempty(nodeNames)
        G = graph(A);
    else
        G = graph(A, string(nodeNames));
    end

    edgeTable = G.Edges;     % table, 变量名通常为 EndNodes
    m = numedges(G);
    if m == 0
        L  = false(0,0);
        GL = graph();
        warning('原图无边，线图为空。');
        p = [];
        return;
    end

    % —— 构造点-边关联矩阵 B —— %
    % B(i,e) = 1 若结点 i 与边 e 关联
    B = zeros(n, m);
    ends = G.Edges.EndNodes;  % m×2，存储每条边的两个端点名（字符串）
    % 将端点名映射回索引
    if isempty(nodeNames)
        % 节点默认编号为 1..n 的字符串
        name2idx = containers.Map(string(1:n), 1:n);
    else
        names = string(G.Nodes.Name);
        name2idx = containers.Map(names, 1:n);
    end
    for e = 1:m
        u = name2idx(string(ends(e,1)));
        v = name2idx(string(ends(e,2)));
        B(u,e) = 1;
        B(v,e) = 1;
    end

    % —— 线图邻接 —— %
    Lw = B.' * B;            % 对角线上=2；非对角=共享端点的次数（0或1）
    Lw = Lw - 2*eye(m);      % 去掉对角 2
    L  = Lw > 0;             % 二值化
    L  = logical(L);

    % —— 线图结点命名（用“u-v”） —— %
    ln = strings(m,1);
    for e = 1:m
        ln(e) = ends(e,1) + "-" + ends(e,2);
    end

    GL = graph(L, ln);

    % —— 绘图 —— %
    figure('Name','Line Graph');
    p = plot(GL, 'Layout','force');
    title('Line Graph');

end





function [Ac, Gc, p] = complement_graph(A, nodeNames)
% COMPLEMENT_GRAPH  由邻接矩阵 A 生成补图的邻接矩阵 Ac，并绘图。
% 假设 A 表示无向简单图（无自环、无多重边）。函数会做合理化处理：
%   - 去掉对角元
%   - 二值化
%   - 对称化
%
% 输入：
%   A          (n×n) 邻接矩阵（数值或逻辑）
%   nodeNames  (可选) 1×n 的字符串/字符元胞数组，用作结点标签
%
% 输出：
%   Ac  (n×n) 补图的邻接矩阵（逻辑）
%   Gc  MATLAB graph 对象
%   p   plot 句柄

    arguments
        A {mustBeNumericOrLogical, mustBeSquare(A)}
        nodeNames = []
    end

    n = size(A,1);

    % —— 规范化 A 为无向简单图 —— %
    A = double(A ~= 0);
    A(1:n+1:end) = 0;           % 去对角元（自环）
    A = triu(A,1);              % 只保留上三角，去掉可能的多重与方向
    A = A + A.';                % 对称化
    
    % —— 生成补图 —— %
    Ac = ~A & ~eye(n);          % 非边且非自环
    Ac = logical(Ac);

    % —— 绘图 —— %
    if isempty(nodeNames); Gc = graph(Ac);
    else;                   Gc = graph(Ac, string(nodeNames));
    end
    figure('Name','Complement Graph');
    p = plot(Gc, 'Layout','force');
    title('Complement Graph');

end

function mustBeSquare(A)
    if size(A,1) ~= size(A,2)
        error('A 必须是方阵。');
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
