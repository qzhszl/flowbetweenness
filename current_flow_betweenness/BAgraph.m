function G = BAgraph(n, m0, m)
%BAGRAPH_WEIGHTED 生成带权 BA 网络
%   n  : 节点总数
%   m0 : 初始完全图节点数
%   m  : 每次新加入的节点与已有节点连边数
%
% 输出:
%   G : graph 对象，Edges.Weight 为 [0,1] 均匀分布随机数

% ---- 初始完全图 ----
A = ones(m0) - eye(m0);
G = graph(A);

% 初始边赋随机权重
numE = numedges(G);
G.Edges.Weight = rand(numE,1);

% 度向量
deg = degree(G);

% ---- 逐步加入新节点 ----
for newNode = (m0+1):n
    % 按度分布概率选择 m 个已有节点
    prob = deg / sum(deg);
    targets = randsample(newNode-1, m, true, prob);
    
    % 加边
    G = addedge(G, repmat(newNode,1,m), targets,rand);
    
    % 更新度
    deg = degree(G);
end
end

