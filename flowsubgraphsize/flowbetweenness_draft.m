clear,clc
    
n = 50; 
p = 0.03;
A = rand(n) < p;
A = triu(A,1); A = A + A';

src = 1; 
tgt = 2;

[num_nodes, num_edges, nodes_flow, edges_flow] = count_flow_plot(A, src, tgt, 1e-12);

fprintf('有电流的节点数: %d\n', num_nodes);
fprintf('有电流的边数  : %d\n', num_edges);

G = graph(A);
[flowsubgraphlink,lsg] = flowsubgraph(G,src,tgt)
if lsg ~= 0
    FB = tril(A);
    FB(find(FB)) = flowsubgraphlink;
    FB = FB+FB.';
    flowsubgraphnode = sum(abs(FB)).';
    flowsubgraphnode = find(abs(flowsubgraphnode)>0.00000001);
    nodesize = size(flowsubgraphnode,1)
end

function [num_nodes_flow, num_edges_flow, nodes_flow, edges_flow] = count_flow_plot(A, src, tgt, tol)
% COUNT_FLOW_ALL_COMPONENTS
% 计算每个连通分量内部电流流经的节点和边
%
% 输入:
%   A    : 邻接矩阵 (n x n, 对称)
%   src  : 源节点编号
%   tgt  : 汇节点编号
%   tol  : 电流阈值 (默认 1e-12)
%
% 输出:
%   num_nodes_flow : 有电流流经节点数
%   num_edges_flow : 有电流流经边数
%   nodes_flow     : n x 1 逻辑向量
%   edges_flow     : m x 2 数组，边为节点编号对

if nargin < 4, tol = 1e-12; end

n = size(A,1);
nodes_flow = false(n,1);
edges_flow = [];

G = graph(A);
bins = conncomp(G);

% 如果源节点和汇节点不连通
if bins(src) ~= bins(tgt)
    num_nodes_flow = 0;
    num_edges_flow = 0;
    return
end

% 只在源-汇所在的连通分量内计算
comp_id = bins(src);
comp_idx = find(bins == comp_id);
A_sub = A(comp_idx, comp_idx);

% 映射源和汇到子矩阵
src_sub = find(comp_idx == src);
tgt_sub = find(comp_idx == tgt);

% 拉普拉斯矩阵
d = sum(A_sub,2);
L = diag(d) - A_sub;

% 注入电流向量
b = zeros(length(comp_idx),1);
b(src_sub) = 1;
b(tgt_sub) = -1;

% 使用伪逆求解电位
v = pinv(L) * b;

% 边电流
[ii, jj, ww] = find(triu(A_sub));
Ivals = (v(ii) - v(jj)) .* ww;
edges_mask = abs(Ivals) > tol;
edges_flow_sub = [ii, jj];
edges_flow_sub = edges_flow_sub(edges_mask,:);

% 节点流量
nodes_flow_sub = false(length(comp_idx),1);
for u = 1:length(comp_idx)
    neigh = find(A_sub(u,:));
    if any(abs(v(u) - v(neigh)) .* A_sub(u,neigh) > tol)
        nodes_flow_sub(u) = true;
    end
end

% 转换回原节点编号
nodes_flow(comp_idx) = nodes_flow_sub;
edges_flow = comp_idx(edges_flow_sub);

num_nodes_flow = sum(nodes_flow_sub);
num_edges_flow = size(edges_flow,1);

% 可视化
figure;
p_layout = plot(G,'Layout','force','NodeColor',[0.7 0.7 0.7],'EdgeColor',[0.8 0.8 0.8]);
hold on
highlight(p_layout, find(nodes_flow), 'NodeColor','r','MarkerSize',6);
for k = 1:size(edges_flow,1)
    highlight(p_layout, edges_flow(k,1), edges_flow(k,2), 'EdgeColor','r','LineWidth',1.5);
end
highlight(p_layout, src, 'NodeColor','g','MarkerSize',8);
highlight(p_layout, tgt, 'NodeColor','b','MarkerSize',8);
title('电流流经的节点和边（每个分量独立计算）');
end


