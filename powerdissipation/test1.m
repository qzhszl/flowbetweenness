function [total_power_givenst, EdgeEnergy_flow] = compute_power_dissipation_eachlink_new(G, s, t) 
% 输入:
%   G: graph 对象, G.Edges.Weight = 导纳 g
%   s, t: 节点编号 (注入/抽出电流)

n = numnodes(G);
m = numedges(G);
g = G.Edges.Weight;   % 边导纳

% 构造稀疏拉普拉斯矩阵
Adj = adjacency(G, 'weighted');
L = spdiags(sum(Adj,2), 0, n, n) - Adj;

% 构造 b 向量
b = sparse(n,1);
b(s) = 1;
b(t) = -1;

% 解电势 (固定最后一个点电势=0, 避免奇异性)
L_mod = L(1:end-1, 1:end-1);
b_mod = b(1:end-1);
v = zeros(n,1);
v(1:end-1) = L_mod \ b_mod;

% 各边电势差 (一次向量化)
u = G.Edges.EndNodes(:,1);
vtx = G.Edges.EndNodes(:,2);
dv = v(u) - v(vtx);

% 每条边的能量 g * (Δv)^2
EdgeEnergy_flow = g .* (dv.^2);

% 总能量
total_power_givenst = sum(EdgeEnergy_flow);
end
