function [total_power_givenst,EdgeEnergy_flow] = compute_flownetwork_power_dissipation(G, s, t)
%EDGEENERGYCURRENTFLOW 计算电流定律模型下的能量消耗 (权重=导纳)
%
% 输入:
%   G (graph 对象, G.Edges.Weight = 导纳 g)
%   s, t: 节点编号 (注入/抽出电流)
%
% 输出:
%   EdgeEnergy_flow 包含每条边的能量消耗
%   total_power_givenst以及总能量 (等于有效电阻)

n = numnodes(G);
m = numedges(G);
g = G.Edges.Weight;    % 导纳
% r = 1./g;              % 电阻

% 构造拉普拉斯矩阵
Adj = adjacency(G,g);
L = diag(sum(Adj,2)) - Adj;

% 伪逆解电势
Lplus = pinv(full(L));

b = zeros(n,1); b(s)=1; b(t)=-1;
v = Lplus*b;

% 计算每条边的能量
EdgeEnergy_flow = zeros(m,1);
for e = 1:m
    u = G.Edges.EndNodes(e,1);
    vtx = G.Edges.EndNodes(e,2);
    dv = v(u) - v(vtx);
    EdgeEnergy_flow(e) = g(e) * (dv^2);  % g * (Δv)^2
end

% 总能量 = 有效电阻
total_power_givenst = sum(EdgeEnergy_flow);
end


