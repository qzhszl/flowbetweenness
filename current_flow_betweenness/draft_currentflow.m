clear,clc
clear,clc
% G = BAgraph(1000, 4, 3);
data = readmatrix('D:\\data\\flow betweenness\\BA_network\\BAnetworkN5000m3.txt');  % 每行: u v w

% 提取三列
s = data(:,1)+1;    % 起点
t = data(:,2)+1;    % 终点
w = data(:,3);    % 权重

% 构造无向图（如果是有向图，就用 digraph）
G = graph(s, t, w);
A = adjacency(G,"weighted")

[row,col,v]=find(G.Edges.Weight<=0)


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




