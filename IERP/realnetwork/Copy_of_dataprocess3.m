% keep_lcc_and_renumber.m
% 读取 edges_final.txt（格式：node1 node2 weight）
% 保留最大连通分量（LCC），将该分量节点重命名为 1..n
% 输出到 edges_lcc_renumbered.txt（三列：newNode1 newNode2 weight）

clear; clc;

inputFile  = 'D:\\data\\flow betweenness\\IERP\\realnetwork\\Newspain_clean.txt';         % 原始文件（由你前一步生成）

% -------- 读取文件（支持多种空白分隔格式） --------
if ~isfile(inputFile)
    error('找不到输入文件：%s', inputFile);
end

M = readmatrix(inputFile);
% readmatrix 可能会返回 Nx3 的数值矩阵
u_raw = M(:,1);
v_raw = M(:,2);
w_raw = M(:,3);

G_all = graph(u_raw, v_raw, w_raw); % 无向带权图
A = full(G_all.adjacency("weighted"))
plot(G_all)