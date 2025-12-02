% keep_lcc_and_renumber.m
% 读取 edges_final.txt（格式：node1 node2 weight）
% 保留最大连通分量（LCC），将该分量节点重命名为 1..n
% 输出到 edges_lcc_renumbered.txt（三列：newNode1 newNode2 weight）

clear; clc;

inputFile  = 'D:\\data\\flow betweenness\\IERP\\realnetwork\\edges_final.txt';         % 原始文件（由你前一步生成）
outputFile = 'edges_lcc_renumbered.txt';% 输出文件（重编号后的最大连通分量）

% -------- 读取文件（支持多种空白分隔格式） --------
if ~isfile(inputFile)
    error('找不到输入文件：%s', inputFile);
end

% 尝试用 readmatrix（MATLAB R2019b+），回退到 textscan 如果失败
try
    M = readmatrix(inputFile);
    % readmatrix 可能会返回 Nx3 的数值矩阵
    if size(M,2) < 3
        error('readmatrix 未能读取三列数值，切换到 textscan 读取。');
    end
    u_raw = M(:,1);
    v_raw = M(:,2);
    w_raw = M(:,3);
catch
    % 更稳健的文本读取：逐行 parse
    fid = fopen(inputFile,'r');
    u_raw = [];
    v_raw = [];
    w_raw = [];
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if isempty(line), continue; end
        % 忽略注释行（以 # 或 % 开头）
        if startsWith(line,'#') || startsWith(line,'%')
            continue;
        end
        toks = regexp(line, '\s+', 'split');
        if numel(toks) < 3
            continue;
        end
        % 尝试把前三项解析为数值
        a = str2double(toks{1});
        b = str2double(toks{2});
        c = str2double(toks{3});
        if isnan(a) || isnan(b) || isnan(c)
            continue;
        end
        u_raw(end+1,1) = a; %#ok<SAGROW>
        v_raw(end+1,1) = b; %#ok<SAGROW>
        w_raw(end+1,1) = c; %#ok<SAGROW>
    end
    fclose(fid);
end

% 检查是否成功读取
if isempty(u_raw)
    error('未能从 %s 读取任何边数据，请检查文件格式。', inputFile);
end

% -------- 将原始节点 ID 映射到连续索引 1..N --------
origIDs = unique([u_raw; v_raw]);   % 原始节点 ID（可能不连续）
N_all = numel(origIDs);
% 建立 id -> index 映射
id2idx = containers.Map('KeyType','double','ValueType','double');
for k = 1:N_all
    id2idx(origIDs(k)) = k;
end

% 将边转换为索引形式
E = numel(u_raw);
u_idx = zeros(E,1);
v_idx = zeros(E,1);
for e = 1:E
    u_idx(e) = id2idx(u_raw(e));
    v_idx(e) = id2idx(v_raw(e));
end
w_vals = w_raw;

% 去自环（虽然前一步已经去过，但再保险）
valid = u_idx ~= v_idx;
u_idx = u_idx(valid);
v_idx = v_idx(valid);
w_vals = w_vals(valid);

% 因为原图无重边并且无向，但为稳健起见合并任何重边（权重相加），并确保无向性
% 我们只在上三角记录边（u<v）
a = min(u_idx, v_idx);
b = max(u_idx, v_idx);
keyStrings = strcat(string(a), '_', string(b));
pairMap = containers.Map('KeyType','char','ValueType','double');
for k = 1:numel(a)
    key = keyStrings{k};
    if isKey(pairMap, key)
        pairMap(key) = pairMap(key) + w_vals(k);
    else
        pairMap(key) = w_vals(k);
    end
end

% 从 pairMap 恢复合并后的边表（基于连续索引 1..N_all）
keysList = keys(pairMap);
M_final = numel(keysList);
u_comb = zeros(M_final,1);
v_comb = zeros(M_final,1);
w_comb = zeros(M_final,1);
for k = 1:M_final
    key = keysList{k};
    nums = sscanf(key, '%d_%d');
    u_comb(k) = nums(1);
    v_comb(k) = nums(2);
    w_comb(k) = pairMap(key);
end

% -------- 构建图并找连通分量 --------
G_all = graph(u_comb, v_comb, w_comb); % 无向带权图
comp = conncomp(G_all);                % 每个节点所属连通分量编号
comp_ids = unique(comp);
% 计算每个分量的大小
counts = histcounts(comp, [comp_ids, comp_ids(end)+1]);
[~, idx_max] = max(counts);
lcc_id = comp_ids(idx_max);

% 找到属于 LCC 的节点（这些节点是 G_all 的节点索引）
nodes_in_lcc = find(comp == lcc_id);
n_lcc = numel(nodes_in_lcc);
fprintf('原图节点数：%d，边数：%d\n', numnodes(G_all), numedges(G_all));
fprintf('最大连通分量节点数：%d\n', n_lcc);

% -------- 从 LCC 提取子图并重编号 1..n_lcc --------
subG = subgraph(G_all, nodes_in_lcc);  % subG 节点索引为 1..n_lcc 在 subG 内部
% 注意：subG 的节点索引对应于 nodes_in_lcc 的顺序，
% 也就是说 subG 中的节点 i 对应原图索引 nodes_in_lcc(i)

% 获取上三角非零元素以列出每条无向边一次
A_sub = adjacency(subG,'weighted');
[i_sub, j_sub, w_sub] = find(triu(A_sub)); % i_sub,j_sub 都是 1..n_lcc

% -------- 将结果写入新文件（节点重编号为 1..n_lcc） --------
fid = fopen(outputFile,'w');
if fid == -1
    error('无法创建输出文件：%s', outputFile);
end

% 写入每条边： newNode1 newNode2 weight
for k = 1:numel(i_sub)
    fprintf(fid, '%d\t%d\t%g\n', i_sub(k), j_sub(k), w_sub(k));
end
fclose(fid);

fprintf('已保存最大连通分量（重编号）到 %s 。共有 %d 个节点，%d 条边。\n', ...
    outputFile, n_lcc, numel(i_sub));
