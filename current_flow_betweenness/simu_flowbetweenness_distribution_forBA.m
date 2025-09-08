function simu_flowbetweenness_distribution_forBA()
% experiment_metrics.m
% 假设 current_flow_betweenness.m 已经在路径中

% rng(1); % 固定随机种子以保证可重复

% parameter
model = 'BA';  % graph type: 'ER' or 'BA'
R = 1;        % 重复次数

% 保存容器
results = struct();
% filefolder_name = "D:\\data\\flow betweenness\\";

input_BAname = 'BAnetworkN1000m3.txt';

for r = 1:R
    fprintf('rep %d/%d\n', r, R);
    % --- 生成网络 ---
    if strcmp(model,'ER')
        % Step 1: Generate a connected network
        ind_conn = 1;
        A = GenerateERfast(N,p,1);
        while(ind_conn) % Until a connected graph is created
            A = GenerateERfast(N,p,1);       
            G = graph(A);
            if(sum(conncomp(G)) == N)
                ind_conn = 0;
            end
        end
        % --- 构造 graph 对象 ---
        G = graph(A);
    elseif strcmp(model,'BA')
        data = readmatrix(input_BAname);  % 每行: u v w
        % 提取三列
        s = data(:,1)+1;    % 起点
        t = data(:,2)+1;    % 终点
        w = data(:,3);    % 权重
        
        % 构造无向图（如果是有向图，就用 digraph）
        G = graph(s, t, w);
    end
    A = adjacency(G,"weighted");
    
    % --- 计算指标 ---
    % 1. & 2. current-flow betweenness
    [nodeCFB, edgeCFB, edges] = current_flow_betweenness(A);

    % 3. node shortest-path betweenness
    nodeSPB = centrality(G,'betweenness','Cost',G.Edges.Weight);
    
    % 4. edge shortest-path betweenness 
    edgeSPB = edgeBetweennessCentrality(G); 

    % 5. degree
    degree_unweighted = degree(G);

    % 6. weighted degree
    degree_weighted = sum(A,2);

    % --- 存储 ---
    results(r).A = A;
    results(r).edges = edges;  % [i j w]
    results(r).nodeCFB = nodeCFB;
    results(r).edgeCFB = edgeCFB;
    results(r).nodeSPB = nodeSPB;
    results(r).edgeSPB = edgeSPB;
    results(r).deg = degree_unweighted;
    results(r).wdeg = degree_weighted;
end
filename = sprintf('bet_cbet_degree_.mat');
[~, insertStrNoTxt, ~] = fileparts(input_BAname);
% 找到 '_' 的位置
pos = strfind(filename,'_');
filename = [input_BAname(1:pos), insertStrNoTxt, input_BAname(pos+1:end)];
% 保存到文件
save(filename,'results');
disp("mission_completed")
end