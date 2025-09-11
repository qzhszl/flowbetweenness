function simu_flowbetweenness_distribution(N,p)
% experiment_metrics.m
% 假设 current_flow_betweenness.m 已经在路径中

% rng(1); % 固定随机种子以保证可重复

% parameter
model = 'ER';  % graph type: 'ER' or 'BA'
n = N;       % 节点数
p = p;      % ER 参数
R = 1000;        % 重复次数

% 保存容器
results = struct();
filefolder_name = "D:\\data\\flow betweenness\\";

for r = 1:R
    fprintf('rep %d/%d\n', r, R);
    % --- 生成网络 ---
    if strcmp(model,'ER')
        % Step 1: Generate a connected network
        ind_conn = 1;
        A = GenerateERfast(N,p,0);
        while(ind_conn) % Until a connected graph is created
            A = GenerateERfast(N,p,0);       
            G = graph(A);
            if(sum(conncomp(G)) == N)
                ind_conn = 0;
            end
        end
    elseif strcmp(model,'BA')
        A = zeros(n);
        A(1:m+1,1:m+1) = 1 - eye(m+1);
        deg = sum(A,2);
        for new = (m+2):n
            probs = deg(1:new-1) / sum(deg(1:new-1));
            idx = randsample(new-1, m, true, probs);
            for t = idx(:)'
                A(new,t) = 1;
                A(t,new) = 1;
            end
            deg = sum(A(1:new,1:new),2);
        end
    end

    % --- 构造 graph 对象 ---
    G = graph(A);

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
filename = sprintf('bet_cbet_degree_N%dp%.2fER_unweighted.mat',N,p)
% 保存到文件
save(filefolder_name+filename,'results');
disp("mission_completed")
end