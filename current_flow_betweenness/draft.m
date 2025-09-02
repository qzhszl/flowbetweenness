%% ER graph betweenness vs current-flow betweenness
% 参数设置
n = 100;       % 节点数
p = 0.05;      % 边的概率 (ER G(n,p))

% 生成ER随机图
A = rand(n) < p;        % adjacency matrix
A = triu(A,1);          % 只保留上三角，避免双向重复
A = A + A';             % 转为无向
G = graph(A);

% 计算 betweenness centrality
betw = centrality(G,'betweenness');

% --- current-flow betweenness ---
L = diag(sum(A,2)) - A;   % Laplacian
cf_betw = zeros(n,1);

% 遍历所有节点对 (s,t)
pairs = nchoosek(1:n,2);
for k = 1:size(pairs,1)
    s = pairs(k,1); t = pairs(k,2);
    
    % 电流注入向量 b
    b = zeros(n,1);
    b(s) = 1; b(t) = -1;
    
    % 解电位 Lv = b（伪逆，因为L奇异）
    v = pinv(L)*b;
    
    % 计算边电流
    for u = 1:n
        for vtx = u+1:n
            if A(u,vtx) == 1
                I = abs(v(u)-v(vtx)); % conductance=1
                cf_betw(u)   = cf_betw(u) + I/2; % 平分到两个节点
                cf_betw(vtx) = cf_betw(vtx) + I/2;
            end
        end
    end
end

% 归一化
cf_betw = cf_betw / size(pairs,1);

% 可视化分布
figure;
subplot(1,2,1);
histogram(betw,20,'FaceColor',[0.2 0.6 0.8]);
title('Betweenness Centrality');
xlabel('Value'); ylabel('Frequency');

subplot(1,2,2);
histogram(cf_betw,20,'FaceColor',[0.8 0.4 0.4]);
title('Current-Flow Betweenness');
xlabel('Value'); ylabel('Frequency');

sgtitle(sprintf('ER Graph G(%d, %.2f)',n,p));
