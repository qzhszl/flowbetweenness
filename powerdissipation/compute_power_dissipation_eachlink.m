function compute_power_dissipation_eachlink(n, p, simu_times,weighted)
%SIMULATE_ER_POWER 重复生成 ER 图，计算所有节点对的能量消耗
%
% 输入:
%   n          节点数
%   p          ER 图连边概率
%   simu_times 重复次数
%
% 输出:
%   Results_SP   每次仿真 shortest-path 模型能量
%   Results_Flow 每次仿真 current-flow 模型能量
%
% 要求：已存在两个函数
%   [power_givenst_SP, EdgeEnergy_SP] = compute_path_power_dissipation(G_path,s,t)
%   [power_givenst_flow, EdgeEnergy_flow] = compute_flownetwork_power_dissipation(G,s,t)
%
% 其中 G_path 的边权重 = G 的边权重的倒数

results = struct();
filefolder_name = "D:\\data\\flow betweenness\\";

for k = 1:simu_times
    fprintf('rep %d/%d\n', k, simu_times);
    % ---- 生成 ER 图 ----
    A =  GenerateERfast(n,p,weighted);
    G = graph(A);  % 无向图
    ind_conn = 1;
    while(ind_conn) % Until a connected graph is created
        A = GenerateERfast(n,p,weighted);       
        G = graph(A);
        if(sum(conncomp(G)) == n)
            ind_conn = 0;
        end
    end

%     figure;
%     plot(G,'EdgeLabel',G.Edges.Weight,'NodeColor',[0.8500 0.3250 0.0980], ...
% 'EdgeAlpha',0.5,'LineWidth',1,'MarkerSize',7,'EdgeLabelColor',[0 0.4470 0.7410],'NodeFontSize',10);

    % 随机生成边权重 (导纳)，避免 0 权
    
    % G_path 的边权重 = 1/g
    % figure;
    G_path = graph(G.Edges.EndNodes(:,1), G.Edges.EndNodes(:,2), 1./G.Edges.Weight,numnodes(G));
%     plot(G_path,'EdgeLabel',G_path.Edges.Weight,'NodeColor',[0.8500 0.3250 0.0980], ...
% 'EdgeAlpha',0.5,'LineWidth',1,'MarkerSize',7,'EdgeLabelColor',[0 0.4470 0.7410],'NodeFontSize',10);

    % ---- 对所有节点对 (i,j) 计算 ----
    n = numnodes(G);
    total_SP   = zeros(nchoosek(n,2),1);   
    total_Flow = zeros(nchoosek(n,2),1);
    linkP_SP = zeros(numedges(G),1);
    linkP_Flow = zeros(numedges(G),1);
    count = 0;
    for s = 1:n
        for t = (s+1):n
            count = count+1;
            % shortest path
            
            [total_SP_fornodepair, link_SP_fornodepair,connected_flag] = compute_path_power_dissipation(G_path,s,t);
            
            if connected_flag==1
                total_SP(count) = total_SP_fornodepair; 
                linkP_SP = linkP_SP+link_SP_fornodepair;
                % current flow
                [total_Flow_fornodepair, link_Flow_fornodepair] = compute_flownetwork_power_dissipation(G,s,t);
                total_Flow(count) = total_Flow_fornodepair;
                linkP_Flow = linkP_Flow+link_Flow_fornodepair;        
            end
        end
    end
    
    % 保存结果
    results(k).edges = G.Edges;  % [i j w]
    results(k).total_Flow = total_Flow;
    results(k).total_SP = total_SP;
    results(k).linkP_Flow = linkP_Flow;
    results(k).linkP_SP = linkP_SP;
end
if weighted ==0
    filename = sprintf('power_dissipation_N%dp%.2fER_unweighted.mat',n,p);
else
    filename = sprintf('power_dissipation_N%dp%.2fER.mat',n,p);
end
% 保存到文件
save(filefolder_name+filename,'results');
disp("mission_completed")
end


