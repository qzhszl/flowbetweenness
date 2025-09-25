clear,clc
% this.m inverstigate the size of the flow subgraph(nodes) in ER graph
% We need 

N = 1000;
pc= log(N)/N;
ave_degree = 1:0.2:4.9;
ave_degree_2 = 5:10;

% ave_degree = 1.2:0.2:4;
% ave_degree_2 = 5:10;
% ave_degree = [ave_degree,ave_degree_2];
p_list = ave_degree_2/(N-1);
   

% nodep_list =zeros(length(p_list),1);
avesize_fsg_list =zeros(length(p_list),1);
stdsize_fsg_list =zeros(length(p_list),1);
% lcc_diffp = zeros(length(p_list),1);
% slcc_diffp = zeros(length(p_list),1);
% p_ave = zeros(length(p_list),1);
% AnalysisSolution = zeros(length(p_list),1);

count=0;
siumtimes = 1000;
filefolder_name = "D:\\data\\flow betweenness\\sizeofflowsubgraph\\new";


for p=p_list(1:length(p_list))
    disp(p)
    count = count+1;
    nodep=0;
    plist = zeros(siumtimes,1);
    Linknumlist = zeros(siumtimes,1);
    % Lcc_list = zeros(siumtimes,1);
    % sLcc_list = zeros(siumtimes,1);
    flow_subgraph_size_list = zeros(siumtimes,1);
    flow_subgraph_linksize_list = zeros(siumtimes,1);
    real_ave_degree_list = zeros(siumtimes,1);
    for i = 1:siumtimes
        if mod(i,100)==0
            disp(i/100)
        end     
        A = rand(N,N) < p;
        A = triu(A,1);
        A = A + A';
     
   
        real_ave_degree = mean(sum(A));
        real_ave_degree_list(i) = real_ave_degree;
        
        G = graph(A);

    %     plot(G,'NodeColor',[0.8500 0.3250 0.0980], ...
    % 'EdgeAlpha',0.5,'LineWidth',1,'MarkerSize',7,'EdgeLabelColor',[0 0.4470 0.7410],'NodeFontSize',10);

        Linknumlist(i) = numedges(G);
        % [max_size, second_max_size] = getLargestComponentsSize(G);
        % Lcc_list(i)  = max_size;
        % sLcc_list(i) = second_max_size;
        nodei = randi(N);
        nodej = randi(N);
        while nodej == nodei
            nodej = randi(N);
        end
        [flowsubgraphlink,lsg] = flowsubgraph(G,nodei,nodej);
        
        flow_subgraph_linksize_list(i) = lsg;
        
        if lsg ~= 0
            FB = tril(A);
            FB(find(FB)) = flowsubgraphlink;
            FB = FB+FB.';
            flowsubgraphnode = sum(abs(FB)).';
            flowsubgraphnode = find(abs(flowsubgraphnode)>0.00000001);
            nodesize = size(flowsubgraphnode,1);
            flow_subgraph_size_list(i) = nodesize;
            
            % nodek = randi(N);
            % while nodek == nodei || nodek == nodej
            %     nodek = randi(N);
            % end
            % if ismember(nodek,flowsubgraphnode)
            %     nodep = nodep+1;
            % end
        end
%         p_nodebelongtoFSG = length(flowsubgraphnode)/N;
%         plist(i) = p_nodebelongtoFSG;
               
    end
%     p_ave(count) = mean(plist);  % ave probability of 
    % lcc_diffp(count) = mean(Lcc_list);
    % slcc_diffp(count) = mean(sLcc_list);
    % nodep_list(count) = nodep/siumtimes;
    avesize_fsg_list(count) = mean(flow_subgraph_size_list);
    stdsize_fsg_list(count) = std(flow_subgraph_size_list);
    
    filename_flow_subgraph_size = sprintf("%dnode\\size_fsg_p%.5f.txt",N,p);
    filename_flow_subgraph_size = fullfile(filefolder_name, filename_flow_subgraph_size);
    writematrix(flow_subgraph_size_list,filename_flow_subgraph_size)

    filename_real_ave_degree = sprintf("%dnode\\real_ave_degree_p%.5f.txt",N,p);
    filename_real_ave_degree = fullfile(filefolder_name, filename_real_ave_degree);
    writematrix(real_ave_degree_list,filename_real_ave_degree)

    filename_linkflow_subgraph_size = sprintf("%dnode\\linksize_fsg_p%.5f.txt",N,p);
    filename_linkflow_subgraph_size = fullfile(filefolder_name, filename_linkflow_subgraph_size);
    writematrix(flow_subgraph_linksize_list,filename_linkflow_subgraph_size)

    filename_real_linknum = sprintf("%dnode\\linknum_p%.5f.txt",N,p);
    filename_real_linknum = fullfile(filefolder_name, filename_real_linknum);
    writematrix(Linknumlist,filename_real_linknum)
        
%     AnalysisSolution(count) = SolutionAnalytic(N,p);

end
% matfilename = sprintf("%dnode\\prob_anodein_fsg_withdiff_p.txt",N);
% matfilename = fullfile(filefolder_name, matfilename);
% writematrix(nodep_list,matfilename);
% 
% avefsg_size_filename = sprintf("%dnode\\ave_fsg_withdiff_p.txt",N);
% avefsg_size_filename = fullfile(filefolder_name, avefsg_size_filename);
% writematrix(avesize_fsg_list,avefsg_size_filename);
% 
% stdfsg_size_filename = sprintf("%dnode\\std_fsg_withdiff_p.txt",N);
% stdfsg_size_filename = fullfile(filefolder_name, stdfsg_size_filename);
% writematrix(stdsize_fsg_list,stdfsg_size_filename);


% figure
figure()
plot(ave_degree,avesize_fsg_list/N)
hold on
% plot(ave_degree,lcc_diffp/N)
% hold on
% plot(ave_degree,slcc_diffp/N)
% hold on

 


function [max_size, second_max_size] = getLargestComponentsSize(G)
% 输入：
%   G - graph 或 digraph 对象
%
% 输出：
%   max_size - 最大联通子图的节点数
%   second_max_size - 第二大联通子图的节点数（如果只有一个，则为 0）

    % 获取每个节点所属的联通分量编号
    comp_ids = conncomp(G);

    % 统计每个联通分量的节点数量
    comp_sizes = histcounts(comp_ids, 1:(max(comp_ids)+1));

    % 如果只存在一个分量，返回第二个为0
    if isscalar(comp_sizes)
        max_size = comp_sizes(1);
        second_max_size = 0;
    else
        % 排序，找到最大和第二大
        sorted_sizes = sort(comp_sizes, 'descend');
        max_size = sorted_sizes(1);
        second_max_size = sorted_sizes(2);
    end
end
