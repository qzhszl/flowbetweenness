function SizeofFlowSubgraph_server(inputpara)
    rng(inputpara*10)
    N = 100000;
%     pc= log(N)/N;
    ave_degree = 1:0.1:4.9;
    ave_degree_2 = 5:10;
    ave_degree = [ave_degree,ave_degree_2];
    p_list = ave_degree/(N-1);
    
    siumtimes = 10;
    nodep_list =zeros(length(p_list),1);
    avesize_fsg_list =zeros(length(p_list),1);
    stdsize_fsg_list =zeros(length(p_list),1);
    % p_ave = zeros(length(p_list),1);
    % AnalysisSolution = zeros(length(p_list),1);
    % LCC_list = zeros(length(p_list),1);


%         plist = zeros(siumtimes,1);
    Linknum_mat = zeros(length(p_list),siumtimes);
%         Lcc = zeros(siumtimes,1);
    flow_subgraph_nodesize_mat = zeros(length(p_list),siumtimes);
    flow_subgraph_linksize_mat = zeros(length(p_list),siumtimes);
    real_ave_degree_mat = zeros(length(p_list),siumtimes);

    count=0;
    
    
    for p=p_list(1:length(p_list))
        disp(p)
        count = count+1;
        nodep=0;
        %         plist = zeros(siumtimes,1);
        Linknumlist = zeros(1,siumtimes);
    %         Lcc = zeros(siumtimes,1);
        flow_subgraph_nodesize_list = zeros(1,siumtimes);
        flow_subgraph_linksize_list = zeros(1,siumtimes);
        real_ave_degree_list = zeros(1,siumtimes);


        for i = 1:siumtimes
            disp(i)
            A = rand(N,N) < p;
            A = triu(A,1);
            A = A + A';
            % writematrix(A,"D:\\data\\flow betweenness\\analysisbranchingprocess\\100nodeER")
            % A = readmatrix("D:\\data\\flow betweenness\\analysisbranchingprocess\\100nodeER");
            
            real_ave_degree = mean(sum(A));
            real_ave_degree_list(i) = real_ave_degree;
            
            G = graph(A);
            Linknumlist(i) = numedges(G);
        %     plot(G)
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
                flow_subgraph_nodesize_list(i) = nodesize-2;
                nodek = randi(N);
                while nodek == nodei || nodek == nodej
                    nodek = randi(N);
                end
                if ismember(nodek,flowsubgraphnode)
                    nodep = nodep+1;
                end
            end             
        end
        nodep_list(count) = nodep/siumtimes;
        avesize_fsg_list(count) = mean(flow_subgraph_nodesize_list);
        Linknum_mat(count,:) = Linknumlist;
        flow_subgraph_nodesize_mat(count,:) = flow_subgraph_nodesize_list;
        flow_subgraph_linksize_mat(count,:) = flow_subgraph_linksize_list;
        real_ave_degree_mat(count,:) = real_ave_degree_list;

    end
    filename_flow_subgraph_linksize = sprintf("linksize_fsg_N%d_simu%d.txt",N,inputpara);
    writematrix(flow_subgraph_linksize_mat,filename_flow_subgraph_linksize)        
    filename_flow_subgraph_nodesize = sprintf("nodesize_fsg_N%d_simu%d.txt",N,inputpara);
    writematrix(flow_subgraph_nodesize_mat,filename_flow_subgraph_nodesize)
    filename_real_ave_degree = sprintf("real_ave_degree_N%d_simu%d.txt",N,inputpara);
    writematrix(real_ave_degree_mat,filename_real_ave_degree)
    filename_linknum = sprintf("linknum_N%d_simu%d.txt",N,inputpara);
    writematrix(flow_subgraph_linksize_mat,filename_linknum)


    matfilename = sprintf("prob_anodein_fsg_withdiff_p_simu%d.txt",N,inputpara);
    writematrix(nodep_list,matfilename);
    avefsg_size_filename = sprintf("ave_fsg_withdiff_p_simu%d.txt",N,inputpara);
    writematrix(avesize_fsg_list,avefsg_size_filename);
%     stdfsg_size_filename = sprintf("std_fsg_withdiff_p_%d.txt",N,inputpara);
%     writematrix(stdsize_fsg_list,stdfsg_size_filename);
end

 

