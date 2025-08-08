function SizeofFlowSubgraph_server(inputpara)
    rng(inputpara*10)
    N = 10000;
    pc= log(N)/N;
    ave_degree = 1:0.1:4.9;
    ave_degree_2 = 5:10;
    ave_degree = [ave_degree,ave_degree_2];
    p_list = ave_degree/(N-1);
    
    
    nodep_list =zeros(length(p_list),1);
    avesize_fsg_list =zeros(length(p_list),1);
    stdsize_fsg_list =zeros(length(p_list),1);
    % p_ave = zeros(length(p_list),1);
    % AnalysisSolution = zeros(length(p_list),1);
    % LCC_list = zeros(length(p_list),1);
    count=0;
    siumtimes = 10;
    
    for p=p_list(1:length(p_list))
        disp(p)
        count = count+1;
        nodep=0;
        plist = zeros(siumtimes,1);
        Linknumlist = zeros(siumtimes,1);
        Lcc = zeros(siumtimes,1);
        flow_subgraph_size_list = zeros(siumtimes,1);
        real_ave_degree_list = zeros(siumtimes,1);
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
            
            if lsg ~= 0
                FB = tril(A);
                FB(find(FB)) = flowsubgraphlink;
                FB = FB+FB.';
                flowsubgraphnode = sum(abs(FB)).';
                flowsubgraphnode = find(abs(flowsubgraphnode)>0.00000001);
                nodesize = size(flowsubgraphnode,1);
                flow_subgraph_size_list(i) = nodesize;
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
        avesize_fsg_list(count) = mean(flow_subgraph_size_list);
        stdsize_fsg_list(count) = std(flow_subgraph_size_list);
        filename_flow_subgraph_size = sprintf("D:\\data\\flow betweenness\\sizeofflowsubgraph\\nodesize10000node\\%dnode\\size_fsg_p%.5f_%d.txt",N,p,inputpara);
        writematrix(flow_subgraph_size_list,filename_flow_subgraph_size)
        filename_real_ave_degree = sprintf("D:\\data\\flow betweenness\\sizeofflowsubgraph\\nodesize10000node\\%dnode\\real_ave_degree_p%.5f_%d.txt",N,p,inputpara);
        writematrix(real_ave_degree_list,filename_real_ave_degree)
    end
    matfilename = sprintf("D:\\data\\flow betweenness\\sizeofflowsubgraph\\nodesize10000node\\%dnode\\prob_anodein_fsg_withdiff_p_%d.txt",N,inputpara);
    writematrix(nodep_list,matfilename);
    avefsg_size_filename = sprintf("D:\\data\\flow betweenness\\sizeofflowsubgraph\\nodesize10000node\\%dnode\\ave_fsg_withdiff_p_%d.txt",N,inputpara);
    writematrix(avesize_fsg_list,avefsg_size_filename);
    stdfsg_size_filename = sprintf("D:\\data\\flow betweenness\\sizeofflowsubgraph\\nodesize10000node\\%dnode\\std_fsg_withdiff_p_%d.txt",N,inputpara);
    writematrix(stdsize_fsg_list,stdfsg_size_filename);
end

 

