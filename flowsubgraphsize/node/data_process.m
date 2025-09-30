clear,clc

N = 10000;

% node subgraph size
% avgfsg_size_node_mat=[];
% for inputpara =1:10
%     avefsg_size_filename = sprintf("D:\\data\\flow betweenness\\sizeofflowsubgraph\\nodesize10000node\\%dnode\\nodesize_fsg_N%d_simu%d.txt",N,N,inputpara);
%     avgfsg_size_node_vec = readmatrix(avefsg_size_filename);
%     avgfsg_size_node_mat = [avgfsg_size_node_mat, avgfsg_size_node_vec];
% end
% avgfsg_size_node_res = mean(avgfsg_size_node_mat,2);
% ave_degree = 1:0.1:4.9;
% ave_degree_2 = 5:10;
% ave_degree = [ave_degree,ave_degree_2].';
% 
% avefsg_size_filename = sprintf("D:\\data\\flow betweenness\\sizeofflowsubgraph\\nodesize10000node\\%dnode\\ave_fsg_withdiff_avg.txt",N);
% T = table(ave_degree, avgfsg_size_node_res, 'VariableNames', {'avg','Mean'});
% % 保存为制表符分隔的 txt 文件
% writetable(T,avefsg_size_filename , 'Delimiter', '\t');


% link subgraph size
% avgfsg_size_node_mat=[];
% for inputpara =1:10
%     avefsg_size_filename = sprintf("D:\\data\\flow betweenness\\sizeofflowsubgraph\\nodesize10000node\\%dnode\\linksize_fsg_N%d_simu%d.txt",N,N,inputpara);
%     avgfsg_size_node_vec = readmatrix(avefsg_size_filename);
%     avgfsg_size_node_mat = [avgfsg_size_node_mat, avgfsg_size_node_vec];
% end
% avgfsg_size_node_res = mean(avgfsg_size_node_mat,2);
% ave_degree = 1:0.1:4.9;
% ave_degree_2 = 5:10;
% ave_degree = [ave_degree,ave_degree_2].';
% 
% avefsg_size_filename = sprintf("D:\\data\\flow betweenness\\sizeofflowsubgraph\\nodesize10000node\\%dnode\\ave_fsglink_withdiff_avg.txt",N);
% T = table(ave_degree, avgfsg_size_node_res, 'VariableNames', {'avg','Mean'});
% % 保存为制表符分隔的 txt 文件
% writetable(T,avefsg_size_filename , 'Delimiter', '\t');


% Previous data
% link number size
% avgfsg_size_node_mat=[];
% for inputpara =1:10
%     avefsg_size_filename = sprintf("D:\\data\\flow betweenness\\sizeofflowsubgraph\\nodesize10000node\\%dnode\\linknum_N%d_simu%d.txt",N,N,inputpara);
%     avgfsg_size_node_vec = readmatrix(avefsg_size_filename);
%     avgfsg_size_node_mat = [avgfsg_size_node_mat, avgfsg_size_node_vec];
% end
% avgfsg_size_node_res = mean(avgfsg_size_node_mat,2);
% ave_degree = 1:0.1:4.9;
% ave_degree_2 = 5:10;
% ave_degree = [ave_degree,ave_degree_2].';
% 
% avefsg_size_filename = sprintf("D:\\data\\flow betweenness\\sizeofflowsubgraph\\nodesize10000node\\%dnode\\link_num_withdiff_avg.txt",N);
% T = table(ave_degree, avgfsg_size_node_res, 'VariableNames', {'avg','Mean'});
% % 保存为制表符分隔的 txt 文件
% writetable(T,avefsg_size_filename , 'Delimiter', '\t');


N = 10000;
ave_degree = 1:0.2:4.9;
ave_degree_2 = 5:10;
ave_degree = [ave_degree,ave_degree_2];
p_list = ave_degree/(N-1);

% 假设已有 N, p_list, filefolder_name
result = table;   % 初始化一个空表

% 定义文件类别和列名
file_types = { ...
    'size_fsg', ...
    'real_ave_degree', ...
    'linksize_fsg', ...
    'linknum'};

col_names = {'SizeFSG', 'RealAveDegree', 'LinkSizeFSG', 'LinkNum'};

filefolder_name = "D:\\data\\flow betweenness\\sizeofflowsubgraph\\new";

% 遍历 p_list
for i = 1:length(p_list)
    p = p_list(i);

    % 当前行的数据
    row_data = zeros(1, length(file_types));

    % 遍历四类文件
    for j = 1:length(file_types)
        data =[];
        for inputpara =1:100
            switch file_types{j}
                case 'size_fsg'
                    filename = sprintf("%dnode\\size_fsg_p%.5f_%d.txt", N, p,inputpara);
                case 'real_ave_degree'
                    filename = sprintf("%dnode\\real_ave_degree_p%.5f_%d.txt", N, p,inputpara);
                case 'linksize_fsg'
                    filename = sprintf("%dnode\\linksize_fsg_p%.5f_%d.txt", N, p,inputpara);
                case 'linknum'
                    filename = sprintf("%dnode\\linknum_p%.5f_%d.txt", N, p,inputpara);
            end
            filepath = fullfile(filefolder_name, filename);
            % concat all the data
            try
                data_onefile = readmatrix(filepath);
                data = [data; data_onefile];
            catch
                % 读取失败就跳过
                fprintf('文件无法读取，已跳过: %s\n', filepath);
                continue;
            end
        end   
        row_data(j) = mean(data(:));
    end

    % 添加到结果表
    result = [result; array2table([p, row_data], ...
                'VariableNames', ['p', col_names])];
end


% 保存结果表
outname = fullfile(filefolder_name, sprintf('%dnode_results_summary.csv', N));
writetable(result, outname);
