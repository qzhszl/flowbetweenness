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



% link number size
avgfsg_size_node_mat=[];
for inputpara =1:10
    avefsg_size_filename = sprintf("D:\\data\\flow betweenness\\sizeofflowsubgraph\\nodesize10000node\\%dnode\\linknum_N%d_simu%d.txt",N,N,inputpara);
    avgfsg_size_node_vec = readmatrix(avefsg_size_filename);
    avgfsg_size_node_mat = [avgfsg_size_node_mat, avgfsg_size_node_vec];
end
avgfsg_size_node_res = mean(avgfsg_size_node_mat,2);
ave_degree = 1:0.1:4.9;
ave_degree_2 = 5:10;
ave_degree = [ave_degree,ave_degree_2].';

avefsg_size_filename = sprintf("D:\\data\\flow betweenness\\sizeofflowsubgraph\\nodesize10000node\\%dnode\\link_num_withdiff_avg.txt",N);
T = table(ave_degree, avgfsg_size_node_res, 'VariableNames', {'avg','Mean'});
% 保存为制表符分隔的 txt 文件
writetable(T,avefsg_size_filename , 'Delimiter', '\t');


