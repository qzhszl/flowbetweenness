clear,clc
% plot the relation between the flow subgraph node size with the real ave degree
N = 10000

avefsg_size_filename = sprintf("D:\\data\\flow betweenness\\sizeofflowsubgraph\\nodesize10000node\\%dnode\\ave_fsg_withdiff_avg.txt",N);
% avefsg_size_filename = sprintf("D:\\data\\flow betweenness\\sizeofflowsubgraph\\nodesize10000node\\%dnode\\ave_fsg_withdiff_avg.txt",N);
result_table = readtable(avefsg_size_filename);
plot(result_table.avg,result_table.Mean/N)
hold on

% p_vals = result_table.avg/(N-1);
% s_vals = zeros(size(p_vals));
avg = 0:0.1:8;
p_vals = avg/(N-1);
s_vals = zeros(size(p_vals));
for i = 1:length(p_vals)
    s_vals(i) = compute_s_from_ER(N, p_vals(i));
end
plot(avg, s_vals, 'LineWidth', 2)


N = 1000

avefsg_size_filename = sprintf("D:\\data\\flow betweenness\\sizeofflowsubgraph\\nodesize10000node\\%dnode\\ave_fsg_withdiff_avg.txt",N);
% avefsg_size_filename = sprintf("D:\\data\\flow betweenness\\sizeofflowsubgraph\\nodesize10000node\\%dnode\\ave_fsg_withdiff_avg.txt",N);
result_table = readtable(avefsg_size_filename);
plot(result_table.avg,result_table.Mean/N)
hold on

N = 100

avefsg_size_filename = sprintf("D:\\data\\flow betweenness\\sizeofflowsubgraph\\nodesize10000node\\%dnode\\ave_fsg_withdiff_avg.txt",N);
% avefsg_size_filename = sprintf("D:\\data\\flow betweenness\\sizeofflowsubgraph\\nodesize10000node\\%dnode\\ave_fsg_withdiff_avg.txt",N);
result_table = readtable(avefsg_size_filename);
plot(result_table.avg,result_table.Mean/N)
hold on

lgd=legend("10000","ana","1000","100")


% for large network
%__________________________________________________________________________
% N = 10000;
% pc= log(N)/N;
% ave_degree = 1:10;
% p_list = ave_degree/(N-1);
% size_fsubgraph = zeros(10,1);
% nodesize_data =[];
% for i = 1:10
%     filename = sprintf("D:\\data\\flow betweenness\\sizeofflowsubgraph\\nodesize10000node\\task%d.txt",i);
%     nodesizep = readmatrix(filename);
%     nodesize_data = [nodesize_data,nodesizep] 
% end
% ave = mean(nodesize_data,2)
% std = std(nodesize_data,0,2)
% plot(ave_degree,ave)
