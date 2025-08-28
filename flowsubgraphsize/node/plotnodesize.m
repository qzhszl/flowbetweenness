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


x = [1	1.10000000000000	1.20000000000000	1.30000000000000	1.40000000000000	1.50000000000000	1.60000000000000	1.70000000000000	1.80000000000000	1.90000000000000	2	2.10000000000000	2.20000000000000	2.30000000000000	2.40000000000000	2.50000000000000	2.60000000000000	2.70000000000000	2.80000000000000	2.90000000000000	3	3.10000000000000	3.20000000000000	3.30000000000000	3.40000000000000	3.50000000000000	3.60000000000000	3.70000000000000	3.80000000000000	3.90000000000000	4	4.10000000000000	4.20000000000000	4.30000000000000	4.40000000000000	4.50000000000000	4.60000000000000	4.70000000000000	4.80000000000000	4.90000000000000	5	6	7	8	9	10
]
y = [0.5012
0.8618
1.7104
3.1167
4.9296
7.8625
11.2092
15.8194
20.4515
25.1236
30.4038
34.8763
40.0192
44.5795
49.1496
53.0437
57.1763
59.841
63.3747
66.2774
69.2454
71.9736
74.1801
76.3822
78.1056
80.7693
81.763
83.364
85.0065
86.1047
87.1364
88.6903
89.8674
90.4199
91.2431
92.0713
92.8297
93.3789
94.0907
94.528
95.0764
98.0062
99.174
99.7519
99.9073
99.9685
]

plot(x,y/100)
hold on
lgd=legend("10000","ana","1000","100","100_2")


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
