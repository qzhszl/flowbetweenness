clear,clc

N = 100;
ave_degree = 0.2:0.2:4.9;
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
        switch file_types{j}
            case 'size_fsg'
                filename = sprintf("%dnode\\size_fsg_p%.5f.txt", N, p);
            case 'real_ave_degree'
                filename = sprintf("%dnode\\real_ave_degree_p%.5f.txt", N, p);
            case 'linksize_fsg'
                filename = sprintf("%dnode\\linksize_fsg_p%.5f.txt", N, p);
            case 'linknum'
                filename = sprintf("%dnode\\linknum_p%.5f.txt", N, p);
        end
        filepath = fullfile(filefolder_name, filename);

        % 读入数据并取平均
        data = readmatrix(filepath);
        row_data(j) = mean(data(:));
    end

    % 添加到结果表
    result = [result; array2table([p, row_data], ...
                'VariableNames', ['p', col_names])];
end


% 保存结果表
outname = fullfile(filefolder_name, sprintf('%dnode_results_summary.csv', N));
writetable(result, outname);



