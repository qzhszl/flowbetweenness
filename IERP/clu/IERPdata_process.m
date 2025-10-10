% 合并并删除 IERP 文件
clear; clc;

for N = [100,200] 
    x = log(N)/N;
    y = ceil(x * 1e4) / 1e4;  % round 4 decimal
    p_start = y;
        
    p_vec = linspace(p_start, 1, 15);
    % 前两个点
    p1 = p_vec(1);
    p2 = p_vec(2);  
    % 在 p1 和 p2 之间插入两个点
    extra_points = linspace(p1, p2, 4);  % 生成4个点
    extra_points = extra_points(2:3);    % 去掉第一个和最后一个（原本已有）
    % 合并
    p_vec = [p_vec(1), extra_points, p_vec(2:end)];
    p_vec = round(p_vec,4);

    for p =p_vec
        process_data(N,p)
    end

end

function process_data(N,p)
    folder = 'D:\data\flow betweenness\IERP\';
    
    % 存储合并数据
    mergedData = [];
    
    % 遍历 1~20 的 inputpara
    for inputpara = 1:20
        % 构造文件名
        filename = sprintf('%sIERP_N%dERp%.4f_weight_exp_simu%d.txt', folder, N, p, inputpara);
    
        % 检查文件是否存在
        if isfile(filename)
            % 读取数据
            data = readmatrix(filename);
    
            % 合并
            mergedData = [mergedData; data];
    
            % 删除原文件
            delete(filename);
            fprintf('已合并并删除: %s\n', filename);
        else
            fprintf('未找到文件: %s\n', filename);
        end
    end
    
    % 保存合并结果
    outputFile = sprintf('%sIERP_N%dERp%.4f_weight_exp.txt', folder, N, p);
    writematrix(mergedData, outputFile, 'Delimiter', ',');
    
    % fprintf('所有文件已合并并保存至: %s\n', outputFile);
end


