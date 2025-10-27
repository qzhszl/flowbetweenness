clear,clc
% colors = ["#D08082", "#C89FBF", "#62ABC7", "#7A7DB1", "#6FB494", "#D9B382"];
% count = 1;
% for N = [100,1000,10000]
%     avg = 0:0.1:10;
%     p_vals = avg/(N-1);
%     s_vals = zeros(size(p_vals));
%     for i = 1:length(p_vals)
%         s_vals(i) = compute_s_from_ER(N, p_vals(i));
%     end
%     % s_vals = s_vals.^1.5
%     % plot(avg, s_vals, 'LineWidth', 4, Color=colors(count))
% 
%     filename = sprintf('analytic_N%d.txt', N);
% 
%     % 打开文件写入（覆盖写）
%     fid = fopen(filename, 'w');
% 
%     % 写入 avg（第一行，用逗号分隔）
%     fprintf(fid, '%.4f', avg(1));
%     fprintf(fid, ',%.4f', avg(2:end));
%     fprintf(fid, '\n');
% 
%     % 写入 s_vals（第二行，用逗号分隔）
%     fprintf(fid, '%.4f', s_vals(1));
%     fprintf(fid, ',%.4f', s_vals(2:end));
%     fprintf(fid, '\n');
% 
%     % 关闭文件
%     fclose(fid);
% 
%     count = count+1;
%     hold on
% end

% count = 1;
% for N = [10,20,30,50,80,100,100,1000,10000]
%     filefolder_name = "D:\\data\\flow betweenness\\sizeofflowsubgraph\\new";
%     outname = fullfile(filefolder_name, sprintf('%dnode_results_summary.csv', N));
%     result_table = readtable(outname);
%     % plot(result_table.RealAveDegree,result_table.SizeFSG/N,'o--','LineWidth', 2, 'MarkerSize', 10,Color=colors(5))
%     avg = result_table.RealAveDegree;
%     s_vals = result_table.SizeFSG/N;
% 
%     filename = sprintf('simu_N%d.txt', N);
% 
%     % 打开文件写入（覆盖写）
%     fid = fopen(filename, 'w');
% 
%     % 写入 avg（第一行，用逗号分隔）
%     fprintf(fid, '%.4f', avg(1));
%     fprintf(fid, ',%.4f', avg(2:end));
%     fprintf(fid, '\n');
% 
%     % 写入 s_vals（第二行，用逗号分隔）
%     fprintf(fid, '%.4f', s_vals(1));
%     fprintf(fid, ',%.4f', s_vals(2:end));
%     fprintf(fid, '\n');
% 
%     % 关闭文件
%     fclose(fid);
% 
%     count = count+1;
%     hold on
% end



% for links:
% count = 1;
% for N = [10,20,30,50,80,100,100,1000,10000]
%     filefolder_name = "D:\\data\\flow betweenness\\sizeofflowsubgraph\\new";
%     outname = fullfile(filefolder_name, sprintf('%dnode_results_summary.csv', N));
%     result_table = readtable(outname);
%     plot(result_table.RealAveDegree,result_table.LinkSizeFSG./result_table.LinkNum,'o--','LineWidth', 2, 'MarkerSize', 10)
%     hold on
%     avg = result_table.RealAveDegree;
%     s_vals = result_table.LinkSizeFSG./result_table.LinkNum;
% 
% 
%     filename = sprintf('simu_link_N%d.txt', N);
% 
%     % 打开文件写入（覆盖写）
%     fid = fopen(filename, 'w');
% 
%     % 写入 avg（第一行，用逗号分隔）
%     fprintf(fid, '%.4f', avg(1));
%     fprintf(fid, ',%.4f', avg(2:end));
%     fprintf(fid, '\n');
% 
%     % 写入 s_vals（第二行，用逗号分隔）
%     fprintf(fid, '%.4f', s_vals(1));
%     fprintf(fid, ',%.4f', s_vals(2:end));
%     fprintf(fid, '\n');
% 
%     % 关闭文件
%     fclose(fid);
% 
%     count = count+1;
%     hold on
% end

count = 1;
for N = [10,20,30,50,80,100,100,1000,10000]
    avg = 0:0.1:10;
    p_vals = avg/(N-1);
    s_vals = zeros(size(p_vals));
    for i = 1:length(p_vals)
        s_vals(i) = compute_S_link_from_ER(N, p_vals(i));
    end
    % s_vals = s_vals.^1.5
    % plot(avg, s_vals, 'LineWidth', 4, Color=colors(count))

    filename = sprintf('analytic_link_N%d.txt', N);

    % 打开文件写入（覆盖写）
    fid = fopen(filename, 'w');

    % 写入 avg（第一行，用逗号分隔）
    fprintf(fid, '%.4f', avg(1));
    fprintf(fid, ',%.4f', avg(2:end));
    fprintf(fid, '\n');

    % 写入 s_vals（第二行，用逗号分隔）
    fprintf(fid, '%.4f', s_vals(1));
    fprintf(fid, ',%.4f', s_vals(2:end));
    fprintf(fid, '\n');

    % 关闭文件
    fclose(fid);

    count = count+1;
    hold on
end


