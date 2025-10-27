clear,clc
% plot the relation between the flow subgraph node size with the real ave degree


fig = figure; 
fig.Position = [100 100 900 600]; 
hold on;
colors = ["#D08082", "#C89FBF", "#62ABC7", "#7A7DB1", "#6FB494", "#D9B382"];
% colors2 = ["#D08082","#6FB494","#D9B382","#C89FBF",];
count = 1;
for N = [10,50,100,10000]
    avg = 0:0.1:10;
    p_vals = avg/(N-1);
    s_vals = zeros(size(p_vals));
    for i = 1:length(p_vals)
        s_vals(i) = compute_S_link_from_ER(N, p_vals(i));
    end
    % s_vals = s_vals.^1.5
    plot(avg, s_vals, 'LineWidth', 4, Color=colors(count))
    count = count+1;
    hold on
end


count = 1;
% for N = [10,20,30,50,80,100,100,1000,10000]
for N = [10,50,100,10000]
    filefolder_name = "D:\\data\\flow betweenness\\sizeofflowsubgraph\\new";
    outname = fullfile(filefolder_name, sprintf('%dnode_results_summary.csv', N));
    result_table = readtable(outname);
    plot(result_table.RealAveDegree,result_table.LinkSizeFSG./result_table.LinkNum,'LineStyle', 'none', 'Marker', '*','LineWidth', 2, 'MarkerSize', 10, Color=colors(count))
    count = count+1;
    hold on
end

legend("10","50","100","10000")
