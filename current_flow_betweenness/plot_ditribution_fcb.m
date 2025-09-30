clear,clc
filefolder_name = "D:\\data\\flow betweenness\\";

N = 1000;
avg = 50
p = avg/(N-1)
if N == 1000
    resname  = sprintf('bet_cbet_degree_N%dp%.4fER_unweighted.mat',N,p);
else
    resname  = sprintf('bet_cbet_degree_N%dp%.2fER_unweighted.mat',N,p);
end
% resname = "bet_cbet_degree_BAnetworkN1000m3.mat"
pos = strfind(resname,"_");
pos  =pos(end);
pos2 = strfind(resname,".");
pos2 = pos2(end);
file_network_name = resname(pos:pos2-1);

filename = filefolder_name+resname;

% load_and_plot_nodeCFB.m
S = load(filename);
results = S.results;

% 2. 拼接所有实验的 nodeCFB
allNodeCFB = [];
for r = 1:numel(results)
    allNodeCFB = [allNodeCFB; results(r).nodeCFB(:)./nchoosek(N,2)];
end


% alllinkCFB = [];
% for r = 1:numel(results)
%     alllinkCFB = [alllinkCFB; results(r).edgeCFB(:)];
% end

allnodeSPB = [];
for r = 1:numel(results)
    allnodeSPB = [allnodeSPB; results(r).nodeSPB(:)./nchoosek(N,2)];
end

allnodedegree = [];
for r = 1:numel(results)
    allnodedegree = [allnodedegree; results(r).deg(:)];
end


% % 3. plot distribution画直方图 (分布图)
figure;
% allnodedegree_normalized = allnodedegree./N;
% h0 = histogram(allnodedegree_normalized, 50, 'Normalization', 'pdf');
% hold on

% colors = ["#D08082", "#C89FBF", "#62ABC7", "#7A7DB1", "#6FB494", "#D9B382"];

colors = ["#1F77B4","#FF7F0E", "#2CA02C","#D62728"]

h1 = histogram(allNodeCFB, 50, 'Normalization', 'pdf','FaceColor', colors(1)); % 50 bins
hold on
binCenters1 = h1.BinEdges(1:end-1) + diff(h1.BinEdges)/2;
counts1 = h1.Values / trapz(binCenters1,h1.Values); % 归一化

ft = fittype('a*exp(-((x-b)^2)/(2*c^2))', ...
             'independent','x','coefficients',{'a','b','c'});
% 初值猜测
startA = max(counts1);
startB = sum(binCenters1 .* counts1) / sum(counts1); % 加权平均估计均值
startC = std(binCenters1); % 粗略估计

% 拟合
[fitresult, gof] = fit(binCenters1.', counts1.', ft, ...
                       'StartPoint', [startA, startB, startC]);

% 提取参数
a = fitresult.a;
b = fitresult.b;
c = fitresult.c;
x_fit_normal = linspace(min(binCenters1),max(binCenters1),100);
y_fit_normal = a*exp(-((x_fit_normal-b).^2)/(2*c^2));
normStr = sprintf('$y = %.4f e^{-\\frac{(x - %.4f)^2}{%.5f}}$', a, b,2*c^2);


% h2 = histogram(allnodeSPB, 50, 'Normalization', 'pdf','FaceColor', colors(2)); % 50 bins
% hold on
% binCenters = h2.BinEdges(1:end-1) + diff(h2.BinEdges)/2;
% counts = h2.Values / trapz(binCenters,h2.Values); % 归一化
% counts = log(counts);
% 
% idx = find(~isinf(counts));
% binCenters = binCenters(idx);
% counts = counts(idx);
% 
% 
% f = polyfit(binCenters,counts,1);
% b = f(1);
% loga = f(2);
% a = exp(loga);
% x_fit = linspace(min(binCenters),max(binCenters),100);
% y_fit = a*exp(b*x_fit);
% 
% % 绘制拟合曲线
% plot(x_fit,y_fit,'r-','LineWidth',2,'Color', colors(3));
% hold on

plot(x_fit_normal, y_fit_normal, 'LineWidth', 2,'Color', colors(4))         
hold on

ylabel('$f_b(x)$','interpreter','latex','FontSize',30)
xlabel('$x$','interpreter','latex','FontSize',30);
% ylim([0.0001,10000])

eqnStr = sprintf('$y = %.4f e^{%.4f x}$', a, b);

% lgd = legend('current-flow betweenness','path betweenness',eqnStr,normStr, 'Interpreter','latex','FontSize',18,'Location','northeast');
lgd = legend('current-flow betweenness',normStr, 'Interpreter','latex','FontSize',18,'Location','northeast');
pos = lgd.Position;    % [x y width height]
pos(1) = pos(1) + 0.014; % 向右挪一点
pos(2) = pos(2) + 0.02; % 向上挪一点
lgd.Position = pos;
set(gca, 'FontSize', 22, ...   % 坐标轴字体大小
         'TickDir', 'in', ...
         'TickLength', [0.02 0.02]);     % 刻度线朝里
% set(gca,"XScale", "log")
% set(gca,"YScale", "log")
figure_name = filefolder_name+"distribution_fcb_pb"+file_network_name+".pdf";
% print(gcf,figure_name, '-dpdf', '-r600');  % 保存为 scatter_plot.pdf




% 4. 如果想要核密度估计 (平滑曲线)
% figure;
% ksdensity(allNodeCFB);
% xlabel('Node Current-Flow Betweenness');
% ylabel('Density');
% title('KDE of Node Current-Flow Betweenness');


% 5. 平均度和flow-current betweenness的关系
figure;
scatter(allnodedegree,allNodeCFB,'Color',"#FF7F0E")
R = corrcoef(allnodedegree, allNodeCFB);
r_value = R(1,2);

x_pos = min(allnodedegree) + 0.05*(max(allnodedegree)-min(allnodedegree));
y_pos = max(allNodeCFB) - 0.02*(max(allNodeCFB)-min(allNodeCFB));
text(x_pos, y_pos, sprintf('r = %.2f', r_value), ...
     'FontSize', 24);

xlabel('Degree',Interpreter='latex',FontSize=24);
ylabel('Current-flow betweenness',Interpreter='latex',FontSize=24);
set(gca, 'FontSize', 22, ...   % 坐标轴字体大小
         'TickDir', 'in', ...
         'TickLength', [0.02 0.02]);     % 刻度线朝里
box on

figure_name = filefolder_name+"scatter_degree_cfb"+file_network_name+".pdf";
print(gcf,figure_name, '-dpdf', '-r600');  % 保存为 scatter_plot.pdf



uniqueDegrees = unique(allnodedegree);
meanCFB = zeros(size(uniqueDegrees));
stdCFB = zeros(size(uniqueDegrees));
for i = 1:length(uniqueDegrees)
    deg = uniqueDegrees(i);
    meanCFB(i) = mean(allNodeCFB(allnodedegree == deg));
    stdCFB(i) = std(allNodeCFB(allnodedegree == deg));
end

figure;
h = errorbar(uniqueDegrees, meanCFB, stdCFB, '-o', 'LineWidth', 3, ...    % 主线加粗
    'MarkerSize', 10, ...                            % 点大小
    'CapSize', 12, ...
    'Color',"#1F77B4");                          
hold on
p = polyfit(uniqueDegrees, meanCFB, 1);   % p(1)=斜率, p(2)=截距
y_fit = polyval(p, uniqueDegrees);
plot(uniqueDegrees, y_fit, 'r-', 'LineWidth', 2);
sprintf('Linear Fit: y = %.6fx + %.6f', p(1), p(2))
target_a = 1/avg^(1.5)/sqrt(pi)
target_a/p(1)
% set(gca,"YScale", "log")
% set(gca,"XScale", "log")
% ylim([0,1200])
xlabel('Degree',Interpreter='latex',FontSize=24);
ylabel('Current-flow betweenness',Interpreter='latex',FontSize=24);
set(gca, 'FontSize', 22, ...   % 坐标轴字体大小
         'TickDir', 'in', ...
         'TickLength', [0.02 0.02]);     % 刻度线朝里

figure_name = filefolder_name+"degree_avecfb"+file_network_name+".pdf";
print(gcf,figure_name, '-dpdf', '-r600');  % 保存为 scatter_plot.pdf




% % 6. 平均度和path betweenness的关系
% figure;
% scatter(allnodedegree,allnodeSPB,'MarkerEdgeColor',"#FF7F0E")
% R = corrcoef(allnodedegree, allnodeSPB);
% r_value = R(1,2);
% 
% x_pos = min(allnodedegree) + 0.05*(max(allnodedegree)-min(allnodedegree));
% y_pos = max(allnodeSPB) - 0.02*(max(allnodeSPB)-min(allnodeSPB));
% text(x_pos, y_pos, sprintf('r = %.2f', r_value), ...
%      'FontSize', 24);
% 
% xlabel('Degree',Interpreter='latex',FontSize=24);
% ylabel('Path betweenness',Interpreter='latex',FontSize=24);
% set(gca, 'FontSize', 22, ...   % 坐标轴字体大小
%          'TickDir', 'in', ...
%          'TickLength', [0.02 0.02]);     % 刻度线朝里
% box on
% 
% figure_name = filefolder_name+"scatter_degree_spb"+file_network_name+".pdf";
% print(gcf,figure_name, '-dpdf', '-r600');  % 保存为 scatter_plot.pdf
% 
% 
% 
% uniqueDegrees = unique(allnodedegree);
% meanCFB = zeros(size(uniqueDegrees));
% stdCFB = zeros(size(uniqueDegrees));
% for i = 1:length(uniqueDegrees)
%     deg = uniqueDegrees(i);
%     meanCFB(i) = mean(allnodeSPB(allnodedegree == deg));
%     stdCFB(i) = std(allnodeSPB(allnodedegree == deg));
% end
% 
% figure;
% h = errorbar(uniqueDegrees, meanCFB, stdCFB, '-o', 'LineWidth', 3, ...    % 主线加粗
%     'MarkerSize', 10, ...                            % 点大小
%     'CapSize', 12, ...
%     'Color',"#FF7F0E");                          
% 
% 
% % set(gca,"YScale", "log")
% % set(gca,"XScale", "log")
% % ylim([0,1200])
% xlabel('Degree',Interpreter='latex',FontSize=24);
% ylabel('Path betweenness',Interpreter='latex',FontSize=24);
% set(gca, 'FontSize', 22, ...   % 坐标轴字体大小
%          'TickDir', 'in', ...
%          'TickLength', [0.02 0.02]);     % 刻度线朝里
% 
% figure_name = filefolder_name+"degree_avespb"+file_network_name+".pdf";
% print(gcf,figure_name, '-dpdf', '-r600');  % 保存为 scatter_plot.pdf


