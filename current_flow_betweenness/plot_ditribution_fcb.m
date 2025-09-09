clear,clc
filefolder_name = "D:\\data\\flow betweenness\\";

N = 50;
p = 0.2;
resname  = sprintf('bet_cbet_degree_N%dp%.2fER.mat',N,p);
filename = filefolder_name+resname;

% load_and_plot_nodeCFB.m
S = load(filename);
results = S.results;

% 2. 拼接所有实验的 nodeCFB
allNodeCFB = [];
for r = 1:numel(results)
    allNodeCFB = [allNodeCFB; results(r).nodeCFB(:)];
end


% alllinkCFB = [];
% for r = 1:numel(results)
%     alllinkCFB = [alllinkCFB; results(r).edgeCFB(:)];
% end

allnodeSPB = [];
for r = 1:numel(results)
    allnodeSPB = [allnodeSPB; results(r).nodeSPB(:)];
end

allnodedegree = [];
for r = 1:numel(results)
    allnodedegree = [allnodedegree; results(r).deg(:)];
end


% 3. 画直方图 (分布图)
figure;
histogram(allNodeCFB, 50, 'Normalization', 'pdf'); % 50 bins
set(gca,"YScale", "log")
hold on

histogram(allnodeSPB, 50, 'Normalization', 'pdf'); % 50 bins

ylabel('$f_b(x)$','interpreter','latex','FontSize',30)
xlabel('$x$','interpreter','latex','FontSize',30);




% 4. 如果想要核密度估计 (平滑曲线)
% figure;
% ksdensity(allNodeCFB);
% xlabel('Node Current-Flow Betweenness');
% ylabel('Density');
% title('KDE of Node Current-Flow Betweenness');


% 4. 平均度和flow-current betweenness的关系
figure;
scatter(allnodedegree,allNodeCFB)

uniqueDegrees = unique(allnodedegree);
meanCFB = zeros(size(uniqueDegrees));
stdCFB = zeros(size(uniqueDegrees));
for i = 1:length(uniqueDegrees)
    deg = uniqueDegrees(i);
    meanCFB(i) = mean(allNodeCFB(allnodedegree == deg));
    stdCFB(i) = std(allNodeCFB(allnodedegree == deg));
end

figure;
errorbar(uniqueDegrees,meanCFB,stdCFB)
% set(gca,"YScale", "log")
% set(gca,"XScale", "log")
xlabel('Node Degree');
ylabel('Mean CFB');