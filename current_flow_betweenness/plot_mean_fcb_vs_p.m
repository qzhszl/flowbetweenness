clear,clc
filefolder_name = "D:\\data\\flow betweenness\\";

N = 1000;
flow_current_node_bet_mean =[];
flow_current_node_bet_std =[];
if N == 1000
    avg = [8,10,20,50,100,200,500];
    p_vec = avg./(N-1);
elseif N ==200  
    p_vec = [0.03,0.05,0.1,0.2,0.3,0.4];
elseif N ==100
    p_vec = [0.05,0.09,0.1,0.2,0.3,0.4];
elseif N == 50
    p_vec = [0.07,0.08,0.1,0.2,0.3,0.4];
end

count = 0;
for p =p_vec
    count = count+1;
    if N>999
        resname  = sprintf('bet_cbet_degree_N%dp%.4fER_unweighted.mat',N,p);
    else
        resname  = sprintf('bet_cbet_degree_N%dp%.2fER.mat',N,p);
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

    flow_current_node_bet_mean(count) = mean(allNodeCFB);
    flow_current_node_bet_std(count) = std(allNodeCFB);
    % alllinkCFB = [];
    % for r = 1:numel(results)
    %     alllinkCFB = [alllinkCFB; results(r).edgeCFB(:)];
    % end

end

figure;
h = errorbar(p_vec, flow_current_node_bet_mean, flow_current_node_bet_std, '-o', 'LineWidth', 3, ...    % 主线加粗
    'MarkerSize', 10, ...                            % 点大小
    'CapSize', 12, ...
    'Color',"#1F77B4");                          

xtest = log(p_vec(1:6));
ytest = log(flow_current_node_bet_mean(1:6));
set(gca,"YScale", "log")
set(gca,"XScale", "log")
% ylim([0,1200])
xlabel('$p$',Interpreter='latex',FontSize=24);
ylabel('Current-flow betweenness',Interpreter='latex',FontSize=24);
set(gca, 'FontSize', 22, ...   % 坐标轴字体大小
         'TickDir', 'in', ...
         'TickLength', [0.02 0.02]);     % 刻度线朝里

% figure_name = filefolder_name+"degree_avecfb"+file_network_name+".pdf";
% print(gcf,figure_name, '-dpdf', '-r600');  % 保存为 scatter_plot.pdf




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
% 
% plot(x_fit_normal, y_fit_normal, 'LineWidth', 2,'Color', colors(4))         
% hold on
% 
% 
% ylabel('$f_b(x)$','interpreter','latex','FontSize',30)
% xlabel('$x$','interpreter','latex','FontSize',30);
% ylim([0.0001,10000])
% 
% eqnStr = sprintf('$y = %.4f e^{%.4f x}$', a, b);
% 
% lgd = legend('current-flow betweenness','path betweenness',eqnStr,normStr, 'Interpreter','latex','FontSize',18,'Location','northeast');
% pos = lgd.Position;    % [x y width height]
% pos(1) = pos(1) + 0.014; % 向右挪一点
% pos(2) = pos(2) + 0.02; % 向上挪一点
% lgd.Position = pos;
% set(gca, 'FontSize', 22, ...   % 坐标轴字体大小
%          'TickDir', 'in', ...
%          'TickLength', [0.02 0.02]);     % 刻度线朝里
% % set(gca,"XScale", "log")
% set(gca,"YScale", "log")
% figure_name = filefolder_name+"distribution_fcb_pb"+file_network_name+".pdf";
% print(gcf,figure_name, '-dpdf', '-r600');  % 保存为 scatter_plot.pdf


