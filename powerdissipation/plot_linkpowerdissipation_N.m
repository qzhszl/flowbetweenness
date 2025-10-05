clear,clc
filefolder_name = "D:\\data\\flow betweenness\\";

weigthed_flag = 0;
n_vec = [77,100,120,287,444,686,1062,1643];
n_vec = [77,120,287,444,686,1062];


ave_link_power_vec = zeros(length(n_vec),1);
std_link_power_vec = zeros(length(n_vec),1);

count = 0;
for n = n_vec
    if  n == 200
        p =0.11;
    else
        p =0.1;
    end
    count = count+1;
    if weigthed_flag==0
        resname  = sprintf('power_dissipation_N%dp%.2fER_unweighted.mat',n,p);
    else
        resname  = sprintf('power_dissipation_N%dp%.2fER.mat',n,p);
    end
    filename = filefolder_name+resname;
    
    % load_and_plot_powerdissipation.m
    S = load(filename);
    results = S.results;
    
    % % 2. (1)load and combine data
    total_energy_path = [];
    for r = 1:numel(results)
        total_energy_path = [total_energy_path; results(r).total_SP(:)];
    end
     
    total_energy_flow = [];
    for r = 1:numel(results)
        total_energy_flow = [total_energy_flow; results(r).total_Flow(:)];
    end
    
    % % (2). load power dissipation for each link
    link_energy_path = [];
    for r = 1:numel(results)
        link_energy_path = [link_energy_path; results(r).linkP_SP(:)];
    end
    link_energy_path = link_energy_path./nchoosek(n,2);
    
    link_energy_flow = [];
    for r = 1:numel(results)
        link_energy_flow = [link_energy_flow; results(r).linkP_Flow(:)];
    end
    link_energy_flow = link_energy_flow./nchoosek(n,2);
    

    ave_link_power_vec(count) = mean(link_energy_flow);
    std_link_power_vec(count)  =std(link_energy_flow);
end


% % 3. 画直方图 (分布图)
fig = figure; 
fig.Position = [100 100 900 600]; 
colors = ["#D08082", "#C89FBF", "#62ABC7", "#7A7DB1", "#6FB494", "#D9B382"];

% errorbar(n_vec,ave_link_power_vec,std_link_power_vec,"-o",'LineWidth', 4, 'MarkerSize', 12,color = "#D08082")
plot(n_vec,ave_link_power_vec,"-o",'LineWidth', 4, 'MarkerSize', 12,color = "#D08082")
hold on

% curve fit................
xdata = n_vec(:);
ydata = ave_link_power_vec(:);
% 定义幂律拟合模型：y = a*x^b
X = log(xdata(4:length(xdata)));
Y = log(ydata(4:length(xdata)));

coeff = polyfit(X, Y, 1);   % 拟合直线 Y = b*X + log(a)

b = coeff(1);
loga = coeff(2);
a = exp(loga);

plot(xdata, a*xdata.^b, '--','Color','#7A7DB1','LineWidth',5)


set(gca,"XScale", "log")
set(gca,"YScale", "log")

xlabel('$N$','interpreter','latex','FontSize',50)
ylabel('$E[\Lambda_l]$','interpreter','latex','FontSize',50);


ax = gca;  % Get current axis
ax.FontSize = 30;  % Set font size for tick label
% xlim([0.01 0.55])
% ylim([0.05 0.25])
% xticks([1 2 3 4])
% xticklabels({'10','20','50','100'})
lgd = legend({'simultion:$\Lambda_l$', sprintf('fit: $E[\\Lambda_l] = %.2g N^{%.3g}$', a, b)}, 'interpreter','latex','Location', 'best',FontSize=30);
% lgd.NumColumns = 2;
% set(legend, 'Position', [0.446, 0.73, 0.2, 0.1]);
box on
% set(fig, 'Color', 'none');              % figure 背景透明
% set(gca,  'Color', 'none');             % 坐标轴区域背景透明
hold off
if weigthed_flag ==1
    picname = sprintf("D:\\data\\flow betweenness\\power_dissipation\\link_power_distribution_p%.2f_withdiffN.pdf",p);
else
    picname = sprintf("D:\\data\\flow betweenness\\power_dissipation\\link_power_distribution_p%.2f_withdiffN_unweighted.pdf",p);
end
% exportgraphics(fig, picname,'BackgroundColor', 'none','Resolution', 600);
print(fig, picname, '-dpdf', '-r600', '-bestfit');
    
 



