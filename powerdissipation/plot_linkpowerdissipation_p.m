clear,clc
filefolder_name = "D:\\data\\flow betweenness\\";

n = 200;
% p_vec = [0.08,0.1,0.15,0.2,0.28,0.39,0.5,0.66];
% p_vec=[0.15,0.28,0.39,0.66,0.88]
% p_vec = [0.05,0.1,0.2,0.5]
if n == 200
    p_vec = [0.03,0.04,0.06,0.11,0.15,0.28,0.39,0.66,0.88];
end

ave_link_power_vec = zeros(length(p_vec),1);
std_link_power_vec = zeros(length(p_vec),1);

count = 0;
for p = p_vec
    count = count+1;
    resname  = sprintf('power_dissipation_N%dp%.2fER.mat',n,p);
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

% errorbar(p_vec,ave_link_power_vec,std_link_power_vec,"-o",'LineWidth', 4, 'MarkerSize', 12,color = "#D08082")
plot(p_vec,ave_link_power_vec,"-o",'LineWidth', 4, 'MarkerSize', 12,color = "#D08082")
hold on

% curve fit................
xdata = p_vec(:);
ydata = ave_link_power_vec(:);
% 定义幂律拟合模型：y = a*x^b
X = log(xdata(3:length(xdata)));
Y = log(ydata(3:length(xdata)));

coeff = polyfit(X, Y, 1);   % 拟合直线 Y = b*X + log(a)

b = coeff(1);
loga = coeff(2);
a = exp(loga);

plot(xdata, a*xdata.^b, '--','Color','#7A7DB1','LineWidth',5)


set(gca,"XScale", "log")
set(gca,"YScale", "log")

xlabel('$p$','interpreter','latex','FontSize',30)
ylabel('$E[E_l]$','interpreter','latex','FontSize',30);


ax = gca;  % Get current axis
ax.FontSize = 20;  % Set font size for tick label
% xlim([0.01 0.55])
% ylim([0.05 0.25])
% xticks([1 2 3 4])
% xticklabels({'10','20','50','100'})
lgd = legend({'link power dissipation', sprintf('fit: $y = %.2g x^{%.3g}$', a, b)}, 'interpreter','latex','Location', 'best',FontSize=30);
% lgd.NumColumns = 2;
% set(legend, 'Position', [0.446, 0.73, 0.2, 0.1]);
box on
% set(fig, 'Color', 'none');              % figure 背景透明
% set(gca,  'Color', 'none');             % 坐标轴区域背景透明
hold off

picname = sprintf("D:\\data\\flow betweenness\\power_dissipation\\link_power_distribution_N%d_withdiffp.pdf",n);
% exportgraphics(fig, picname,'BackgroundColor', 'none','Resolution', 600);
print(fig, picname, '-dpdf', '-r600', '-bestfit');
    
 



