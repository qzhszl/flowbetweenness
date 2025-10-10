clear,clc

weighted_flag =0;
color_count = 1;


fig = figure; 
fig.Position = [100 100 900 600];
nvec = [100,200,1000]
for n = nvec
    if n==1000
        s = plot_elvsp(n,weighted_flag,1,color_count);
    else
        plot_elvsp(n,weighted_flag,0,color_count)
    end
    color_count = color_count+1;
end


set(gca,"XScale", "log")
set(gca,"YScale", "log")

xlabel('$p$','interpreter','latex','FontSize',50)
ylabel('$E[\Lambda_l]$','interpreter','latex','FontSize',50);


ax = gca;  % Get current axis
ax.FontSize = 30;  % Set font size for tick label
ylim([1e-12 5e-2])
% xlim([0.05 0.25])
% xticks([1 2 3 4])
% xticklabels({'10','20','50','100'})
if length(nvec)>1
    lgd = legend({'simultion:$\Lambda_l$, $N=100$', 'simultion:$\Lambda_l$, $N=200$', 'simultion:$\Lambda_l$, $N=1000$',s}, 'interpreter','latex','Location', 'best',FontSize=26);
    lgd.Position = [0.39, 0.29, 0.15, 0.15];
    lgd.ItemTokenSize = [40,10];
else
    lgd = legend({sprintf('simultion:$\\Lambda_l$, $N=%d$',n), s}, 'interpreter','latex','Location', 'best',FontSize=30);
end
% lgd.NumColumns = 2;
% set(legend, 'Position', [0.446, 0.73, 0.2, 0.1]);
box on
% set(fig, 'Color', 'none');              % figure 背景透明
% set(gca,  'Color', 'none');             % 坐标轴区域背景透明
hold off
if weighted_flag ==0
    picname = sprintf("D:\\data\\flow betweenness\\power_dissipation\\link_power_N%d_withdiffp_unweighted.pdf",n);
else
    picname = sprintf("D:\\data\\flow betweenness\\power_dissipation\\link_power_N%d_withdiffp.pdf",n);
end
exportgraphics(fig, picname,'BackgroundColor', 'none','Resolution', 600);
% print(fig, picname, '-dpdf', '-r600', '-bestfit');



function s = plot_elvsp(n,weighted_flag,curve_fit_flag,color_count)
filefolder_name = "D:\\data\\flow betweenness\\";
colors = ["#D08082", "#C89FBF", "#62ABC7", "#7A7DB1", "#6FB494", "#D9B382"];
% p_vec = [0.08,0.1,0.15,0.2,0.28,0.39,0.5,0.66];
% p_vec=[0.15,0.28,0.39,0.66,0.88]
% p_vec = [0.05,0.1,0.2,0.5]
if n == 200
    % p_vec = [0.03,0.04,0.06,0.11,0.15,0.28,0.39,0.66,0.88];
    p_vec = [0.03, 0.05, 0.07, 0.11, 0.16, 0.25, 0.38, 0.58, 0.88];
elseif n ==100
    % p_vec = [0.05,0.06,0.08,0.11,0.15,0.28,0.39,0.66,0.88];
    p_vec = [0.05,0.07, 0.10, 0.15, 0.21, 0.30, 0.43, 0.61, 0.88];
elseif n == 1000
    p_vec = [0.007,0.013,0.0269,0.0529,0.1037,0.2034,0.3990,0.88];
end

ave_link_power_vec = zeros(length(p_vec),1);
std_link_power_vec = zeros(length(p_vec),1);

count = 0;
for p = p_vec
    count = count+1;
    if weighted_flag ==0
        resname  = sprintf('power_dissipation_N%dp%.2fER_unweighted.mat',n,p);
        filename = filefolder_name+resname;
        % load_and_plot_powerdissipation.m
        try
            S = load(filename);
        catch
            resname  = sprintf('power_dissipation_N%dp%.4fER_unweighted.mat',n,p);
            filename = filefolder_name+resname;
            S = load(filename);
        end
        results = S.results;
    else
        resname  = sprintf('power_dissipation_N%dp%.2fER.mat',n,p);
        filename = filefolder_name+resname;
        % load_and_plot_powerdissipation.m
        S = load(filename);
        results = S.results;
    end
    
    
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
errorbar(p_vec,ave_link_power_vec,std_link_power_vec,"-o",'LineWidth', 4, 'MarkerSize', 12,color = colors(color_count))
% plot(p_vec,ave_link_power_vec,"-o",'LineWidth', 4, 'MarkerSize', 12,color = colors(color_count))
hold on

if curve_fit_flag == 1
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
    plot(xdata, a*xdata.^b, '--','Color','#D9B382','LineWidth',5)
    s = sprintf('fit: $E[\\Lambda_l] = %.2e p^{%.3g}$', a, b);
    s = regexprep(s, 'e([+-]?\d+)', '\\times 10^{$1}');
    s = regexprep(s, '\^\{\+(\d+)\}', '^{$1}') ; % 去掉正号
    s = regexprep(s, '\^\{(-?)0+(\d+)\}', '^{$1$2}');
end
end

    
 



