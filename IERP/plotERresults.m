clear,clc

N_vec = [20,50,100,200];
% N_vec = [20,50]
p_start_vec = zeros(length(N_vec),1);
countN = 1;

for N = N_vec
    x = log(N)/N;
    y = ceil(x * 1e4) / 1e4;  % round 4 decimal
    p_start_vec(countN) = y;
    countN = countN+1;
end

data_mean = zeros(length(N_vec),1);
data_std = zeros(length(N_vec),1);
countN = 1;

for N = N_vec
    % p_vec = linspace(p_start_vec(countN), 1, 15);
    % p_vec = round(p_vec,4);
    
    % more points in the beginning
    p_vec = linspace(p_start_vec(countN), 1, 15);
    % 前两个点
    p1 = p_vec(1);
    p2 = p_vec(2);  
    % 在 p1 和 p2 之间插入两个点
    extra_points = linspace(p1, p2, 4);  % 生成4个点
    extra_points = extra_points(2:3);    % 去掉第一个和最后一个（原本已有）
    % 合并
    p_vec = [p_vec(1), extra_points, p_vec(2:end)];
    p_vec = round(p_vec,4);
    if N ==100
        p_vec = [p_vec(1), 0.0500, p_vec(2:end)];
    elseif N ==200
        p_vec = [p_vec(1), 0.0300, p_vec(2:end)];
    end
    p_vec
    % p_vec =p_vec(1:6)

    countp = 1;
    for p= p_vec
        filename = sprintf("D:\\data\\flow betweenness\\IERP\\IERP_N%dERp%.4f.txt",N,p);
        results = readmatrix(filename);
        results = results(:,4);
        mean_values = mean(results);
        std_values = std(results);
        data_mean(countp,countN) = mean_values;
        data_std(countp,countN) = std_values;
        countp =countp+1;
    end
    countN = countN+1;
end

fig = figure; 
fig.Position = [100 100 900 600]; 
hold on;
colors = ["#D08082", "#C89FBF", "#62ABC7", "#7A7DB1", "#6FB494", "#D9B382"];


for countplot = 1:length(N_vec)
    N = N_vec(countplot);
    % p_vec = linspace(p_start_vec(countplot), 1, 15);
    % p_vec = round(p_vec,4);


    % more points in the beginning
    p_vec = linspace(p_start_vec(countplot), 1, 15);
    % 前两个点
    p1 = p_vec(1);
    p2 = p_vec(2);  
    % 在 p1 和 p2 之间插入两个点
    extra_points = linspace(p1, p2, 4);  % 生成4个点
    extra_points = extra_points(2:3);    % 去掉第一个和最后一个（原本已有）
    % 合并
    p_vec = [p_vec(1), extra_points, p_vec(2:end)];
    p_vec = round(p_vec,4);
    if N ==100
        p_vec = [p_vec(1), 0.0500, p_vec(2:end)];
    elseif N ==200
        p_vec = [p_vec(1), 0.0300, p_vec(2:end)];
    end
    % if N<100
    %     p_vec = p_vec
    % end
    y_plot = data_mean(:,countplot);
    std_plot = data_std(:,countplot);
    errorbar(p_vec, y_plot(1:length(p_vec)), std_plot(1:length(p_vec)), 'o-', 'Color', colors(countplot), 'LineWidth', 4, 'MarkerSize', 10,'CapSize',8);
end


% 图像美化
ax = gca;  % Get current axis
ax.FontSize = 30;  % Set font size for tick label
% xlim([0.01 0.55])
% ylim([0.05 0.3])
% xticks([1 2 3 4])
% xticklabels({'10','20','50','100'})
xlabel('$p$',Interpreter='latex',FontSize=40);
% ylabel('$\frac{1}{N(N-1)}\sum_i\sum_j\frac{|d_{ij}-\omega_{ij}|}{d_{ij}}$','interpreter','latex',FontSize=40)
ylabel('$\|D-\Omega\|$','interpreter','latex',FontSize=40)
% lgd = legend({'$N = 20$', '$N = 50$', '$N = 100$', '$N = 200$'}, 'interpreter','latex','Location', 'northwest',FontSize=30);
% lgd.NumColumns = 2;
% set(legend, 'Position', [0.446, 0.73, 0.2, 0.1]);
box on
hold off

picname = sprintf("D:\\data\\flow betweenness\\IERP\\IERP_ER_norm.pdf");
exportgraphics(fig, picname,'BackgroundColor', 'none','Resolution', 600);