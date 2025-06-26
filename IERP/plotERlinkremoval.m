clear,clc
N_vec = [20, 50, 100, 200];
p_start_vec = zeros(length(N_vec),1);
countN = 1; 
for N = N_vec
    p_start_vec(countN) = round(log(N)/N,4);
    countN = countN+1;
end

data_mean = zeros(length(N_vec),1);
data_std = zeros(length(N_vec),2);
countN = 1;

for N = N_vec
    p_vec = linspace(p_start_vec(countN), 1, 15);
    p_vec = round(p_vec,4);
    countp = 1;
    for p= p_vec
        filename = sprintf("D:\\data\\flow betweenness\\IERP\\IERP_N%dERp%.4f.txt",N,p);
        results = readmatrix(filename);
        results = results(:,1);
        results = results*(2/(N*(N-1)));
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
    p_vec = linspace(p_start_vec(countplot), 1, 15);
    p_vec = round(p_vec,4);
    errorbar(p_vec, data_mean(:,countplot), data_std(:,countplot), 'o-', 'Color', colors(countplot), 'LineWidth', 4, 'MarkerSize', 10,'CapSize',8);
end


% 图像美化
ax = gca;  % Get current axis
ax.FontSize = 20;  % Set font size for tick label
% xlim([0.01 0.55])
% ylim([0.05 0.25])
% xticks([1 2 3 4])
% xticklabels({'10','20','50','100'})
xlabel('$p$',Interpreter='latex',FontSize=24);
ylabel('$\frac{2}{N(N-1)}(L_H-L_G)$','interpreter','latex',FontSize=30)
lgd = legend({'$N = 20$', '$N = 50$', '$N = 100$', '$N = 200$'}, 'interpreter','latex','Location', 'northeast',FontSize=30);
% lgd.NumColumns = 2;
% set(legend, 'Position', [0.446, 0.73, 0.2, 0.1]);
box on
hold off

picname = sprintf("D:\\data\\flow betweenness\\IERP\\IERP_ER_linkremovalnum.pdf");
exportgraphics(fig, picname,'BackgroundColor', 'none','Resolution', 600);