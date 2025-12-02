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
colors = ["#D08082","#6FB494","#D9B382","#7A7DB1","#62ABC7","#A0A0A0"];

p_vec = linspace(0,1,10);
h_star = plot(p_vec, -5*ones(length(p_vec),1), 'kp', 'MarkerSize', 10, 'LineWidth',3,'MarkerFaceColor', "k");  % 黑星
h_bar  = errorbar(p_vec, -5*ones(length(p_vec),1), zeros(length(p_vec),1), 'o-', 'Color', 'black', 'LineWidth', 1, 'MarkerSize', 14,'CapSize',12);                        % 黑线


for countplot = 1:length(N_vec)
    h_star = plot(p_vec, -5*ones(length(p_vec),1), 'Color', colors(countplot),'LineWidth', 5);  % 黑星
end


for countplot = 1:length(N_vec)
    p_vec = linspace(p_start_vec(countplot), 1, 15);
    p_vec = round(p_vec,4);
    errorbar(p_vec, data_mean(:,countplot), data_std(:,countplot), 'o-', 'Color', colors(countplot), 'LineWidth', 4, 'MarkerSize', 15,'CapSize',10);

    errorbar(p_vec, zeros(length(p_vec),1), zeros(length(p_vec),1), 'p', 'Color', colors(countplot), 'LineWidth', 4, 'MarkerSize', 10,'CapSize',8);
end


% 图像美化
ax = gca;  % Get current axis
ax.FontSize = 30;  % Set font size for tick label
xlim([0 1.02])
ylim([-1.05 0.05])
% xticks([1 2 3 4])
% xticklabels({'10','20','50','100'})
xlabel('$p$',Interpreter='latex',FontSize=36);
ylabel('$(L_H-L_G)/{N \choose 2}$','interpreter','latex',FontSize=36)
lgd = legend({'Fiedler','RGP','$N = 20$', '$N = 50$', '$N = 100$', '$N = 200$'}, 'interpreter','latex','Location', 'southwest',FontSize=30);
lgd.NumColumns = 2;
lgd.ItemTokenSize = [25, 40];
box on
hold off

picname = sprintf("D:\\data\\flow betweenness\\IERP\\IERP_ER_linkremovalnum_compare2alg.pdf");
exportgraphics(fig, picname,'BackgroundColor', 'none','Resolution', 600);

picname = sprintf("D:\\data\\flow betweenness\\IERP\\IERP_ER_linkremovalnum_compare2alg.svg");
set(gcf, 'Renderer', 'painters');   % 强制矢量方式
print(gcf, '-dsvg', picname); % 输出真正的 svg 文件