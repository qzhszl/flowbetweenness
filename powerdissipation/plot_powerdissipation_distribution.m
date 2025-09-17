clear,clc
filefolder_name = "D:\\data\\flow betweenness\\";

n = 100;
p = 0.2;
resname  = sprintf('power_dissipation_N%dp%.2fER.mat',n,p);
filename = filefolder_name+resname;

% load_and_plot_powerdissipation.m
S = load(filename);
results = S.results;

% % 2. 拼接所有实验的 energy
total_energy_path = [];
for r = 1:numel(results)
    total_energy_path = [total_energy_path; results(r).total_SP(:)];
end


total_energy_flow = [];
for r = 1:numel(results)
    total_energy_flow = [total_energy_flow; results(r).total_Flow(:)];
end


% % 3. Figure. 1 plot the distribution of total energy 
figure;
histogram(total_energy_path, 30, 'Normalization', 'pdf'); % 50 bins
% set(gca,"YScale", "log")
hold on

histogram(total_energy_flow, 30, 'Normalization', 'pdf'); % 50 bins

ylabel('$f_e(x)$','interpreter','latex','FontSize',30)
xlabel('$x$','interpreter','latex','FontSize',30);



% % 2. load power dissipation for each link
total_link_energy_path = [];
for r = 1:numel(results)
    total_link_energy_path = [total_link_energy_path; results(r).linkP_SP(:)];
end


total_link_energy_flow = [];
for r = 1:numel(results)
    total_link_energy_flow = [total_link_energy_flow; results(r).linkP_Flow(:)];
end

% % 3. 画直方图 (分布图)
figure;
histogram(total_link_energy_path, 60, 'Normalization', 'pdf'); % 50 bins
% set(gca,"YScale", "log")
hold on

histogram(total_link_energy_flow, 60, 'Normalization', 'pdf'); % 50 bins

ylabel('$f_e(x)$','interpreter','latex','FontSize',30)
xlabel('$x$','interpreter','latex','FontSize',30);



