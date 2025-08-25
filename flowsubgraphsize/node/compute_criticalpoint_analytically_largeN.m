clear,clc
% This .m is for compute the critical point of the  node size of the flow subgraph
% In this case, we use phi(z) = exp(-u(1-z)) 

% Result1: the critical point of the p* is the same as the S

% show the pic of p* given diff p
N = 1000
n=1000
ave_degree = 1:0.2:2.9;
ave_degree_2 = 3:3:20;
ave_degree = [ave_degree,ave_degree_2];
% ave_degree=[2.2,2.3]
p_vals = ave_degree/(N-1);
figure;
i = 1;
colors = lines(length(p_vals)); 
for p = p_vals
    h = show_rootpic_pstar_largenetwork(n,p);
%     legendEntries{i} = sprintf('p=%.3f', p);
    set(h, 'Color', colors(i,:));
    i=i+1;
%     p_star = obtain_pstar_largenetwork(n, p);
end
% legend(legendEntries, 'Location','best');

% Input network size and ave_degree(connection probability p of G_ER)

pc= log(N)/N;
ave_degree = 1:0.1:4.9;
ave_degree_2 = 5:10;
ave_degree = [ave_degree,ave_degree_2];
p_vals = ave_degree/(N-1);
figure;
%Compute the critical point of S
y_vec =zeros(length(p_vals),1);
i=1;
for p = p_vals
    x = obtain_pstar_largenetwork(N, p);
    if length(x)>1
        x = x(3);
    else
        x = x(1);
    end
    y_vec(i) = x;
    i = i+1;
end
plot(ave_degree,y_vec)
hold on
[ave_degree_c,y_vec_c] = obtain_crtical_point(ave_degree,y_vec)

s_vals = zeros(size(p_vals));

for i = 1:length(p_vals)
    s_vals(i) = compute_s_from_largeER(N, p_vals(i));
end

plot(ave_degree, s_vals, 'LineWidth', 2)
legend("P*","S")
[ave_degree_c2,y_vec_c2] = obtain_crtical_point(ave_degree,s_vals)



function x = obtain_pstar_largenetwork(N, p)
    % obtain_pstar: 
    % We denote the probability that following one random link $l$, 
    % a node in the flow subgraph is found in one of its end $l^+$ as $p^*$
    % Inputs:
    %   N - number of nodes
    %   p - link connection probability of ER
    % Output:
    %   x - solution to the fixed-point equation

    c = (N-1)*p; 

    % 定义方程：f(x) = LHS - RHS
    f = @(x) x - (1 - exp(-c*x) - c*x.*(1-x).*exp(-c*x));

    % 初始区间 [eps, 1-eps] 避免边界不适
    eps = 1e-4;
    x0 = 0-eps;
    x1 = 1+eps;

    % 数值求解
    x = find_all_roots(f, x0, x1);
    if isempty(x)
        x=0;
    end
end

function s = compute_s_from_largeER(N, p)
% compute_s_from_largeER: 计算 ER 随机图中的 s = 1 - φ(1 - x) - x φ'(1 - x)
% Inputs:
%   N - 节点数
%   p - 连边概率
% Output:
%   s - 表达式值

    % 第一步：用 PGF 方程求 x = p*
    x = obtain_pstar_largenetwork(N, p);
    if length(x)>1
        x = x(3);
    else
        x = x(1);
    end
    % 第二步：计算 φ(1 - x) 和 φ'(1 - x)
    c = (N-1)*p;
    
    phi = exp(-c*x);
    phi_prime = x*c*phi;

    % 第三步：代入表达式计算 s
    s = 1 - phi - x * phi_prime;
end



function [x_last,y_last] = obtain_crtical_point(x,y)
    tol = 1e-3;  % 容差（认为小于这个值就算“接近0”）
    idx_all = find(abs(y) < tol);
    
    if ~isempty(idx_all)
        idx_last = idx_all(end);   % 最后一个满足条件的点
        x_last = x(idx_last);
        y_last = y(idx_last);
        fprintf('最后一个接近0的点在 x = %.4f, y = %.4e\n', x_last, y_last);
    else
        disp('没有找到接近0的点');
    end  
end


function h = show_rootpic_pstar_largenetwork(n,p)
%     n = 100
%     p = 0.5
    c = (n-1)*p; 
    f2 = @(x) x - (1 - exp(-c*x) - c*x.*(1-x).*exp(-c*x));
    
    p_vals = linspace(-0.01, 1.1, 10000);
    % x_vals = arrayfun(f1, p_vals);
    s_vals = arrayfun(f2, p_vals);
    
    % plot(p_vals, x_vals, 'b', 'LineWidth', 2); hold on;
    h = plot(p_vals, s_vals, 'LineWidth', 2);
    hold on
    plot(zeros(10),linspace(min(s_vals),max(s_vals),10),"Color","black")
    hold on
    plot(linspace(min(p_vals),max(p_vals),10),zeros(10),"Color","black")
    hold on
    legend('p*');
    xlabel('p');
    ylabel('Function value');
%     title('x(p) and s(p) from ER Random Graph');
end