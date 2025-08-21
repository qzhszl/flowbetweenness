clear,clc
% the simulation results should be not in the largest LCC

N = 100;

pc= log(N)/N
% ave_degree = 1:0.1:4.9;
% ave_degree_2 = 5:10;
% ave_degree = [ave_degree,ave_degree_2];
% p_vals = ave_degree/(N-1);
% p_vals = pc:0.01:1;
% s_vals = zeros(size(p_vals));
% 
% for i = 1:length(p_vals)
%     s_vals(i) = compute_s_from_ER(N, p_vals(i));
% end
% 
% plot(p_vals, s_vals, 'LineWidth', 2)
% xlabel('p')
% ylabel('s')
% title(sprintf('s vs p for ER random graph (N = %d)', N))
% grid on
n = 100
% p = 0.5
% k = n-1
% % f1= @(x) x
% f2 = @(x) x-(1 - (1 - p * x)^(k - 1) - x .* (1 - x) .* (k - 1) .* p .* (1 - p * x).^(k - 2));
% 
% p_vals = linspace(0.001, 1, 2000);
% % x_vals = arrayfun(f1, p_vals);
% s_vals = arrayfun(f2, p_vals);
% 
% % plot(p_vals, x_vals, 'b', 'LineWidth', 2); hold on;
% plot(p_vals, s_vals, 'r--', 'LineWidth', 2);
% legend('x(p)', 's(p)');
% xlabel('p');
% ylabel('Function value');
% title('x(p) and s(p) from ER Random Graph');

ave_degree = 1:0.1:4.9;
ave_degree_2 = 5:10;
ave_degree = [ave_degree,ave_degree_2];
p_vals = ave_degree/(N-1);
y_vec =zeros(length(p_vals),1);
i=1;
for p = p_vals
    x = obtain_pstar(N, p);
    if length(x)>1
        x = x(2);
    else
        x = x(1);
    end
    y_vec(i) = x;
    i = i+1;
end
plot(ave_degree,y_vec)

function x = obtain_pstar(N, p)
% obtain_pstar: 
% We denote the probability that following one random link $l$, 
% a node in the flow subgraph is found in one of its end $l^+$ as $p^*$
% Inputs:
%   N - number of nodes
%   p - link connection probability of ER
% Output:
%   x - solution to the fixed-point equation

    k = N - 1;

    % 定义方程：f(x) = LHS - RHS
    f = @(x) x - (1 - (1 - p * x)^(k - 1) - x .* (1 - x) .* (k - 1) .* p .* (1 - p * x).^(k - 2));

    % 初始区间 [eps, 1-eps] 避免边界不适
    x0 = 0;
    x1 = 1.1;

    % 数值求解
    x = find_all_roots(f, x0, x1);
    if isempty(x)
        x=0;
    end
end


function s = compute_s_from_ER(N, p)
% compute_s_from_ER: 计算 ER 随机图中的 s = 1 - φ(1 - x) - x φ'(1 - x)
% Inputs:
%   N - 节点数
%   p - 连边概率
% Output:
%   s - 表达式值

    % 第一步：用 PGF 方程求 x = p*
    x = obtain_pstar(N, p);
    if length(x)>1
        x = x(2);
    else
        x = x(1);
    end

    % 第二步：计算 φ(1 - x) 和 φ'(1 - x)
    k = N - 1;
    s1 = 1 - x;
    
    phi = (1 - p + p * s1)^k;
    phi_prime = k * p * (1 - p + p * s1)^(k - 1);

    % 第三步：代入表达式计算 s
    s = 1 - phi - x * phi_prime;
end



