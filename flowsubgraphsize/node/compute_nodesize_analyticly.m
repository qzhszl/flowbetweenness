clear,clc
N = 100;

pc= log(N)/N;
ave_degree = 1:0.1:4.9;
ave_degree_2 = 5:10;
ave_degree = [ave_degree,ave_degree_2];
p_vals = ave_degree/(N-1);

s_vals = zeros(size(p_vals));

for i = 1:length(p_vals)
    s_vals(i) = compute_s_from_ER(N, p_vals(i));
end

plot(ave_degree, s_vals, 'LineWidth', 2)
xlabel('p')
ylabel('s')
title(sprintf('s vs p for ER random graph (N = %d)', N))
grid on



function s = compute_s_from_ER(N, p)
% compute_s_from_ER: 计算 ER 随机图中的 s = 1 - φ(1 - x) - x φ'(1 - x)
% Inputs:
%   N - 节点数
%   p - 连边概率
% Output:
%   s - 表达式值

    % 第一步：用 PGF 方程求 x = p*
    x = obtain_pstar(N, p);

    % 第二步：计算 φ(1 - x) 和 φ'(1 - x)
    k = N - 1;
    s1 = 1 - x;
    
    phi = (1 - p + p * s1)^k;
    phi_prime = k * p * (1 - p + p * s1)^(k - 1);

    % 第三步：代入表达式计算 s
    s = 1 - phi - x * phi_prime;
end




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
    x0 = 1e-6;
    x1 = 1 - 1e-6;

    % 数值求解
    try
        x = fzero(f, [x0, x1]);
    catch
        warning('fzero failed to converge. Returning NaN.');
        x = NaN;
    end
end



