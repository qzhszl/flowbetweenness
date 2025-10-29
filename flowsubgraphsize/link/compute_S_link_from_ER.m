function s = compute_S_link_from_ER(N, p)
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

    % % 第二步：计算 φ(1 - x) 和 φ'(1 - x)
    % k = N - 1;
    % s1 = 1 - x;
    % 
    % phi = (1 - p + p * s1)^k;
    % phi_prime = k * p * (1 - p + p * s1)^(k - 1);
    % 
    % % 第三步：代入表达式计算 s
    % s = 1 - phi - x * phi_prime;
    % try new
    s = x^4;
end
