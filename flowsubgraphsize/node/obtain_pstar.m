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
    c = (N-1)*p;
    % 定义方程：f(x) = LHS - RHS
    f = @(x) x - (1 - (1 - p * x)^(k - 1) - x .* (1 - x) .* (k - 1) .* p .* (1 - p * x).^(k - 2));
    f = @(x) x - (1-exp(-c*x));

    % 初始区间 [eps, 1-eps] 避免边界不适
    x0 = 0;
    x1 = 1.1;

    % 数值求解
    x = find_all_roots(f, x0, x1);
    if isempty(x)
        x=0;
    end
end

