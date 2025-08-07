function roots = find_all_roots(f, a, b, N)
% 自动搜索函数 f 在区间 [a, b] 内的所有根
% f : 待求根的函数句柄 @(x) ...
% a, b : 搜索区间
% N : 网格点数量（越大越精细，默认10000）

    if nargin < 4
        N = 10000;
    end

    x_grid = linspace(a, b, N);
    f_vals = arrayfun(f, x_grid);

    roots = [];
    for i = 1:length(x_grid) - 1
        if f_vals(i) * f_vals(i+1) < 0
            try
                x_root = fzero(f, [x_grid(i), x_grid(i+1)]);
                % 避免重复根（由于数值误差）
                if isempty(roots) || all(abs(roots - x_root) > 1e-6)
                    roots(end+1) = x_root; %#ok<AGROW>
                end
            catch
                % fzero失败，跳过
            end
        end
    end
end


