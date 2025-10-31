clear,clc

N = 100;

p_vec = linspace(0.7,1,100); 
y_vec = zeros(length(p_vec),1);
count = 1;
for p  = [0.6]
    y_vec(count) = obtain_ex(N, p);
    count = count+1;
end
plot(p_vec,y_vec)


function x_solution = obtain_ex(N, p)
    % obtain_pstar: 
    EL = nchoosek(N, 2) * p;

    % 3. 定义函数 f(x)
    % f(x) = x - EL + ((x/EL*p)^2+(1-x/EL*p)^2)^(N-2)*p;
    f = @(x) x - EL + (((x/EL*p).^2) + ((1 - x/EL*p).^2)).^(N-2)* nchoosek(N-2, 2)* p;
    
    % 4. 使用 fzero 求解
    % 选取一个初始猜测值，例如 EL
    x0 = EL; 
    x_solution = fzero(f, x0);

    figure;
    fplot(f, [0, EL]); % 绘制范围从 0 到 2*EL
    hold on;
    plot(x_solution, 0, 'ro', 'MarkerSize', 8); % 标记找到的解
    yline(0, 'k--'); % 绘制 x 轴
    xlabel('x');
    ylabel('f(x)');
    title('Function f(x) Plot');
    grid on;

end
