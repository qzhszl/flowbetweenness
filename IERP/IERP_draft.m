clear,clc

% this .m is to verify Inverse effecitve resitance problem:
% Given an effective resitance matrix D, we are planning to obtain a network
% whose effective resitance is similar with given demand matrix


% the main idea is :1. start with a complete graph whose link weight = demand
% 2. R = \Omega - W_tilde
% 3. find i,j where reaches the largest R
% 4. remove that link
% 5. renormalize the link
% 6. update the \Omega
% 7. repeat 2-6 until Diff between Omega and the given demand is minimum: Return lastly removed link

N = 1000
p_start_vec = zeros(4,1);
count = 1; 
p = 0.3

simutimes = 1;
count =1;

result = zeros(simutimes,4);
p_vec = linspace(p_start_vec(count), 1, 15);
p_vec = round(p_vec,4);
    
% 1. generate a graph
% _________________________________________________________________________
% (b) ER:
A_input = GenerateERfast(N,p,10);
% check connectivity
connect_flag = network_isconnected(A_input);
while ~connect_flag
    A_input = GenerateERfast(N,p,10);
    % check connectivity
    connect_flag = network_isconnected(A_input);
end
% A_input = [0,0.5,0
%     0.5 0 1
%     0 1 0
%     ]
% Omega = EffectiveResistance(A_input)
% A = ISPP_tree(Omega)
% diff = find(abs(A-A_input)>0.000001)

% 2. run simulations
A_input(A_input ~= 0) = 1 ./ A_input(A_input ~= 0);
% [L_add_output,L_ouput,L_comm_output,Norm_output] = experiment_on_ER(A_input);

Omega = EffectiveResistance(A_input);
A = ISPP_tree(Omega);
find(abs(A - A_input)>0.000001)



function A = ISPP_tree(Omega)
    m = size(Omega,1);
    u = ones(m,1);
    p = 1/(u.'*inv(Omega)*u)*inv(Omega)*u;
    Q_tilde = 2*(u.'*inv(Omega)*u)*p*p.'-2*inv(Omega);
    A = diag(diag(Q_tilde)) - Q_tilde;
    A = round(A, 10);
    % A(A ~= 0) = 1 ./ A(A ~= 0);
end



function [output_Atilde,output_Omega] = IERP_2(D)
    N = size(D,1);
    Input_Omega = D
    W_tilde  = Input_Omega;
    W_tilde(W_tilde ~= 0) = 1 ./ W_tilde(W_tilde ~= 0);
    Omega_new = EffectiveResistance(W_tilde)  ;          % Compute the effective resistance Omega
    
    
    diff_change = sum(sum((abs(D - Omega_new))./(D+(D==0))))/(N*(N-1));

    % val = 1;
    flag = 1;
    A = (W_tilde > 0);
    % Gnow = graph(W_tilde,"upper");
    while(flag==1)                     % Remove links one by one until we exceed the constraints                            
        previous_change = diff_change;
        % method 1 R = A.*((Omega_new+eye(N)).^-1 - W_tilde)*(D-Omega_new) 
        R = A.*((Omega_new+eye(N)).^-1 - W_tilde).*(D-Omega_new);     % Compute R
        [val,~] = max(max(R));                              % Identify the maximum element
        [row,col] = find(R == val);                        % Identify the link
        A(row(1),col(1)) = 0; A(col(1),row(1)) = 0;        % Remove the link
        
        %test: update here
%         
%         W_tilde = A.*D
%         W_tilde(W_tilde ~= 0) = 1 ./ W_tilde(W_tilde ~= 0)     % Compute W tilde
%         
%         Omega_new = EffectiveResistance(W_tilde)              % Update the shortest path weight matrix
%         
%         alpha = alpha_l1_global2(Omega_new,D)
%         
% 
%         Omega_new = alpha*Omega_new
%         W_tilde = 1/alpha* W_tilde
        

        diff_change = sum(sum((abs(D - Omega_new))./(D+(D==0))))/(N*(N-1));

        if abs(diff_change) > abs(previous_change)
            flag=0;
        end

    end
    
    if(row(1) ~= col(1))
        A(row(1),col(1)) = 1; A(col(1),row(1)) = 1;         % Return lastly removed link
        W = A.*D;
    %     W = 2*1/sum(sum(A.*D)).*A.*D;                       % Update the weighted adjacency matrix
    end
    % output_Atilde = W;
    output_Omega = EffectiveResitance_withinverseA(W);
    alpha = alpha_l1_global2(output_Omega,D)
    alpha2 = alpha_l1_global3(output_Omega,D)
    output_Omega = alpha*output_Omega
    

%     output_Atilde2 = ISPP_tree(output_Omega);
    % output_Atilde(output_Atilde ~= 0) = 1 ./ output_Atilde(output_Atilde ~= 0);
    % output_Atilde
    output_Atilde = W;
    % output_Atilde2(output_Atilde2 ~= 0) = 1 ./ output_Atilde2(output_Atilde2 ~= 0);
    output_Atilde = alpha*output_Atilde
%     diff = find(abs(output_Atilde2-output_Atilde)>0.00001)
end



function [L_add_output,L_ouput,L_comm_output_ratio,Norm_output] = experiment_on_ER(A_input)
    % Input_Omega = EffectiveResitance_withinverseA(A_input);
    Input_Omega = EffectiveResistance(A_input);
    % Generate a demand matrix: the effective resistance matrix of the
    % input network
    D = Input_Omega;
    N = size(D,1);
    tic
    [output_Atilde,output_Omega] = IERP(D);
    t3 = toc;
    tic
    [output_Atilde3,output_Omega3] = IERP_test(D);
   

    [output_Atilde2,output_Omega2] = IERP(D);
    output_Omega2
    EffectiveResistance(output_Atilde3)
    EffectiveResistance(output_Atilde2)


    
    t4 = toc;
    data3diff = find(abs(output_Atilde2-output_Atilde)>0.00001);
    data4diff = find(abs(output_Omega2-output_Omega)>0.00001);
    

    % Store the results
    % 1. The number of links added in the graph
    L_add_output = 0.5*(nnz(output_Atilde)-nnz(A_input));
    % 2. The number of links in the obtained graph
    L_ouput = 0.5*nnz(output_Atilde);  
    % 3. The number of common links between two graphs
    L_comm_output_ratio = nnz(A_input.*output_Atilde)/nnz(output_Atilde); 
    % 4. The norm of common links between two graphs
    Norm_output = sum(sum((abs(D - output_Omega))./(D+(D==0))))/(N*(N-1))
    Norm_output_2 = sum(sum((abs(D - output_Omega2))./(D+(D==0))))/(N*(N-1))
    Norm_output_3 = sum(sum((abs(D - output_Omega3))./(D+(D==0))))/(N*(N-1))
end



function A = GenerateERfast(n,p,weighted)
    A = rand(n,n) < p;
    A = triu(A,1);
    if weighted == 0
     
    elseif weighted == 1
        linkweight_matrix = rand(n,n);
        A = A.*linkweight_matrix;
    else
        linkweight_matrix = randi(weighted,n,n);
        A = A.*linkweight_matrix;
    end
    
    A = A + A';
end


function [isConnected] = network_isconnected(adj)
    G = graph(adj);
    components = conncomp(G);
    % 判断图是否连通
    isConnected = (max(components) == 1);
end

function alpha = alpha_l1_global2(A, D)

    r = D(A>0)./A(A>0);
    w = A(A>0)./D(A>0);
    
    % Sort and find weighted median
    [rsort, idx] = sort(r(:));
    w = w(idx);
    cw = cumsum(w)/sum(w);
    alpha = rsort(find(cw>=0.5, 1));
end

function alpha = alpha_l1_global3(A, D)

    r = D(A>0)./A(A>0);
    
    alpha = mean(r);
end

function [output_Atilde,output_Omega] = IERP_test(D)
    N = size(D,1);
    Input_Omega = D;
    W_tilde  = Input_Omega;
    W_tilde(W_tilde ~= 0) = 1 ./ W_tilde(W_tilde ~= 0);
    Omega_new = EffectiveResistance(W_tilde);            % Compute the effective resistance Omega
    diff_change = sum(sum((abs(D - Omega_new))./(D+(D==0))))/(N*(N-1));

    % val = 1;
    flag = 1;
    A = (W_tilde > 0);
    % Gnow = graph(W_tilde,"upper");
    while(flag==1)                     % Remove links one by one until we exceed the constraints                            
        previous_change = diff_change;
        % method 1 R = A.*((Omega_new+eye(N)).^-1 - W_tilde)*(D-Omega_new) 
        R = A.*((Omega_new+eye(N)).^-1 - W_tilde).*(D-Omega_new);     % Compute R
        [val,~] = max(max(R));                              % Identify the maximum element
        [row,col] = find(R == val);                        % Identify the link
        A(row(1),col(1)) = 0; A(col(1),row(1)) = 0;        % Remove the link
        W_tilde = A.*D;
        W_tilde(W_tilde ~= 0) = 1 ./ W_tilde(W_tilde ~= 0);      % Compute W tilde
        % Gnow = graph(W_tilde,"upper");
        Omega_new = EffectiveResistance(W_tilde);               % Update the shortest path weight matrix
        
        diff_change = sum(sum((abs(D - Omega_new))./(D+(D==0))))/(N*(N-1));

        if abs(diff_change) > abs(previous_change)
            flag=0;
        end    
    end
    
    if(row(1) ~= col(1))
        A(row(1),col(1)) = 1; A(col(1),row(1)) = 1;         % Return lastly removed link
        W = A.*D;
    %     W = 2*1/sum(sum(A.*D)).*A.*D;                       % Update the weighted adjacency matrix
    end
    W(W ~= 0) = 1 ./ W(W ~= 0);
    output_Atilde = W;
  
    output_Omega = EffectiveResistance(output_Atilde);
    alpha = alpha_l1_global_para(output_Omega,D);
    output_Omega = alpha*output_Omega;
    output_Atilde = 1/alpha*output_Atilde;
end



function alpha = alpha_l1_global(A, D)
% ALPHA_L1_GLOBAL  找到标量 alpha 以最小化 || alpha*A - D ||_1
%   alpha = alpha_l1_global(A,D)
% 返回加权中位数解，处理 A==0 的条目（这部分不影响 alpha）。
%
% 若所有 A(:)==0，则任意 alpha 都相同，这里返回 0 并发出警告。

a = A(:);
d = D(:);

% 只保留 a~=0 的项（a==0 的项与 alpha 无关）
mask = (a ~= 0);
if ~any(mask)
    warning('All entries of A are zero. Objective independent of alpha. Return alpha = 0.');
    alpha = 0;
    return;
end

b = d(mask) ./ a(mask);   % 比值 d_i / a_i
w = abs(a(mask));         % 权重 |a_i|

% 排序 b 并把权重一并排序
[bs, idx] = sort(b);
ws = w(idx);

W = sum(ws);
cumw = cumsum(ws);

% 找到第一个使 cumw >= W/2 的索引
k = find(cumw >= W/2, 1, 'first');

% 若恰好等于 W/2 且下一个 b 不相同，则任意在 [bs(k), bs(k+1)] 都是最优解。
% 这里我们返回中点以保证确定性。
if cumw(k) == W/2 && k < numel(bs)
    % 若下一个不同则取区间中点，否则等于时取 bs(k)
    if bs(k) ~= bs(k+1)
        alpha = 0.5*(bs(k) + bs(k+1));
    else
        alpha = bs(k);
    end
else
    alpha = bs(k);
end
end
