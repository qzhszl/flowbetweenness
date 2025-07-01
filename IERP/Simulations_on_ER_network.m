clear, clc
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

N_vec = [10, 20, 50, 100, 200];
N_vec = [10];
p_start_vec = zeros(4,1);
count = 1; 
for N = N_vec
    p_start_vec(count) = round(log(N)/N,4);
    count = count+1;
end

simutimes = 1000;
count =1;
for N = N_vec
    N
    result = zeros(simutimes,4);
    p_vec = linspace(p_start_vec(count), 1, 15);
    p_vec = round(p_vec,4);
    for p= p_vec
        p
        for simu_time = 1:simutimes
            simu_time
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

            % 2. run simulations
            [L_add_output,L_ouput,L_comm_output,Norm_output] = experiment_on_ER(A_input);
            result(simu_time,:) = [L_add_output,L_ouput,L_comm_output,Norm_output];
        end
        filename = sprintf("D:\\data\\flow betweenness\\IERP\\IERP_N%dERp%.4f_weight01.txt",N,p);
        writematrix(result,filename)
    end
    count = count+1;
end




function [L_add_output,L_ouput,L_comm_output_ratio,Norm_output] = experiment_on_ER(A_input)
    Input_Omega = EffectiveResitance_withinverseA(A_input);
    % Generate a demand matrix: the effective resistance matrix of the
    % input network
    D = Input_Omega;
    N = size(D,1);
    [output_Atilde,output_Omega] = IERP(D);
    % Store the results
    % 1. The number of links added in the graph
    L_add_output = 0.5*(nnz(output_Atilde)-nnz(A_input));
    % 2. The number of links in the obtained graph
    L_ouput = 0.5*nnz(output_Atilde);  
    % 3. The number of common links between two graphs
    L_comm_output_ratio = nnz(A_input.*output_Atilde)/nnz(output_Atilde); 
    % 4. The norm of common links between two graphs
    Norm_output = sum(sum((abs(D - output_Omega))./(D+(D==0))))/(N*(N-1));
end




function T = generate_a_tree(N,minlinkweight,maxlinkweight)
% 生成完全连接的随机加权图
W = randi([minlinkweight,maxlinkweight], N, N);  % 生成 1-10 之间的随机整数
W = triu(W,1);            % 仅保留上三角部分以避免重复
W = W + W';               % 生成对称矩阵，表示无向图
% 计算最小生成树
G = graph(W);             % 生成图
T = minspantree(G);       % 计算最小生成树
T.Edges.Weight = randi([minlinkweight,maxlinkweight], numedges(T), 1);
end

function A = ISPP_tree(Omega)
m = size(Omega,1);
u = ones(m,1);
p = 1/(u.'*inv(Omega)*u)*inv(Omega)*u;
Q_tilde = 2*(u.'*inv(Omega)*u)*p*p.'-2*inv(Omega);
A = diag(diag(Q_tilde)) - Q_tilde;
A = round(A, 10);
A(A ~= 0) = 1 ./ A(A ~= 0);
end

function Omega = EffectiveResitance_withinverseA(A)
% the resitance is the inverse of the link weight
    Aforomega = A;
    Aforomega(Aforomega ~= 0) = 1 ./ Aforomega(Aforomega ~= 0);
    Omega = EffectiveResistance(Aforomega);
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
