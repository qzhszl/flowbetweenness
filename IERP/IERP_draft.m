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

N_vec = [100];
p_start_vec = zeros(4,1);
count = 1; 
for N = N_vec
    p_start_vec(count) = round(log(N)/N,4);
    count = count+1;
end

simutimes = 1;
count =1;
for N = N_vec
    N
    result = zeros(simutimes,4);
    p_vec = linspace(p_start_vec(count), 1, 15);
    p_vec = round(p_vec,4);
    for p= p_vec
        p
        for simu_time = 1:simutimes
%             simu_time
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
%         writematrix(result,filename)
    end
    count = count+1;
end




function [L_add_output,L_ouput,L_comm_output_ratio,Norm_output] = experiment_on_ER(A_input)
    Input_Omega = EffectiveResitance_withinverseA(A_input);
    % Generate a demand matrix: the effective resistance matrix of the
    % input network
    D = Input_Omega;
    N = size(D,1);
    tic
    [output_Atilde,output_Omega] = IERP(D);
    t3 = toc
    tic
    [output_Atilde2,output_Omega2] = IERP_speedtest(D);
    t4 = toc
    data3diff = find(abs(output_Atilde2-output_Atilde)>0.00001)
    data4diff = find(abs(output_Omega2-output_Omega)>0.00001)
    

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
