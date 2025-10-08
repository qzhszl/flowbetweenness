function Simulations_on_ER_network_exponentiallinkweight_server(N,inputpara)
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


rng(inputpara)


x = log(N)/N;
y = ceil(x * 1e4) / 1e4;  % round 4 decimal
p_start = y;
    
simutimes = 50;

result = zeros(simutimes,4);
p_vec = linspace(p_start, 1, 15);
% 前两个点
p1 = p_vec(1);
p2 = p_vec(2);  
% 在 p1 和 p2 之间插入两个点
extra_points = linspace(p1, p2, 4);  % 生成4个点
extra_points = extra_points(2:3);    % 去掉第一个和最后一个（原本已有）
% 合并
p_vec = [p_vec(1), extra_points, p_vec(2:end)];
p_vec = round(p_vec,4);

for p= p_vec
    for simu_time = 1:simutimes
        current_output_s = sprintf("N%d, p%.4f: %d/%d",N,p,simu_time,simutimes);
        disp(current_output_s)
        % 1. generate a graph
        % _________________________________________________________________________
        % (b) ER:
        A_input = GenerateERfastexp(N,p,0.5);
        % check connectivity
        connect_flag = network_isconnected(A_input);
        while ~connect_flag
            A_input = GenerateERfastexp(N,p,0.5);
            % check connectivity
            connect_flag = network_isconnected(A_input);
        end

        % 2. run simulations
        [L_add_output,L_ouput,L_comm_output,Norm_output] = experiment_on_ER_exp(A_input);
        
        result(simu_time,:) = [L_add_output,L_ouput,L_comm_output,Norm_output];
    end
    filename = sprintf("IERP_N%dERp%.4f_weight_exp_simu%d.txt",N,p,inputpara);
    writematrix(result,filename)
end

end


