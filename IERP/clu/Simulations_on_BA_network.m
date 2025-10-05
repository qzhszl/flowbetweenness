function Simulations_on_BA_network()
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
    
    filefolder_name_input = "/home/zqiu1/flowbetweeness/IERP/";
    N = 1000;
    m=2;
    input_BAname = sprintf('BAnetworkN%dm%d.txt',N,m);
    simutimes = 1;
    data = readmatrix(input_BAname);  % 每行: u v w
    
     % 提取三列
    s = data(:,1)+1;    % 起点
    t = data(:,2)+1;    % 终点
    w = data(:,3);    % 权重
    
    % 构造无向图（如果是有向图，就用 digraph）
    G = graph(s, t, w);
    A_input = full(adjacency(G,"weighted"));
    % 2. run simulations
    tic
    result = zeros(simutimes,4);
    [L_add_output,L_ouput,L_comm_output,Norm_output] = experiment_on_BA(A_input);
    result(1,:) = [L_add_output,L_ouput,L_comm_output,Norm_output];
    filename = sprintf("IERP_N%dBAm%d_weight01.txt",N,m);
    writematrix(result,filename)
end
