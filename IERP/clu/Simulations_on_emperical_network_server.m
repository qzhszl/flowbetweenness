function Simulations_on_emperical_network_server()
    % this .m is to verify Inverse effecitve resitance problem:
    % Given an effective resitance matrix D, we are planning to obtain a network
    % whose effective resitance is similar with given demand matrix
    
    
    % the main idea is :1. start with a complete graph whose link weight = demand
    % 2. R = \Omega - W_tilde
    % 3. find i,j where reaches tshe largest R
    % 4. remove that link
    % 5. renormalize the link
    % 6. update the \Omega
    % 7. repeat 2-6 until Diff between Omega and the given demand is minimum: Return lastly removed link
    
    % A_input = [];
    % folder_name = "/home/zqiu1/flowbetweeness/IERP/";
    filename = "soc-wiki-Vote.mtx";
    totalname = filename;
    fid = fopen(totalname);
    if fid == -1
        error('无法打开文件：%s，请检查路径和文件名。', totalname);
    end
    % 跳过头部注释行
    line = fgetl(fid);
    while ischar(line) && (startsWith(line, '%') || isempty(line))
        line = fgetl(fid);
    end
    
    % 第一行非注释行：行数、列数、非零元数量
    dims = sscanf(line, '%d %d %d');
    m = dims(1);
    n = dims(2);
    nz = dims(3);
    
    % 读取后续 (i,j,val)
    data = fscanf(fid, '%d %d %f', [2 nz])';
    fclose(fid);
    
    % 构造稀疏矩阵
    A = sparse(data(:,1), data(:,2), 1, m, n);
    A = A + A.';
    A = A - diag(diag(A));
    A_input = full(A);
    
    [L_add_output,L_ouput,L_comm_output,Norm_output] = experiment_on_emperical_network(A_input);
    result = [L_add_output,L_ouput,L_comm_output,Norm_output]
    filename = sprintf("IERPsoc-wiki-Vote.txt");
    writematrix(result,filename)
end

