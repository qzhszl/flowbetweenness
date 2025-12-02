clear, clc
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

A_input = [];
folder_name = "D:\\data\\flow betweenness\\IERP\\realnetwork\\";
filename = "soc-wiki-Vote.mtx";
filename = "soc-dolphins.mtx";
filename = "soc-karate.mtx";
totalname = folder_name+filename;
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
NODENUM = size(A,1)
Linknum = nnz(A)/2



% [L_add_output,L_ouput,L_comm_output,Norm_output] = experiment_on_emperical_network(A_input);
% result = [L_add_output,L_ouput,L_comm_output,Norm_output]
% filename = sprintf("D:\\data\\flow betweenness\\IERP\\realnetwork\\IERPsoc-wiki-Vote.txt");
% writematrix(result,filename)

function [L_add_output,L_ouput,L_comm_output_ratio,Norm_output] = experiment_on_emperical_network(A_input)
    % % ensure the input link weight change from 0,1
    % Input_Omega = EffectiveResitance_withinverseA(A_input);

    Input_Omega = EffectiveResistance(A_input);
    
    % Generate a demand matrix: the effective resistance matrix of the
    % input network
    D = Input_Omega;
    N = size(D,1);
    % tic
    [output_Atilde,output_Omega] = IERP(D);
    % t3 = toc
    % tic
    % [output_Atilde2,output_Omega2] = IERP_speedtest(D);
    % t4 = toc
    % data3diff = find(abs(output_Atilde2-output_Atilde)>0.00001)
    % data4diff = find(abs(output_Omega2-output_Omega)>0.00001)
    

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



function Omega = EffectiveResitance_withinverseA(A)
% the resitance is the inverse of the link weight
    Aforomega = A;
    Aforomega(Aforomega ~= 0) = 1 ./ Aforomega(Aforomega ~= 0);
    Omega = EffectiveResistance(Aforomega);
end

