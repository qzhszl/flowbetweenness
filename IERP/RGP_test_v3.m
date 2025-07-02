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
% 7. repeat 2-6 until minimum Diff between i and j: Return lastly removed link


% This.m is for understand for large dense network, if it is possible to
% reducing the time complexity： yes it is


% 1. generate a graph
% (a) tree:
% N = 10;
% T = generate_a_tree(N,1,10);
% A = full(adjacency(T,"weighted"));
% subplot(2,2,1)
% plot(T,'EdgeLabel',T.Edges.Weight,'NodeColor',[0.8500 0.3250 0.0980], ...
% 'EdgeAlpha',0.5,'LineWidth',1,'MarkerSize',7,'EdgeLabelColor',[0 0.4470 0.7410],'NodeFontSize',10);

% 1. generate a graph
% _________________________________________________________________________
% (a) tree:
% N = 20;
% T = generate_a_tree(N,1,10);
% A_input = full(adjacency(T,"weighted"));

% A_input = [0	0	0	0	0	0	0	0	3	6
%     0	0	0	0	0	0	6	10	0	0
%     0	0	0	0	2	0	0	2	0	0
%     0	0	0	0	0	0	0	4	0	0
%     0	0	2	0	0	10	0	0	4	0
%     0	0	0	0	10	0	0	0	0	0
%     0	6	0	0	0	0	0	0	0	0
%     0	10	2	4	0	0	0	0	0	0
%     3	0	0	0	4	0	0	0	0	3
%     6	0	0	0	0	0	0	0	3	0];

% (b) ER network
% A_input = [0	5	0	4	0	5
% 5	0	0	5	0	2
% 0	0	0	9	0	0
% 4	5	9	0	5	8
% 0	0	0	5	0	0
% 5	2	0	8	0	0];
% A_input = [0	7	4	9	7	0	0	0	7	2	7	1	8	8	1	8	10	0	6	0
% 7	0	0	2	9	1	0	0	6	0	0	8	2	6	5	7	3	0	0	0
% 4	0	0	8	0	0	2	9	10	0	8	7	0	4	9	0	3	0	4	9
% 9	2	8	0	10	4	4	0	0	6	0	1	0	0	2	0	0	0	8	0
% 7	9	0	10	0	0	1	6	0	0	1	0	5	0	0	10	6	1	0	0
% 0	1	0	4	0	0	4	0	0	5	6	9	5	0	0	0	3	10	0	2
% 0	0	2	4	1	4	0	8	6	9	0	0	3	2	10	4	0	1	0	1
% 0	0	9	0	6	0	8	0	0	0	0	5	4	0	0	0	1	2	7	0
% 7	6	10	0	0	0	6	0	0	0	0	0	9	9	0	0	0	6	0	0
% 2	0	0	6	0	5	9	0	0	0	0	0	9	0	4	5	6	0	10	0
% 7	0	8	0	1	6	0	0	0	0	0	8	3	5	9	0	1	4	0	10
% 1	8	7	1	0	9	0	5	0	0	8	0	7	4	0	4	6	3	5	0
% 8	2	0	0	5	5	3	4	9	9	3	7	0	0	8	10	0	2	0	0
% 8	6	4	0	0	0	2	0	9	0	5	4	0	0	0	0	0	10	5	8
% 1	5	9	2	0	0	10	0	0	4	9	0	8	0	0	0	5	3	4	5
% 8	7	0	0	10	0	4	0	0	5	0	4	10	0	0	0	0	4	0	6
% 10	3	3	0	6	3	0	1	0	6	1	6	0	0	5	0	0	0	5	8
% 0	0	0	0	1	10	1	2	6	0	4	3	2	10	3	4	0	0	8	0
% 6	0	4	8	0	0	0	7	0	10	0	5	0	5	4	0	5	8	0	0
% 0	0	9	0	0	2	1	0	0	0	10	0	0	8	5	6	8	0	0	0];

% N = 1000;
% A_input = GenerateERfast(N,0.5,10)
% T = graph(A_input);


% BA network
N = 1000;
filename = 'D:\\data\\flow betweenness\\BAnetworks\\BAnetworkN1000m2.txt';
data = readmatrix(filename);

% 拆分数据
sources = data(:,1);
sources = sources+1;
targets = data(:,2);
targets = targets+1;
weights = data(:,3);

% 构建带权无向图
G = graph(sources, targets, weights);
A_input  = full(adjacency(G,"weighted"));

tic
% subplot(1,2,1)
% h1 = plot(T,'EdgeLabel',T.Edges.Weight,'NodeColor',[0.8500 0.3250 0.0980], ...
% 'EdgeAlpha',0.5,'LineWidth',1,'MarkerSize',7,'EdgeLabelColor',[0 0.4470 0.7410],'NodeFontSize',10);
% x = h1.XData;
% y = h1.YData;
% D = distances(T);
Aforomega = A_input;
Aforomega(Aforomega ~= 0) = 1 ./ Aforomega(Aforomega ~= 0);
Omega = EffectiveResistance(Aforomega);
% Generate a demand matrix
D = Omega;

Input_Omega = D;
W_tilde  = Input_Omega;
W_tilde(W_tilde ~= 0) = 1 ./ W_tilde(W_tilde ~= 0);
Omega_new = EffectiveResistance(W_tilde);            % Compute the effective resistance Omega
diff_change = sum(sum((abs(D - Omega_new))./(D+(D==0))))/(N*(N-1));

val = 1;
flag = 1;
A = (W_tilde > 0);
Gnow = graph(W_tilde,"upper");
% resistormatrix = W_tilde;
printcount = 0

while(flag==1 && val > 0 && all(conncomp(Gnow) == 1))                     % Remove links one by one until we exceed the constraints                            
    previous_change = diff_change;
    % method 1 
    R = A.*((Omega_new+eye(N)).^-1 - W_tilde).*(D-Omega_new); 
    % R = A.*((Omega_new+eye(N)).^-1 - W_tilde).*(D-Omega_new).*resistormatrix; 
    % R = A.*(D-Omega_new).*resistormatrix;     % Compute R
    
    ratio = 2 * printcount / (N * (N - 1));
    if ratio <0.95 & N>500
        R_flat = R(:);
        % 找出 R 中最大的10个元素及其索引
        [~, idx] = maxk(R_flat, 200);
        
        % 将线性索引转换为行列索引
        [row_idx, col_idx] = ind2sub(size(R), idx);
        
        % 将 A 中对应位置及其转置位置设为 0
        for k = 1:length(row_idx)
            i = row_idx(k);
            j = col_idx(k);
            A(i,j) = 0;
            A(j,i) = 0;  % 设为转置位置也为0
        end
        W_tilde = A.*D;
        % resistormatrix = W_tilde;
        W_tilde(W_tilde ~= 0) = 1 ./ W_tilde(W_tilde ~= 0);      % Compute W tilde
        Gnow = graph(W_tilde,"upper");
        Omega_new = EffectiveResistance(W_tilde);               % Update the shortest path weight matrix
                                
        diff_change = sum(sum((abs(D - Omega_new))./(D+(D==0))))/(N*(N-1));
        
        fprintf('diff_change: %.8f\n', diff_change);
        if abs(diff_change) > abs(previous_change)
            flag=0;
        end
        printcount = printcount+ 100;
        
    else
        [val,~] = max(max(R));                              % Identify the maximum element
        [row,col] = find(R == val);                        % Identify the link
        A(row(1),col(1)) = 0; A(col(1),row(1)) = 0;        % Remove the link
        W_tilde = A.*D;
        % resistormatrix = W_tilde;
        W_tilde(W_tilde ~= 0) = 1 ./ W_tilde(W_tilde ~= 0);      % Compute W tilde
        Gnow = graph(W_tilde,"upper");
        Omega_new = EffectiveResistance(W_tilde);               % Update the shortest path weight matrix
                               
        diff_change = sum(sum((abs(D - Omega_new))./(D+(D==0))))/(N*(N-1));

        fprintf('diff_change: %.8f\n', diff_change);
        if abs(diff_change) > abs(previous_change)
            flag=0;
        end
        printcount = printcount+ 1;
    end

    if mod(printcount, 1000) == 0
        ratio = 2 * printcount / (N * (N - 1));
        fprintf('Progress: %.2f%%\n', ratio * 100);
    end
end

if(row(1) ~= col(1))
    A(row(1),col(1)) = 1; A(col(1),row(1)) = 1;         % Return lastly removed link
    W = A.*D;
%     W = 2*1/sum(sum(A.*D)).*A.*D;                       % Update the weighted adjacency matrix
end

final_result = EffectiveResitance_withinverseA(W);
toc
OmegaDiff = round(D-Omega_new,10);
diff_change = sum(sum(OmegaDiff));


A_OLR = A;
% Store the results
% 1. The number of links added in the graph
L_add_OLR = 0.5*(nnz(A_OLR)-nnz(A_input))

% 2. The number of links in the obtained graph

L_DOR = 0.5*nnz(A_OLR)

% 3. The number of common links between two graphs
L_comm_OLR = nnz(A_input.*A_OLR)/nnz(A_input); 
% 4. The norm of common links between two graphs
Norm_OLR = sum(sum((abs(D - final_result))./(D+(D==0))))/(N*(N-1))

G_OLR = graph(W,"upper");
% subplot(1,2,2)
% h2 = plot(G_OLR,'EdgeLabel',G_OLR.Edges.Weight,'XData', x, 'YData', y,'NodeColor',[0.8500 0.3250 0.0980], ...
% 'EdgeAlpha',0.5,'LineWidth',1,'MarkerSize',7,'EdgeLabelColor',[0 0.4470 0.7410],'NodeFontSize',10);




% Aforomega = A;
% Aforomega(Aforomega ~= 0) = 1 ./ Aforomega(Aforomega ~= 0);
% Omega = EffectiveResistance(Aforomega);
% Omega - D
% 
% A_DOR = DOR(D);
% G_DOR = graph(A_DOR,'upper');
% subplot(2,2,2)
% h2 = plot(G_DOR,'EdgeLabel',G_DOR.Edges.Weight,'XData', x, 'YData', y,'NodeColor',[0.8500 0.3250 0.0980], ...
% 'EdgeAlpha',0.5,'LineWidth',1,'MarkerSize',7,'EdgeLabelColor',[0 0.4470 0.7410],'NodeFontSize',10);
% 
% 
% A_P = ISPP_tree(D);
% G_P = graph(A_P,'upper');
% subplot(2,2,3)
% plot(G_P,'EdgeLabel',G_P.Edges.Weight,'XData', x, 'YData', y,'NodeColor',[0.8500 0.3250 0.0980], ...
% 'EdgeAlpha',0.5,'LineWidth',1,'MarkerSize',7,'EdgeLabelColor',[0 0.4470 0.7410],'NodeFontSize',10);




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