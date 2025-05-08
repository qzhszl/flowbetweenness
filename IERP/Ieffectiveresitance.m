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
% 7. repeat 2-6 until \Omega>D : Return lastly removed link

% 1. generate a graph
% (a) tree:
% N = 10;
% T = generate_a_tree(N,1,10);
% A = full(adjacency(T,"weighted"));
% subplot(2,2,1)
% plot(T,'EdgeLabel',T.Edges.Weight,'NodeColor',[0.8500 0.3250 0.0980], ...
% 'EdgeAlpha',0.5,'LineWidth',1,'MarkerSize',7,'EdgeLabelColor',[0 0.4470 0.7410],'NodeFontSize',10);

% 1. generate a graph
% (a) tree:
N = 10;
% T = generate_a_tree(N,1,10);
% A = full(adjacency(T,"weighted"));
A_input = [0	0	0	0	0	0	0	0	3	0
    0	0	0	0	0	0	6	10	0	0
    0	0	0	0	2	0	0	2	0	0
    0	0	0	0	0	0	0	4	0	0
    0	0	2	0	0	10	0	0	4	0
    0	0	0	0	10	0	0	0	0	0
    0	6	0	0	0	0	0	0	0	0
    0	10	2	4	0	0	0	0	0	0
    3	0	0	0	4	0	0	0	0	3
    0	0	0	0	0	0	0	0	3	0];

T = graph(A_input);
% subplot(2,2,1)
% h1 = plot(T,'EdgeLabel',T.Edges.Weight,'NodeColor',[0.8500 0.3250 0.0980], ...
% 'EdgeAlpha',0.5,'LineWidth',1,'MarkerSize',7,'EdgeLabelColor',[0 0.4470 0.7410],'NodeFontSize',10);
% x = h1.XData;
% y = h1.YData;
% D = distances(T);
Aforomega = A_input;
Aforomega(Aforomega ~= 0) = 1 ./ Aforomega(Aforomega ~= 0);
Omega = EffectiveResistance(Aforomega);
% Generate a demand matrix
D = Omega

Input_Omega = D;
W_tilde  = Input_Omega
W_tilde(W_tilde ~= 0) = 1 ./ W_tilde(W_tilde ~= 0);
Omega_new = EffectiveResistance(W_tilde)            % Compute the effective resistance Omega
D-Omega_new
sum(sum(D-Omega_new))
val = 1;
A = (W_tilde > 0);
while(nnz(Omega_new > D) == 0 && val > 0)                     % Remove links one by one until we exceed the constraints                            
    % method 1 R = A.*((Omega_new+eye(N)).^-1 - W_tilde)*(D-Omega_new) 
    R = A.*((Omega_new+eye(N)).^-1 - W_tilde).*(D-Omega_new)     % Compute R
    [val,~] = max(max(R));                              % Identify the maximum element
    [row,col] = find(R == val);                        % Identify the link
    A(row(1),col(1)) = 0; A(col(1),row(1)) = 0;        % Remove the link
    W_tilde = A.*D;
    W_tilde(W_tilde ~= 0) = 1 ./ W_tilde(W_tilde ~= 0);      % Compute W tilde
    Omega_new = EffectiveResistance(W_tilde);               % Update the shortest path weight matrix
    D-Omega_new                             
    sum(sum(D-Omega_new))
end

if(row(1) ~= col(1))
    A(row(1),col(1)) = 1; A(col(1),row(1)) = 1;         % Return lastly removed link
    W = A.*D
%     W = 2*1/sum(sum(A.*D)).*A.*D;                       % Update the weighted adjacency matrix
end
final_result = EffectiveResistance(W);
A_OLR = A;
% Store the results
% 1. The number of links added in the graph
L_add_OLR = 0.5*(nnz(A_OLR)-nnz(A_input))

% 2. The number of links in the obtained graph

L_DOR = 0.5*nnz(A_OLR)

% 3. The number of common links between two graphs
L_comm_OLR = nnz(A_input.*A_OLR)/nnz(A_input); 

Norm_OLR = sum(sum((D - final_result)./(D+(D==0))))/(N*(N-1))

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
