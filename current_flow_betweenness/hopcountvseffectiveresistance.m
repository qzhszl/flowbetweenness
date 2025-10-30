clear,clc

% 1. Generate an effecitive resistance matrix 
% 2. regard this as an demand matrix
% 3. Use DOR and to solve it.

% 1. generate a graph
% (a) tree:
N = 10;
T = generate_a_tree(N,1,10);
A = full(adjacency(T));

% A = GenerateERfast(N,0.5,0);

A = [0	0	1	1	0	1	0	1	1	0
0	0	1	0	0	0	1	1	0	1
1	1	0	0	0	1	0	0	1	1
1	0	0	0	0	1	1	1	0	1
0	0	0	0	0	0	0	1	0	0
1	0	1	1	0	0	0	0	1	0
0	1	0	1	0	0	0	1	0	1
1	1	0	1	1	0	1	0	0	0
1	0	1	0	0	1	0	0	0	0
0	1	1	1	0	0	1	0	0	0]

subplot(2,2,1)

G = graph(A);
H = distances(G);
plot(G,'EdgeLabel',G.Edges.Weight,'NodeColor',[0.8500 0.3250 0.0980], ...
'EdgeAlpha',0.5,'LineWidth',1,'MarkerSize',7,'EdgeLabelColor',[0 0.4470 0.7410],'NodeFontSize',10);



% A sometimes is not invertable
det(H)
rank(H)



% % (b) tree and a link:
% N = 10;
% A = [0	0	0	0	0	0	0	0	3	6
%     0	0	0	0	0	0	6	10	0	0
%     0	0	0	0	2	0	0	2	0	0
%     0	0	0	0	0	0	0	4	0	0
%     0	0	2	0	0	10	0	0	4	0
%     0	0	0	0	10	0	0	0	0	0
%     0	6	0	0	0	0	0	0	0	0
%     0	10	2	4	0	0	0	0	0	0
%     3	0	0	0	4	0	0	0	0	3
%     6	0	0	0	0	0	0	0	3	0];
% T = graph(A);
% subplot(2,2,1)
% plot(T,'EdgeLabel',T.Edges.Weight,'NodeColor',[0.8500 0.3250 0.0980], ...
% 'EdgeAlpha',0.5,'LineWidth',1,'MarkerSize',7,'EdgeLabelColor',[0 0.4470 0.7410],'NodeFontSize',10);

% (c) ER:
% N = 10;
% A = GenerateERfastinteger(N,0.5,1)
% T = graph(A);
% subplot(2,2,1)
% plot(T,'EdgeLabel',T.Edges.Weight,'NodeColor',[0.8500 0.3250 0.0980], ...
% 'EdgeAlpha',0.5,'LineWidth',1,'MarkerSize',7,'EdgeLabelColor',[0 0.4470 0.7410],'NodeFontSize',10);

% 1. Compute the effective resistance matrix
% Aforomega = A;
% Aforomega(Aforomega ~= 0) = 1 ./ Aforomega(Aforomega ~= 0);
% Omega = EffectiveResistance(Aforomega);
% D = Omega;





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


function A = GenerateERfastinteger(n,p,weighted)
    A = rand(n,n) < p;
    A = triu(A,1);
    if weighted == 1
        linkweight_matrix = randi(10,n,n);
        A = A.*linkweight_matrix;
    end
    
    A = A + A';
end