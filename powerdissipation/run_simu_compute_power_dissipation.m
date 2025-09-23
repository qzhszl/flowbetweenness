clear,clc
% A = [0 1 0 1
%     1 0 2 0
%     0 2 0 0
%     1 0 0 0];
% G = graph(A);
% plot(G,'EdgeLabel',G.Edges.Weight,'NodeColor',[0.8500 0.3250 0.0980], ...
% 'EdgeAlpha',0.5,'LineWidth',1,'MarkerSize',7,'EdgeLabelColor',[0 0.4470 0.7410],'NodeFontSize',10);
% 
% 
% A_path = zeros(size(A));
% nz = A~=0;        % 非零位置
% A_path(nz) = 1./A(nz);
% G_path = graph(A_path);
% 
% 
% [power_givenst_SP,EdgeEnergy_SP] = compute_path_power_dissipation(G_path,1,3)
% [power_givenst_flow,EdgeEnergy_flow] = compute_flownetwork_power_dissipation(G,1,3)

% EdgeEnergy_SP = addvars(G.Edges, EdgeEnergy_SP, 'NewVariableNames',"Energy_SP")


for N = 200
    % for p=[0.15,0.28,0.39,0.66,0.88]
    for p=[0.88]
        compute_power_dissipation_eachlink(N,p,100)
    end
end


