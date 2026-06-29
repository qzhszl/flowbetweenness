clear,clc
for i = 1:10
    N = 20;
    p= 0.4;
    noise_amplitude = 0.6;
    
    
    A_input = GenerateERfast(N,p,10);
    % check connectivity
    connect_flag = network_isconnected(A_input);
    while ~connect_flag
        A_input = GenerateERfast(N,p,10);
        % check connectivity
        connect_flag = network_isconnected(A_input);
    end
    
    % 2. run simulations
    A_input(A_input ~= 0) = 1 ./ A_input(A_input ~= 0);
    Input_Omega = EffectiveResistance(A_input);
        
    % Generate a demand matrix: the effective resistance matrix of the
    % input network
    D = Input_Omega;


    N = size(D,1);
    % Generate random noise matrix
    Z = randn(N, N);
    % Make Z symmetric
    Z = (Z + Z') / 2;
    % Set diagonal to zero
    Z(1:N+1:end) = 0;
    
    % Generate perturbed demand matrix
    D = D .* (1 + noise_amplitude * Z);
    D = D  + noise_amplitude * Z;

    D(D < 0.01) = Input_Omega(D < 0.01);
    
    % tic
    [output_Atilde,output_Omega] = IERP(D);
    % t3 = toc
    % tic
    [output_Atilde2,output_Omega2] = GL_leastsquares(D);
    % t4 = toc
    % data3diff = find(abs(output_Atilde2-output_Atilde)>0.00001)
    % data4diff = find(abs(output_Omega2-output_Omega)>0.00001)
    % output_Atilde2
    % D1 = (abs(D - output_Omega))./(D+(D==0))
    % D2 = (abs(D - output_Omega2))./(D+(D==0))
    
    % Store the results
    % 1. The number of links added in the graph
    L_add_output = 0.5*(nnz(output_Atilde)-nnz(A_input));
    % 2. The number of links in the obtained graph
    L_ouput = 0.5*nnz(output_Atilde);  
    % 3. The number of common links between two graphs
    L_comm_output_ratio = nnz(A_input.*output_Atilde)/nnz(output_Atilde); 
    % 4. The norm of common links between two graphs
    Norm_output = sum(sum((abs(D - output_Omega))./(D+(D==0))))/(N*(N-1))
    
    
    % Store the results FOR the benchmark
    L_add_output2 = 0.5*(nnz(output_Atilde2)-nnz(A_input));
    % 2. The number of links in the obtained graph
    L_ouput2 = 0.5*nnz(output_Atilde2);  
    % 3. The number of common links between two graphs
    L_comm_output_ratio2 = nnz(A_input.*output_Atilde2)/nnz(output_Atilde2); 
    % 4. The norm of common links between two graphs
    Norm_output2 = sum(sum((abs(D - output_Omega2))./(D+(D==0))))/(N*(N-1))
    if Norm_output2>10
        D
        output_Atilde2
        D2 = (abs(D - output_Omega2))./(D+(D==0))
        % A1, A2 是两个加权邻接矩阵
        G0 = graph(A_input,'upper');
        G2 = graph(output_Atilde2,'upper');
        G1 = graph(output_Atilde,'upper');

        % 用 G1 生成固定节点坐标
        figure;
        p0 = plot(G1,'Layout','force');
        x = p0.XData;
        y = p0.YData;
        close;

        figure;
        subplot(1,3,1)
        p0 = plot(G0,'XData',x,'YData',y);
        title('Graph 0');
        p0.EdgeLabel = compose('%.2f', G0.Edges.Weight);


        subplot(1,3,2)
        p1 = plot(G1,'XData',x,'YData',y);
        title('Graph 1');
        p1.EdgeLabel = compose('%.2f', G1.Edges.Weight);

        lw1 = 1 + 5 * G1.Edges.Weight / max(G1.Edges.Weight);
        p1.LineWidth = lw1;

        subplot(1,3,3)
        p2 = plot(G2,'XData',x,'YData',y);
        title('Graph 2');
        p2.EdgeLabel = compose('%.2f', G2.Edges.Weight);

        lw2 = 1 + 5 * G2.Edges.Weight / max(G2.Edges.Weight);
        p2.LineWidth = lw2;
        break
    end
    


end



% % A1, A2 是两个加权邻接矩阵
% G0 = graph(A_input,'upper');
% G2 = graph(output_Atilde2,'upper');
% G1 = graph(output_Atilde,'upper');
% 
% % 用 G1 生成固定节点坐标
% figure;
% p0 = plot(G1,'Layout','force');
% x = p0.XData;
% y = p0.YData;
% close;

% figure;
% subplot(1,3,1)
% p0 = plot(G0,'XData',x,'YData',y);
% title('Graph 0');
% p0.EdgeLabel = compose('%.2f', G0.Edges.Weight);
% 
% 
% subplot(1,3,2)
% p1 = plot(G1,'XData',x,'YData',y);
% title('Graph 1');
% p1.EdgeLabel = compose('%.2f', G1.Edges.Weight);
% 
% lw1 = 1 + 5 * G1.Edges.Weight / max(G1.Edges.Weight);
% p1.LineWidth = lw1;
% 
% subplot(1,3,3)
% p2 = plot(G2,'XData',x,'YData',y);
% title('Graph 2');
% p2.EdgeLabel = compose('%.2f', G2.Edges.Weight);
% 
% lw2 = 1 + 5 * G2.Edges.Weight / max(G2.Edges.Weight);
% p2.LineWidth = lw2;