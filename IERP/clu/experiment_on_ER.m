function [L_add_output,L_ouput,L_comm_output_ratio,Norm_output] = experiment_on_ER(A_input)
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

