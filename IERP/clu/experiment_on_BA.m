function [L_add_output,L_ouput,L_comm_output_ratio,Norm_output] = experiment_on_BA(A_input)
    Input_Omega = EffectiveResitance_withinverseA(A_input);
    % Generate a demand matrix: the effective resistance matrix of the
    % input network
    D = Input_Omega;
    N = size(D,1);
    [output_Atilde,output_Omega] = IERP(D);  

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
