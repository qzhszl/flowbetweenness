% Most recent one: Include link weighT change

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

function [output_Atilde,output_Omega] = IERP(D)
    N = size(D,1);
    Input_Omega = D;
    W_tilde  = Input_Omega;
    W_tilde(W_tilde ~= 0) = 1 ./ W_tilde(W_tilde ~= 0);
    Omega_new = EffectiveResistance(W_tilde);            % Compute the effective resistance Omega
    diff_change = sum(sum((abs(D - Omega_new))./(D+(D==0))))/(N*(N-1));

    % val = 1;
    flag = 1;
    A = (W_tilde > 0);
    % Gnow = graph(W_tilde,"upper");
    while(flag==1)                     % Remove links one by one until we exceed the constraints                            
        previous_change = diff_change;
        % method 1 R = A.*((Omega_new+eye(N)).^-1 - W_tilde)*(D-Omega_new) 
        R = A.*((Omega_new+eye(N)).^-1 - W_tilde).*(D-Omega_new);     % Compute R
        [val,~] = max(max(R));                              % Identify the maximum element
        [row,col] = find(R == val);                        % Identify the link
        A(row(1),col(1)) = 0; A(col(1),row(1)) = 0;        % Remove the link
        W_tilde = A.*D;
        W_tilde(W_tilde ~= 0) = 1 ./ W_tilde(W_tilde ~= 0);      % Compute W tilde
        % Gnow = graph(W_tilde,"upper");
        Omega_new = EffectiveResistance(W_tilde);               % Update the shortest path weight matrix
        
        diff_change = sum(sum((abs(D - Omega_new))./(D+(D==0))))/(N*(N-1));

        if abs(diff_change) > abs(previous_change)
            flag=0;
        end    
    end
    
    if(row(1) ~= col(1))
        A(row(1),col(1)) = 1; A(col(1),row(1)) = 1;         % Return lastly removed link
        W = A.*D;
    %     W = 2*1/sum(sum(A.*D)).*A.*D;                       % Update the weighted adjacency matrix
    end
    output_Atilde = W;
    output_Omega = EffectiveResitance_withinverseA(W);
    alpha = alpha_l1_global_para(output_Omega,D);
    output_Omega = alpha*output_Omega;
    output_Atilde = alpha*output_Atilde;
end
