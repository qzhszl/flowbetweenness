clear,clc
% In this file, we would like to figure out the confused system introduced
% by link weight = 1/r, where r is the reistance of link


% weigthed adj or link weight matrix W:
A_input = [0,0.5,2
        0.5 0 1
        2 1 0]

% resistance matrix R:
R = A_input;
R(R ~= 0) = 1 ./ R(R ~= 0)

% Laplacian
DEGREE = sum(A_input)
Qtilde = diag(DEGREE)-A_input


% The effective resistance
% if we want the effective resistance scale two times, then the original
% adj should be half of the original one

Omega = EffectiveResistance(A_input)
% Omega2 = EffectiveResistance(0.5*A_input)

% Fielder block matrix
A_2 = Fielder_block_matrix(Omega)


[output_Atilde,output_Omega] = IERP_test(Omega)

[output_Atilde1,output_Omega1] = IERP(Omega)




function [output_Atilde,output_Omega] = IERP_test(D)
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
    W(W ~= 0) = 1 ./ W(W ~= 0);
    output_Atilde = W;
  
    output_Omega = EffectiveResistance(output_Atilde);
    alpha = alpha_l1_global_para(output_Omega,D);
    output_Omega = alpha*output_Omega;
    output_Atilde = 1/alpha*output_Atilde;
end






function Omega = EffectiveResistance(Atilde)
m=size(Atilde,1);

DEGREE = sum(Atilde);
Qtilde = diag(DEGREE)-Atilde;

PQ = pinv(Qtilde);

xita = diag(PQ);
u = ones(m,1);

Omega = u*xita.'+xita*u.'-2*PQ;
end



function A = Fielder_block_matrix(Omega)
    m = size(Omega,1);
    u = ones(m,1);
    p = 1/(u.'*inv(Omega)*u)*inv(Omega)*u;
    Q_tilde = 2*(u.'*inv(Omega)*u)*p*p.'-2*inv(Omega);
    A = diag(diag(Q_tilde)) - Q_tilde;
    A = round(A, 10);
    % A(A ~= 0) = 1 ./ A(A ~= 0);
end