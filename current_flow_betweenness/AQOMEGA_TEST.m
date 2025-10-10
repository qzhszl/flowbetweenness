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
Omega = EffectiveResistance(A_input)

% Fielder block matrix
A_2 = Fielder_block_matrix(Omega)





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