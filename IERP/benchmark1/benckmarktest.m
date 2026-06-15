clear,clc
% TEST 1: When the given demand is a realizable effective resistance
% matrix, the shortest path distance matrix 


N = 50
p=0.4
noise_amplitude = 0.5;

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
D(D < 0) = Input_Omega(D < 0);


[output_Atilde, output_Omega] = GL_leastsquares(D);

[output_Atilde1,output_Omega1] = IERP(D);


% % test when D is realizable
% A_result = output_Atilde - A_input;
% find(abs(A_result)>0.01)
% omega_result = D- output_Omega ;
% find(abs(omega_result)>0.01)

% test for others 

Norm_output = sum(sum((abs(D - output_Omega))./(D+(D==0))))/(N*(N-1))

Norm_output1 = sum(sum((abs(D - output_Omega1))./(D+(D==0))))/(N*(N-1))




