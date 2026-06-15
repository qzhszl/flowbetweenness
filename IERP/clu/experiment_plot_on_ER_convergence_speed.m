clear,clc
N=100
p = 0.5
% input a D
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

% See the difference of the D

[output_Atilde, output_Omega, diff_history] = IERP_forconvergencey(D);
diff_history_accepted = diff_history;

if length(diff_history) >= 2 && diff_history(end) > diff_history(end-1)
    diff_history_accepted = diff_history(1:end-1);
end


figure;
plot(0:length(diff_history_accepted)-1, diff_history_accepted, 'LineWidth', 1.8);
xlabel('Accepted pruning step');
ylabel('\epsilon');
grid on;

