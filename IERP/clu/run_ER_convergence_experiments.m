clear;
clc;

% Run one ER convergence experiment for each network size and save results.
ierpRoot = 'C:\Users\86748\Documents\MATLAB\flowbetweenness\IERP';
cluDir = fullfile(ierpRoot, 'clu');
dataDir = 'D:\data\flow betweenness\IERP\convergence';

addpath(ierpRoot);
addpath(cluDir);

if ~exist(dataDir, 'dir')
    mkdir(dataDir);
end

N_values = [20, 50, 100, 200];
p = 0.5;
baseSeed = 20260615;

for idx = 1:numel(N_values)
    N = N_values(idx);
    rngSeed = baseSeed + N;
    rng(rngSeed, 'twister');

    fprintf('Running convergence experiment for N = %d...\n', N);

    A_input = GenerateERfast(N, p, 10);
    while ~network_isconnected(A_input)
        A_input = GenerateERfast(N, p, 10);
    end

    A_input(A_input ~= 0) = 1 ./ A_input(A_input ~= 0);
    Input_Omega = EffectiveResistance(A_input);
    D = Input_Omega;

    runStart = tic;
    [output_Atilde, output_Omega, diff_history] = ...
        IERP_forconvergencey(D);
    elapsed_seconds = toc(runStart);

    diff_history_accepted = diff_history;
    if numel(diff_history) >= 2 && ...
            abs(diff_history(end)) > abs(diff_history(end - 1))
        diff_history_accepted = diff_history(1:end - 1);
    end

    resultFile = fullfile(dataDir, ...
        sprintf('ER_convergence_N%d.mat', N));

    save(resultFile, ...
        'N', 'p', 'rngSeed', 'A_input', 'Input_Omega', 'D', ...
        'output_Atilde', 'output_Omega', 'diff_history', ...
        'diff_history_accepted', 'elapsed_seconds', '-v7.3');

    fprintf(['Saved %s (%d accepted pruning steps, ', ...
        '%.2f seconds).\n'], resultFile, ...
        numel(diff_history_accepted) - 1, elapsed_seconds);
end

fprintf('All convergence experiments completed.\n');

