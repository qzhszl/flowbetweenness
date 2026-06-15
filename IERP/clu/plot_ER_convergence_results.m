function fig = plot_ER_convergence_results(dataDir)
%PLOT_ER_CONVERGENCE_RESULTS Plot accepted convergence histories by N.
%   The horizontal coordinate is the cumulative number of accepted pruned
%   edges divided by the number of possible undirected edges, nchoosek(N,2).

if nargin < 1 || isempty(dataDir)
    dataDir = 'D:\data\flow betweenness\IERP\convergence';
end

N_values = [20, 50, 100, 200];
colors = ["#D08082", "#C89FBF", "#62ABC7", "#7A7DB1", "#6FB494", "#D9B382"];

fig = figure('Color', 'w', 'Position', [100, 100, 450, 300]);
ax = axes(fig);
hold(ax, 'on');

for idx = 1:numel(N_values)
    N_expected = N_values(idx);
    resultFile = fullfile(dataDir, ...
        sprintf('ER_convergence_N%d.mat', N_expected));

    if ~isfile(resultFile)
        error('Missing result file: %s', resultFile);
    end

    result = load(resultFile, 'N', 'diff_history', ...
        'diff_history_accepted');

    if isfield(result, 'N') && result.N ~= N_expected
        error('Unexpected N in %s: expected %d, found %d.', ...
            resultFile, N_expected, result.N);
    end

    if isfield(result, 'diff_history_accepted')
        acceptedHistory = result.diff_history_accepted;
    else
        acceptedHistory = result.diff_history;
        if numel(acceptedHistory) >= 2 && ...
                abs(acceptedHistory(end)) > abs(acceptedHistory(end - 1))
            acceptedHistory = acceptedHistory(1:end - 1);
        end
    end

    acceptedHistory = acceptedHistory(:);
    prunedEdgeCount = (0:numel(acceptedHistory) - 1).';
    normalizedPrunedEdgeCount = ...
        prunedEdgeCount / nchoosek(N_expected, 2);

    plot(ax, normalizedPrunedEdgeCount, acceptedHistory, ...
        'LineWidth', 2.0, ...
        'Color', colors(idx), ...
        'DisplayName', sprintf('N = %d', N_expected));
end



legend(ax, 'Location', 'best', 'Box', 'on',FontSize=14);
% grid(ax, 'on');
box(ax, 'on');
set(ax, 'FontName', 'Times New Roman', 'FontSize', 13, ...
    'LineWidth', 1.0, 'TickLabelInterpreter', 'latex');
xlabel(ax, ...
    '$L_p/{N \choose 2}$', ...
    'Interpreter', 'latex','FontSize', 16);
ylabel(ax, ...
    '$\epsilon$', ...
    'Interpreter', 'latex','FontSize', 16);

exportgraphics(fig, fullfile(dataDir, ...
    'ER_convergence_normalized_pruning.pdf'), ...
    'ContentType', 'vector');
exportgraphics(fig, fullfile(dataDir, ...
    'ER_convergence_normalized_pruning.png'), ...
    'Resolution', 300);
end
