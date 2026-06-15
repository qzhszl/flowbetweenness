function [output_Atilde, output_Omega] = GL_leastsquares(D)
%GL_LEASTSQUARES Learn an undirected graph from target effective resistances.
%   [A, OMEGA] = GL_LEASTSQUARES(D) minimizes the least-squares objective
%
%       sum_{i<j, D(i,j)>0} (OMEGA(i,j) - D(i,j))^2
%
%   over nonnegative undirected edge weights. D must be a finite,
%   symmetric, nonnegative matrix with a zero diagonal. As in the authors'
%   reference implementation, a zero off-diagonal entry is unconstrained.
%
%   The implementation preserves and calls the original MATLAB code from:
%   Hoskins, Musco, Musco, and Tsourakakis, NeurIPS 2018.
%
%   For complete exact resistance data, the paper's closed-form
%   initialization is already a global least-squares minimizer, so the
%   routine returns it without unnecessary gradient iterations. Otherwise,
%   it runs the paper's projected gradient or coordinate-descent method.

    D = localValidateDemandMatrix(D);
    n = size(D, 1);

    if n == 1
        output_Atilde = 0;
        output_Omega = 0;
        return;
    end

    thisDir = fileparts(mfilename('fullpath'));
    upstreamDir = fullfile(thisDir, 'upstream_graph_similarity_learning');
    utilsDir = fullfile(upstreamDir, 'utils');
    compatibilityDir = fullfile(thisDir, 'compat');
    localRequireUpstream(upstreamDir, utilsDir, compatibilityDir);

    originalPath = path;
    pathCleanup = onCleanup(@() path(originalPath)); %#ok<NASGU>
    addpath(upstreamDir, utilsDir, '-begin');
    if exist('randsample', 'file') == 0
        addpath(compatibilityDir, '-begin');
    end

    [r, u, v] = localMatrixToVector(D);
    constrained = find(r > 0);
    if isempty(constrained)
        error('GL_leastsquares:NoConstraints', ...
            'D must contain at least one positive off-diagonal constraint.');
    end

    warningState = warning;
    warningCleanup = onCleanup(@() warning(warningState)); %#ok<NASGU>
    warning('off', 'MATLAB:singularMatrix');
    warning('off', 'MATLAB:nearlySingularMatrix');

    w0 = localPaperInitialization(D, r, constrained);
    initialError = localNormalizedObjective(w0, r, constrained, n);

    % The exact all-pairs case is already the global minimizer of Problem 2.
    if initialError <= 1e-12
        finalWeights = w0;
    else
        maxIterations = 3000;
        m = numel(r);
        numConstraints = numel(constrained);

        fullGradientBytes = 8 * double(m) * double(numConstraints);
        incidenceBytes = 8 * double(m) * double(n);
        useSmallSolver = fullGradientBytes <= 384e6 && incidenceBytes <= 128e6;

        referenceLaplacian = localWeightsToLaplacian(w0, n);
        temporaryFolder = tempname;
        mkdir(temporaryFolder);
        originalFolder = pwd;
        folderCleanup = onCleanup( ...
            @() localCleanupTemporaryRun(originalFolder, temporaryFolder)); %#ok<NASGU>
        cd(temporaryFolder);

        if useSmallSolver
            [weightHistory, ~, ~, ~] = effResGDSmall( ...
                sparse(r), referenceLaplacian, 0, w0, ...
                maxIterations, 0, 'GDLS');
        else
            batchByMemory = max(1, floor(256e6 / (8 * double(numConstraints))));
            batchSize = min([m, 5000, batchByMemory]);
            [weightHistory, ~, ~, ~] = effResGD( ...
                sparse(r), u, v, referenceLaplacian, 0, w0, ...
                maxIterations, 0, batchSize, 'GDLS');
        end

        finalWeights = full(weightHistory(end, :)).';
        if any(~isfinite(finalWeights))
            error('GL_leastsquares:NumericalFailure', ...
                'The least-squares solver produced non-finite edge weights.');
        end
        finalWeights = max(finalWeights, 0);
    end

    output_Atilde = localVectorToMatrix(finalWeights, n);
    output_Atilde = (output_Atilde + output_Atilde.') / 2;
    output_Atilde(1:n+1:end) = 0;
    output_Atilde(output_Atilde < 0) = 0;

    if ~localIsConnected(output_Atilde)
        output_Atilde = localConnectComponents(output_Atilde, D);
    end

    output_Omega = localEffectiveResistance(output_Atilde);
end

function D = localValidateDemandMatrix(D)
    validateattributes(D, {'numeric'}, ...
        {'2d', 'square', 'real', 'finite', 'nonempty'}, ...
        mfilename, 'D', 1);

    D = double(D);
    scale = max(1, max(abs(D), [], 'all'));
    tolerance = 1e-12 * scale;

    if any(abs(diag(D)) > tolerance)
        error('GL_leastsquares:NonzeroDiagonal', ...
            'D must have a zero diagonal.');
    end
    if any(D < -tolerance, 'all')
        error('GL_leastsquares:NegativeDemand', ...
            'D must be nonnegative.');
    end

    symmetryError = norm(D - D.', 'fro') / max(1, norm(D, 'fro'));
    if symmetryError > 1e-10
        error('GL_leastsquares:AsymmetricDemand', ...
            ['Effective resistance is symmetric. D must be symmetric ', ...
             'to define the undirected problem in the paper.']);
    end

    D = max((D + D.') / 2, 0);
    D(1:size(D, 1)+1:end) = 0;
end

function localRequireUpstream(upstreamDir, utilsDir, compatibilityDir)
    required = {
        fullfile(upstreamDir, 'effResGD.m')
        fullfile(upstreamDir, 'effResGDSmall.m')
        fullfile(upstreamDir, 'exactRecover.m')
        fullfile(utilsDir, 'pair2index.m')
        fullfile(utilsDir, 'w2A.m')
        fullfile(utilsDir, 'w2L.m')
        fullfile(compatibilityDir, 'randsample.m')
    };

    for i = 1:numel(required)
        if ~isfile(required{i})
            error('GL_leastsquares:MissingUpstreamCode', ...
                'Missing preserved upstream file: %s', required{i});
        end
    end
end

function w0 = localPaperInitialization(D, r, constrained)
    n = size(D, 1);
    m = numel(r);
    targets = r(constrained);
    targetScale = median(targets);

    completedD = localShortestPathCompletion(D, targetScale);
    completedR = localMatrixToVector(completedD);

    lambdaCandidates = targetScale * [0, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2];
    bestError = inf;
    w0 = ones(m, 1) / max(targetScale, eps);

    for lambda = lambdaCandidates
        try
            candidate = full(exactRecover(completedR, lambda));
        catch
            continue;
        end

        if any(~isfinite(candidate))
            continue;
        end

        % This is the projection used by the paper: remove negative edges.
        candidate = max(candidate, 0);
        candidate = localEnsureConnectedWeights(candidate, n, targetScale);
        candidate = localScaleWeights(candidate, r, constrained, n);
        candidateError = localNormalizedObjective( ...
            candidate, r, constrained, n);

        if candidateError < bestError
            bestError = candidateError;
            w0 = candidate;
        end
    end

    if ~isfinite(bestError)
        w0 = localEnsureConnectedWeights(w0, n, targetScale);
        w0 = localScaleWeights(w0, r, constrained, n);
    end
end

function completedD = localShortestPathCompletion(D, targetScale)
    n = size(D, 1);
    completedD = D;
    missing = (completedD == 0) & ~eye(n);
    if ~any(missing, 'all')
        return;
    end

    measurementGraph = graph(D, 'upper');
    shortest = distances(measurementGraph);
    finitePositive = shortest(isfinite(shortest) & shortest > 0);

    if isempty(finitePositive)
        fallbackDistance = max(targetScale, eps) * n;
    else
        fallbackDistance = max(2 * max(finitePositive), targetScale * n);
    end
    shortest(~isfinite(shortest)) = fallbackDistance;
    completedD(missing) = shortest(missing);
    completedD = (completedD + completedD.') / 2;
    completedD(1:n+1:end) = 0;
end

function w = localEnsureConnectedWeights(w, n, targetScale)
    A = localVectorToMatrix(w, n);
    components = localComponents(A);
    if max(components) == 1
        return;
    end

    positiveWeights = w(w > 0);
    if isempty(positiveWeights)
        baseWeight = 1 / max(targetScale, eps);
    else
        baseWeight = median(positiveWeights);
    end
    bridgeWeight = max(baseWeight * 1e-10, realmin('double'));

    componentCount = max(components);
    representatives = zeros(componentCount, 1);
    for c = 1:componentCount
        representatives(c) = find(components == c, 1, 'first');
    end
    for c = 1:componentCount-1
        i = representatives(c);
        j = representatives(c + 1);
        if i > j
            [i, j] = deal(j, i);
        end
        w(localPairIndex(n, i, j)) = bridgeWeight;
    end
end

function w = localScaleWeights(w, targetR, constrained, n)
    A = localVectorToMatrix(w, n);
    currentOmega = localEffectiveResistance(A);
    currentR = localMatrixToVector(currentOmega);
    ratios = currentR(constrained) ./ targetR(constrained);
    ratios = ratios(isfinite(ratios) & ratios > 0);
    if ~isempty(ratios)
        w = w * median(ratios);
    end
end

function value = localNormalizedObjective(w, targetR, constrained, n)
    A = localVectorToMatrix(w, n);
    currentOmega = localEffectiveResistance(A);
    currentR = localMatrixToVector(currentOmega);
    residual = currentR(constrained) - targetR(constrained);
    denominator = sum(targetR(constrained).^2);
    value = sum(residual.^2) / max(denominator, eps);
end

function Omega = localEffectiveResistance(A)
    n = size(A, 1);
    L = diag(sum(A, 2)) - A;
    Lplus = pinv((L + L.') / 2);
    diagonal = diag(Lplus);
    Omega = diagonal + diagonal.' - 2 * Lplus;
    Omega = max((Omega + Omega.') / 2, 0);
    Omega(1:n+1:end) = 0;
end

function [vector, u, v] = localMatrixToVector(M)
    n = size(M, 1);
    m = n * (n - 1) / 2;
    vector = zeros(m, 1);
    u = zeros(m, 1);
    v = zeros(m, 1);

    index = 0;
    for i = 1:n
        for j = i+1:n
            index = index + 1;
            vector(index) = M(i, j);
            u(index) = i;
            v(index) = j;
        end
    end
end

function A = localVectorToMatrix(w, n)
    A = zeros(n);
    index = 0;
    for i = 1:n
        for j = i+1:n
            index = index + 1;
            A(i, j) = w(index);
            A(j, i) = w(index);
        end
    end
end

function L = localWeightsToLaplacian(w, n)
    A = localVectorToMatrix(w, n);
    L = diag(sum(A, 2)) - A;
end

function index = localPairIndex(n, i, j)
    index = (i - 1) * (n - i / 2) + (j - i);
end

function connected = localIsConnected(A)
    components = localComponents(A);
    connected = max(components) == 1;
end

function components = localComponents(A)
    n = size(A, 1);
    adjacency = A > 0;
    components = zeros(n, 1);
    componentId = 0;

    for startNode = 1:n
        if components(startNode) ~= 0
            continue;
        end
        componentId = componentId + 1;
        queue = startNode;
        components(startNode) = componentId;
        head = 1;

        while head <= numel(queue)
            node = queue(head);
            head = head + 1;
            neighbors = find(adjacency(node, :) & components.' == 0);
            components(neighbors) = componentId;
            queue = [queue, neighbors]; %#ok<AGROW>
        end
    end
end

function A = localConnectComponents(A, D)
    n = size(A, 1);
    components = localComponents(A);
    componentCount = max(components);
    positiveWeights = A(A > 0);

    if isempty(positiveWeights)
        positiveTargets = D(D > 0);
        bridgeWeight = 1 / max(median(positiveTargets), eps);
    else
        bridgeWeight = median(positiveWeights) * 1e-10;
    end

    representatives = zeros(componentCount, 1);
    for c = 1:componentCount
        representatives(c) = find(components == c, 1, 'first');
    end
    for c = 1:componentCount-1
        i = representatives(c);
        j = representatives(c + 1);
        A(i, j) = bridgeWeight;
        A(j, i) = bridgeWeight;
    end
    A(1:n+1:end) = 0;
end

function localCleanupTemporaryRun(originalFolder, temporaryFolder)
    if isfolder(originalFolder)
        cd(originalFolder);
    end
    if isfolder(temporaryFolder)
        rmdir(temporaryFolder, 's');
    end
end


