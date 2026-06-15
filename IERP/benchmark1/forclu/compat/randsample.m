function sample = randsample(population, sampleSize)
%RANDSAMPLE Minimal sampling-without-replacement compatibility function.
%   This two-input implementation covers the call made by the preserved
%   NeurIPS 2018 effResGD.m when Statistics Toolbox is unavailable.

    if isscalar(population)
        population = 1:population;
    else
        population = population(:).';
    end

    populationSize = numel(population);
    validateattributes(sampleSize, {'numeric'}, ...
        {'scalar', 'integer', 'nonnegative', '<=', populationSize});

    indices = randperm(populationSize, sampleSize);
    sample = population(indices).';
end
