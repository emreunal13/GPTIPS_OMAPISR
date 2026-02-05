function gp = updatestats(gp)
%UPDATESTATS Updates run statistics and stores the BEST individual's physics.

% Sanitize fitness values (tiny imag → real, big imag → penalty)
vals = gp.fitness.values;

if ~isreal(vals)
    imagv     = imag(vals);
    smallImag = abs(imagv) < 1e-12;      % treat as numerical noise
    vals(smallImag) = real(vals(smallImag));

    bigImag = ~smallImag & (imagv ~= 0); % genuinely complex
    vals(bigImag) = 1e6;                 % hard penalty, but *only* for those
end

bad = isnan(vals) | isinf(vals);
vals(bad) = 1e6;

gp.fitness.values = vals;   % commit sanitized values

% Find best in this generation using these cleaned values
[minval, minind] = min(vals);

% Initialize global best on first generation
if gp.state.count == 1
    if ~isfield(gp,'results'); gp.results = []; end
    gp.results.best.fitness = Inf;
end

% If this generation’s best is better, overwrite gp.results.best
if minval <= gp.results.best.fitness
    gp.results.best.fitness         = minval;
    gp.results.best.individual      = gp.pop{minind};
    gp.results.best.eval_individual = tree2evalstr(gp.pop{minind}, gp);
    gp.results.best.index           = minind;

    % theta
    if isfield(gp.fitness,'returnvalues') && numel(gp.fitness.returnvalues) >= minind
        gp.results.best.returnvalues = gp.fitness.returnvalues{minind};
    else
        gp.results.best.returnvalues = [];
    end

    % constraints
    if isfield(gp.fitness,'constraint_values') && numel(gp.fitness.constraint_values) >= minind
        gp.results.best.constraint_values = gp.fitness.constraint_values{minind};
    else
        gp.results.best.constraint_values = [];
    end

    % complexity / nodecount
    if isfield(gp.fitness,'complexityMeasure') && gp.fitness.complexityMeasure
        gp.results.best.complexity = gp.fitness.complexity(minind);
    else
        gp.results.best.nodecount  = getnumnodes(gp.pop{minind});
    end

    gp.best = gp.results.best;
end

%mean fitness & history
gp.state.meanfitness = mean(gp.fitness.values);
gp.results.history.bestfitness(gp.state.count) = gp.results.best.fitness;
gp.results.history.meanfitness(gp.state.count) = gp.state.meanfitness;
