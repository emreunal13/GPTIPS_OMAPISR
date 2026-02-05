function gp = fix_best_after_run(gp)
%FIX_BEST_AFTER_RUN  One-shot repair of gp.results.best from fitness arrays.

% Safety checks
if ~isfield(gp,'results') || ~isfield(gp.results,'best') ...
        || ~isfield(gp.results.best,'index')
    warning('fix_best_after_run: no gp.results.best.index found.');
    return;
end

ci = gp.results.best.index;   % best individual index

% Sanity: within bounds?
if ci < 1 || ci > numel(gp.fitness.values)
    warning('fix_best_after_run: index out of range.'); 
    return;
end

% Pull the authoritative values from the fitness bookkeeping:
fit_i   = gp.fitness.values(ci);
comp_i  = gp.fitness.complexity(ci);
theta_i = gp.fitness.returnvalues{ci};
cvals_i = gp.fitness.constraint_values{ci};

% Overwrite gp.results.best.*
gp.results.best.fitness          = fit_i;
gp.results.best.returnvalues     = theta_i;
gp.results.best.constraint_values = cvals_i;

if isfield(gp.fitness,'complexityMeasure') && gp.fitness.complexityMeasure
    gp.results.best.complexity = comp_i;
else
    gp.results.best.nodecount  = comp_i;
end

% And keep the convenience alias in sync
gp.best = gp.results.best;
end
