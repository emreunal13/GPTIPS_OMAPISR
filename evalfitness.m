function gp = evalfitness(gp)
%EVALFITNESS  Serial fitness evaluation.
%
%   GP = EVALFITNESS(GP) evaluates the fitnesses of individuals stored
%   in the GP structure and updates various other fields of GP accordingly.
%
%   If gp.runcontrol.parallel.enable && gp.runcontrol.parallel.ok is true,
%   this function immediately delegates to EVALFITNESS_PAR and returns.

% If parallel mode is enabled, delegate and exit
if gp.runcontrol.parallel.enable && gp.runcontrol.parallel.ok
    gp = evalfitness_par(gp);
    return;
end

popSize = gp.runcontrol.pop_size;
fitfun  = gp.fitness.fitfun;

% Preallocate
complexities   = nan(popSize,1);
fitvals        = nan(popSize,1);
returnvals     = cell(popSize,1);
constraintvals = cell(popSize,1);
constrScores   = nan(popSize,1);
paramvals      = cell(popSize,1);   % theta + ERCs/exponents (if provided)

useCache = isfield(gp.runcontrol,'usecache') && gp.runcontrol.usecache;
hasCache = isfield(gp.fitness,'cache');

for i = 1:popSize

    gp.state.current_individual = i;

    % ---------------------------------------------------------------------
    % 1. Cache branch (if enabled)
    % ---------------------------------------------------------------------
    if useCache && hasCache && gp.fitness.cache.isKey(i)

        cache = gp.fitness.cache(i);

        % Basic values
        if isfield(cache,'complexity')
            complexities(i) = double(cache.complexity);
        end
        if isfield(cache,'value')
            fitvals(i)      = double(cache.value);
        end

        % Theta
        if isfield(cache,'returnvalues')
            returnvals{i,1} = cache.returnvalues;
        else
            returnvals{i,1} = [];
        end

        % Param vector (theta + ERCs), if cached
        if isfield(cache,'param_values')
            paramvals{i,1} = cache.param_values;
        else
            paramvals{i,1} = [];
        end

        % Constraint struct + Cscore
        if isfield(cache,'constraint_values') && ~isempty(cache.constraint_values)
            constraintvals{i,1} = cache.constraint_values;
            if isstruct(cache.constraint_values) && ...
                    isfield(cache.constraint_values,'Cscore')
                constrScores(i) = cache.constraint_values.Cscore;
            else
                constrScores(i) = 1.0;
            end
        else
            constraintvals{i,1} = [];
            constrScores(i)     = 1.0;
        end

    else
        
        % 2. Fresh evaluation branch
        

        % preprocess cell array of string expressions
        evalstr = tree2evalstr(gp.pop{i},gp);

        % complexity (either node count or expressional complexity)
        if isfield(gp.fitness,'complexityMeasure') && gp.fitness.complexityMeasure
            complexities(i) = double(getcomplexity(gp.pop{i}));
        else
            complexities(i) = double(getnumnodes(gp.pop{i}));
        end

        % Call the fitness function (which also performs NM + Lamarckian update)
        [fitness_i,gp] = feval(fitfun,evalstr,gp);

        % sanitise fitness: keep complex allowed, only kill NaN/Inf
        if ~(isnumeric(fitness_i) && isscalar(fitness_i) && all(isfinite(fitness_i)))
            fitness_i = 1e6;
        end
        fitvals(i) = double(fitness_i);

        % pull out theta that the fitfun stored
        if isfield(gp.fitness,'returnvalues') && numel(gp.fitness.returnvalues) >= i
            returnvals{i} = gp.fitness.returnvalues{i};
        else
            returnvals{i} = [];
        end

        % pull out full param vector (theta + ERCs/exponents), if stored
        if isfield(gp.fitness,'param_values') && numel(gp.fitness.param_values) >= i
            paramvals{i} = gp.fitness.param_values{i};
        else
            paramvals{i} = [];
        end

        % pull out constraint values and Cscore
        rvC    = [];
        cscore = 1.0;
        if isfield(gp.fitness,'constraint_values') && ...
           numel(gp.fitness.constraint_values) >= i
            rvC = gp.fitness.constraint_values{i};
            if isstruct(rvC) && isfield(rvC,'Cscore')
                cscore = rvC.Cscore;
            end
        end
        constraintvals{i} = rvC;
        constrScores(i)   = cscore;

        % 3. Update cache if enabled 
        if useCache && hasCache
            cachedData.complexity        = complexities(i);
            cachedData.value             = fitvals(i);
            cachedData.returnvalues      = returnvals{i};
            cachedData.constraint_values = constraintvals{i};
            cachedData.param_values      = paramvals{i};
            gp.fitness.cache(i)          = cachedData;
        end
    end
end

% Post-pass: sanitise constraint scores (only NaN/Inf -> 1)
bad = ~isfinite(constrScores);
constrScores(bad) = 1.0;

% Commit to main GP structure
gp.fitness.values            = fitvals(:);
gp.fitness.complexity        = complexities(:);
gp.fitness.returnvalues      = returnvals;
gp.fitness.constraint_values = constraintvals;
gp.fitness.param_values      = paramvals;

% MULTI-OBJECTIVE MATRIX: [RMSE, Complexity, ConstraintScore]
gp.fitness.multiobj = [gp.fitness.values, ...
                       gp.fitness.complexity, ...
                       constrScores(:)];

end
