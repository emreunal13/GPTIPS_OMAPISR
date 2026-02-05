function gp = evalfitness_par(gp)
%EVALFITNESS_PAR Parallel fitness evaluation (Lamarckian included)

popSize = gp.runcontrol.pop_size;
fitfun  = gp.fitness.fitfun;

% Preallocate
complexities   = nan(popSize,1);
fitvals        = nan(popSize,1);
returnvals     = cell(popSize,1);
constraintvals = cell(popSize,1);
constrScores   = nan(popSize,1);
paramvals      = cell(popSize,1);   % <- full param vectors (theta + ERCs)
newPop         = gp.pop;            % <- will hold Lamarckian-updated pop

parfor idx = 1:popSize
    tempgp = gp;
    tempgp.state.current_individual = idx;

    try
        % 1. Complexity
        if isfield(tempgp.fitness,'complexityMeasure') && tempgp.fitness.complexityMeasure
            complexities(idx) = double(getcomplexity(tempgp.pop{idx}));
        else
            complexities(idx) = double(getnumnodes(tempgp.pop{idx}));
        end

        % 2. Evaluate Fitness & optimize constants (calls regressmulti_fitfun)
        evalstr = tree2evalstr(tempgp.pop{idx}, tempgp);

        % Force re-computation to ensure Nelder–Mead runs
        tempgp.state.force_compute_theta = true;

        [fitness_i, tempgp] = feval(fitfun, evalstr, tempgp);

        % Sanitize Fitness
        if ~(isnumeric(fitness_i) && isscalar(fitness_i) && isfinite(fitness_i))
            fitness_i = 1e6;
        end
        fitvals(idx) = double(fitness_i);

        % 3. Extract optimized weights (theta)
        if isfield(tempgp.fitness,'returnvalues') && ...
           numel(tempgp.fitness.returnvalues) >= idx && ...
           ~isempty(tempgp.fitness.returnvalues{idx})
            returnvals{idx} = tempgp.fitness.returnvalues{idx};
        else
            returnvals{idx} = [];
        end

        % 4. Extract full param vector (theta + ERCs), if present
        if isfield(tempgp.fitness,'param_values') && ...
           numel(tempgp.fitness.param_values) >= idx && ...
           ~isempty(tempgp.fitness.param_values{idx})
            paramvals{idx} = tempgp.fitness.param_values{idx};
        else
            paramvals{idx} = [];
        end

        % 5. Extract constraint scores
        rvC    = [];
        cscore = 1.0; % Default worst case
        if isfield(tempgp.fitness,'constraint_values') && ...
           numel(tempgp.fitness.constraint_values) >= idx && ...
           ~isempty(tempgp.fitness.constraint_values{idx})
            rvC = tempgp.fitness.constraint_values{idx};
            if isstruct(rvC) && isfield(rvC,'Cscore')
                cscore = rvC.Cscore;
            end
        end
        constraintvals{idx} = rvC;
        constrScores(idx)   = cscore;

        % 6. Lamarckian: grab updated individual coming out of fitfun
        %    (regressmulti_fitfun updated tempgp.pop{idx} via
        %     update_individual_constants_in_pop)
        newPop{idx} = tempgp.pop{idx};

    catch
        % Fallback for crash
        complexities(idx)   = 1e6;
        fitvals(idx)        = 1e6;
        returnvals{idx}     = [];
        constraintvals{idx} = [];
        constrScores(idx)   = 1.0;
        paramvals{idx}      = [];
        % Keep original gp.pop{idx} in newPop (no Lamarckian change)
    end
end

bad = ~isfinite(constrScores);
constrScores(bad) = 1.0;


% Commit to main structure
gp.fitness.values            = double(fitvals(:));
gp.fitness.complexity        = double(complexities(:));
gp.fitness.returnvalues      = returnvals;
gp.fitness.constraint_values = constraintvals;
gp.fitness.param_values      = paramvals;
gp.pop                       = newPop;   % <- Lamarckian population update

% Multi-Objective MATRIX: [RMSE, Complexity, ConstraintScore]
gp.fitness.multiobj          = [gp.fitness.values, ...
                                gp.fitness.complexity, ...
                                double(constrScores(:))];

gp.state.current_individual = popSize;
end
