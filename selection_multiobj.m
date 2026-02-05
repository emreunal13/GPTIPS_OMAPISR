function idx = selection_multiobj(gp)
% SELECTION_MULTIOBJ Tournament selection using 3 objectives:
%   columns: [fitness (RMSE), complexity, constraint_score]
% Lower is better for all objectives.


try
    popSize = gp.runcontrol.pop_size;
catch
    idx = selection(gp); return;
end

if ~isfield(gp.fitness,'multiobj') || isempty(gp.fitness.multiobj) ...
        || size(gp.fitness.multiobj,1) ~= popSize
    idx = selection(gp);
    return;
end

P = gp.fitness.multiobj;  % popSize x 3

% tournament size: prefer gp.selection.tournament_size if available, else 2
if isfield(gp,'selection') && isfield(gp.selection,'tournament_size') ...
        && ~isempty(gp.selection.tournament_size)
    tsize = max(2, round(gp.selection.tournament_size));
else
    tsize = 2;
end

% pick tsize competitors
competitors = randi(popSize, tsize, 1);

% Precompute nondominated ranks for entire pop (cheap for typical pop sizes)
[fronts, rank] = nondominated_sort(P);

% pick best among competitors: smallest rank; if tie -> crowding distance in that front
best = competitors(1);
for k = 2:tsize
    cand = competitors(k);
    if rank(cand) < rank(best)
        best = cand;
    elseif rank(cand) == rank(best)
        fidx = rank(cand);
        members = fronts{fidx};
        localP = P(members,:);
        cd = crowding_distance(localP);
        % find local indices
        loc_map = zeros(max(members),1); loc_map(members) = 1:numel(members);
        local_cand_idx = loc_map(cand);
        local_best_idx = loc_map(best);
        % if either has Inf, prefer Inf (edge)
        cd_c = cd(local_cand_idx);
        cd_b = cd(local_best_idx);
        if cd_c > cd_b
            best = cand;
        end
    end
end

idx = best;
end
