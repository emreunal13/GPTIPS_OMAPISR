function ordered = pareto_order(gp)
% PARETO_ORDER Deterministic ordering of population indices by Pareto front.
%   ordered = pareto_order(gp)
%
%   Requires gp.fitness.multiobj (popSize x 3) ideally. If that's missing,
%   builds a fallback multiobj as [fitness, complexity, zeros].
%
%   Ordering: iterate front = 1 .. end; inside each front sort by:
%      1) fitness (ascending = better)
%      2) complexity (ascending)
%      3) constraint score (ascending)
%
%   Returns ordered: vector of population indices (1..popSize)

popSize = gp.runcontrol.pop_size;

% Obtain objective matrix P (popSize x 3), lower is better.
if isfield(gp.fitness,'multiobj') && ~isempty(gp.fitness.multiobj) ...
        && size(gp.fitness.multiobj,1) == popSize
    P = gp.fitness.multiobj;
else
    % fallback: use available fields and warn
    warning('pareto_order: gp.fitness.multiobj missing or mismatched. Building fallback multiobj using fitness & complexity.');
    fvals = gp.fitness.values(:);
    if isfield(gp.fitness,'complexity')
        comp = gp.fitness.complexity(:);
    else
        comp = zeros(size(fvals));
    end
    % third objective = 0 (no constraint score available)
    P = [fvals, comp, zeros(size(fvals))];
end

% nondominated sort -> get front membership and rank
[fronts, rank] = nondominated_sort(P);

% Build ordered list: for each front, sort its members by the tie-break keys
ordered = zeros(popSize,1);
pos = 1;
for fi = 1:numel(fronts)
    members = fronts{fi};            % global indices in this front
    if isempty(members)
        continue;
    end
    % extract objective columns for members
    subP = P(members, :);           % [nMembers x 3]
    % sort rows by (fitness, complexity, constraint) ascending
    % stable sort: build composite key
    [~, idx] = sortrows(subP, [1,3, 2]); % ascending on each column
    sortedMembers = members(idx);
    n = numel(sortedMembers);
    ordered(pos:pos+n-1) = sortedMembers(:);
    pos = pos + n;
end

% truncated if duplicates/empties
ordered = ordered(1:pos-1);

end
